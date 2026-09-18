from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
from typing import Any


RUNBOOK_SCHEMA = "agentic-rag-production-promotion-runbook-v1"
GAP_SCHEMA = "agentic-rag-production-gap-summary-v1"


def main() -> int:
    parser = argparse.ArgumentParser(description="Render a production promotion runbook from production audit evidence.")
    parser.add_argument("--audit-json", type=Path, help="Path to production_audit.json.")
    parser.add_argument("--write-result", type=Path, help="Path to production_audit.py write-result JSON.")
    parser.add_argument("--format", choices=["markdown", "json"], default="markdown", help="Output format.")
    parser.add_argument("--output", type=Path, help="Optional output path. Defaults to stdout.")
    args = parser.parse_args()

    try:
        audit = load_audit(args)
        runbook = build_runbook(audit)
    except Exception as exc:
        print(f"production_promotion_runbook failed: {exc}", file=sys.stderr)
        return 2

    text = json.dumps(runbook, ensure_ascii=False, indent=2, sort_keys=True) if args.format == "json" else render_markdown(runbook)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text)
    return 0 if runbook["verification"]["status"] != "fail" else 1


def load_audit(args: argparse.Namespace) -> dict[str, Any]:
    if bool(args.audit_json) == bool(args.write_result):
        raise ValueError("Provide exactly one of --audit-json or --write-result.")
    if args.audit_json:
        return read_json(args.audit_json)
    result = read_json(args.write_result)
    path = result.get("json_path")
    if not path:
        raise ValueError("write-result JSON does not contain json_path.")
    return read_json(Path(path))


def read_json(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    for encoding in ("utf-8", "utf-8-sig", "utf-16"):
        try:
            payload = json.loads(raw.decode(encoding))
            if not isinstance(payload, dict):
                raise ValueError(f"{path} does not contain a JSON object.")
            return payload
        except (UnicodeDecodeError, json.JSONDecodeError):
            continue
    raise ValueError(f"Could not decode JSON file: {path}")


def build_runbook(audit: dict[str, Any]) -> dict[str, Any]:
    gap_summary = audit.get("production_gap_summary") if isinstance(audit.get("production_gap_summary"), dict) else {}
    gaps = gap_summary.get("gaps") if isinstance(gap_summary.get("gaps"), list) else []
    verification = verify_gap_summary(gap_summary)
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for gap in gaps:
        if isinstance(gap, dict):
            groups[str(gap.get("resolution_scope") or "unknown")].append(_runbook_item(gap))
    grouped = [
        {
            "resolution_scope": scope,
            "count": len(items),
            "items": sorted(items, key=lambda item: (item["priority"], item["area"])),
        }
        for scope, items in sorted(groups.items())
    ]
    proof_checklist = _proof_checklist(grouped)
    source_proof_checklist = gap_summary.get("proof_checklist") if isinstance(gap_summary.get("proof_checklist"), list) else []
    payload = {
        "runbook_schema": RUNBOOK_SCHEMA,
        "source_audit_hash": audit.get("audit_hash"),
        "source_generated_at": audit.get("generated_at"),
        "status": gap_summary.get("status") or verification.get("status") or "warning",
        "production_ready": ((audit.get("summary") or {}).get("production_ready") is True),
        "gap_count": len(gaps),
        "blocking_count": safe_int(gap_summary.get("blocking_count")),
        "promotion_count": safe_int(gap_summary.get("promotion_count")),
        "resolution_scope_counts": gap_summary.get("resolution_scope_counts") or {},
        "resolution_mode_counts": gap_summary.get("resolution_mode_counts") or {},
        "source_proof_checklist_count": safe_int(gap_summary.get("proof_checklist_count")),
        "source_proof_checklist_hash": gap_summary.get("proof_checklist_hash"),
        "proof_checklist_source_match": proof_checklist == source_proof_checklist,
        "proof_checklist_count": len(proof_checklist),
        "proof_checklist_hash": hash_payload(proof_checklist),
        "proof_checklist": proof_checklist,
        "groups": grouped,
        "verification": verification,
    }
    payload["runbook_hash"] = hash_payload({key: value for key, value in payload.items() if key != "runbook_hash"})
    return payload


def _runbook_item(gap: dict[str, Any]) -> dict[str, Any]:
    proof = _proof_evidence(gap)
    return {
        "area": gap.get("area"),
        "priority": gap.get("priority"),
        "status": gap.get("status"),
        "resolution_mode": gap.get("resolution_mode"),
        "proof_hint": gap.get("proof_hint"),
        "proof_artifact": gap.get("proof_artifact") or proof["artifact"],
        "proof_command": gap.get("proof_command") or proof["command"],
        "action": gap.get("action"),
        "evidence_key": gap.get("evidence_key"),
        "check_detail_hash": gap.get("check_detail_hash"),
        "readiness_detail_hash": gap.get("readiness_detail_hash"),
        "gap_hash": gap.get("gap_hash"),
    }


def _proof_checklist(groups: list[dict[str, Any]]) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    for group in groups:
        for item in group.get("items") or []:
            if not isinstance(item, dict):
                continue
            items.append(
                proof_checklist_item(
                    {
                        "area": item.get("area"),
                        "priority": item.get("priority"),
                        "resolution_scope": group.get("resolution_scope"),
                        "resolution_mode": item.get("resolution_mode"),
                        "proof_artifact": item.get("proof_artifact"),
                        "proof_command": item.get("proof_command"),
                        "gap_hash": item.get("gap_hash"),
                    }
                )
            )
    return sorted(items, key=lambda item: (str(item.get("priority") or ""), str(item.get("area") or "")))


def _gap_proof_checklist(gaps: list[Any]) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    for gap in gaps:
        if not isinstance(gap, dict):
            continue
        items.append(
            proof_checklist_item(
                {
                    "area": gap.get("area"),
                    "priority": gap.get("priority"),
                    "resolution_scope": gap.get("resolution_scope"),
                    "resolution_mode": gap.get("resolution_mode"),
                    "proof_artifact": gap.get("proof_artifact"),
                    "proof_command": gap.get("proof_command"),
                    "gap_hash": gap.get("gap_hash"),
                }
            )
        )
    return sorted(items, key=lambda item: (str(item.get("priority") or ""), str(item.get("area") or "")))


def _proof_evidence(gap: dict[str, Any]) -> dict[str, str]:
    area = str(gap.get("area") or "unknown")
    evidence_key = str(gap.get("evidence_key") or area)
    resolution_mode = str(gap.get("resolution_mode") or "")
    area_commands = {
        "deployment_readiness": "Invoke-RestMethod http://127.0.0.1:8000/api/v1/deployment/readiness",
        "rag_embedding_backend": "Invoke-RestMethod http://127.0.0.1:8000/api/v1/rag/embedding/status",
        "rna_folding_backend": "Invoke-RestMethod http://127.0.0.1:8000/api/v1/optimizer/rna-folding/status",
        "artifact_object_store": "Invoke-RestMethod 'http://127.0.0.1:8000/api/v1/artifacts/object-store/mirror/plan?limit=500'",
        "optimizer_diagnostics": "Invoke-RestMethod http://127.0.0.1:8000/api/v1/optimizer/diagnostics",
        "security": "Invoke-RestMethod http://127.0.0.1:8000/api/v1/security/status",
        "storage": "Invoke-RestMethod http://127.0.0.1:8000/api/v1/storage/migration/sqlite/parity",
    }
    mode_commands = {
        "artifact_refresh": "Invoke-WebRequest http://127.0.0.1:8000/api/v1/deployment/audit/export.zip -OutFile production_audit.zip",
        "data_promotion": "Invoke-RestMethod http://127.0.0.1:8000/api/v1/data/coverage",
        "managed_runtime": "python scripts/production_audit.py --require-api",
        "configuration": "python scripts/production_audit.py --require-api",
        "benchmark_calibration": "Invoke-RestMethod http://127.0.0.1:8000/api/v1/optimizer/benchmark",
        "readiness_rollup": "python scripts/production_audit.py --require-api",
    }
    command = area_commands.get(area) or mode_commands.get(resolution_mode) or "python scripts/production_audit.py --require-api"
    return {
        "artifact": f"evidence/{evidence_key}.json",
        "command": command,
    }


def proof_checklist_item(item: dict[str, Any]) -> dict[str, Any]:
    payload = dict(item)
    payload["proof_item_hash"] = hash_payload(payload)
    return payload


def verify_gap_summary(summary: dict[str, Any]) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    if not summary:
        return {
            "status": "warning",
            "errors": [],
            "warnings": ["production_gap_summary is not present; run a live production audit with --require-api for a promotion-action runbook."],
        }
    if summary.get("gap_schema") != GAP_SCHEMA:
        errors.append("production_gap_summary schema is missing or invalid.")
    gaps = summary.get("gaps") if isinstance(summary.get("gaps"), list) else []
    if safe_int(summary.get("gap_count"), -1) != len(gaps):
        errors.append("gap_count does not match gaps length.")
    if summary.get("resolution_scope_counts") != count_by(gaps, "resolution_scope"):
        errors.append("resolution_scope_counts do not match gaps.")
    if summary.get("resolution_mode_counts") != count_by(gaps, "resolution_mode"):
        errors.append("resolution_mode_counts do not match gaps.")
    proof_checklist = summary.get("proof_checklist") if isinstance(summary.get("proof_checklist"), list) else []
    expected_proof_checklist = _gap_proof_checklist(gaps)
    if safe_int(summary.get("proof_checklist_count"), -1) != len(proof_checklist):
        errors.append("proof_checklist_count does not match proof_checklist.")
    if summary.get("proof_checklist_hash") != hash_payload(proof_checklist):
        errors.append("proof_checklist_hash does not match proof_checklist.")
    if proof_checklist != expected_proof_checklist:
        errors.append("proof_checklist does not match gaps.")
    for item in proof_checklist:
        if not isinstance(item, dict):
            errors.append("proof_checklist contains a non-object item.")
            continue
        expected_hash = hash_payload({key: value for key, value in item.items() if key != "proof_item_hash"})
        if item.get("proof_item_hash") != expected_hash:
            errors.append(f"proof_item_hash does not match proof_checklist item for {item.get('area', 'unknown')}.")
    for gap in gaps:
        if not isinstance(gap, dict):
            errors.append("production_gap_summary contains a non-object gap.")
            continue
        for required in (
            "area",
            "priority",
            "resolution_scope",
            "resolution_mode",
            "proof_hint",
            "proof_artifact",
            "proof_command",
            "evidence_key",
            "action",
        ):
            if not gap.get(required):
                errors.append(f"gap {gap.get('area', 'unknown')!r} is missing {required}.")
        if gap.get("gap_hash") != hash_payload({key: value for key, value in gap.items() if key != "gap_hash"}):
            errors.append(f"gap_hash does not match gap {gap.get('area', 'unknown')!r}.")
    if summary.get("gap_summary_hash") != hash_payload({key: value for key, value in summary.items() if key != "gap_summary_hash"}):
        errors.append("gap_summary_hash does not match production_gap_summary.")
    return {
        "status": "fail" if errors else "pass",
        "errors": errors,
        "warnings": warnings,
    }


def render_markdown(runbook: dict[str, Any]) -> str:
    lines = [
        "# Production Promotion Runbook",
        "",
        f"- Runbook hash: `{runbook['runbook_hash']}`",
        f"- Source audit hash: `{runbook.get('source_audit_hash') or 'n/a'}`",
        f"- Source generated at: `{runbook.get('source_generated_at') or 'n/a'}`",
        f"- Status: `{runbook.get('status')}`",
        f"- Production ready: `{runbook.get('production_ready')}`",
        f"- Gaps: `{runbook.get('gap_count')}`; blocking `{runbook.get('blocking_count')}`; promotion `{runbook.get('promotion_count')}`",
        f"- Verification: `{(runbook.get('verification') or {}).get('status')}`",
        f"- Resolution scopes: `{format_counts(runbook.get('resolution_scope_counts') or {})}`",
        f"- Resolution modes: `{format_counts(runbook.get('resolution_mode_counts') or {})}`",
        f"- Proof checklist: `{runbook.get('proof_checklist_count')}` items; hash `{runbook.get('proof_checklist_hash')}`",
        f"- Source proof checklist: `{runbook.get('source_proof_checklist_count')}` items; hash `{runbook.get('source_proof_checklist_hash') or 'n/a'}`",
        f"- Proof checklist source match: `{runbook.get('proof_checklist_source_match')}`",
        "",
    ]
    if runbook.get("proof_checklist"):
        lines.extend(
            [
                "## Promotion Proof Checklist",
                "",
                "| Area | Priority | Scope | Mode | Artifact | Proof hash | Command |",
                "| --- | --- | --- | --- | --- | --- | --- |",
            ]
        )
        for item in runbook.get("proof_checklist") or []:
            lines.append(
                "| {area} | {priority} | {scope} | {mode} | {artifact} | {proof_hash} | {command} |".format(
                    area=escape_cell(item.get("area")),
                    priority=escape_cell(item.get("priority")),
                    scope=escape_cell(item.get("resolution_scope")),
                    mode=escape_cell(item.get("resolution_mode")),
                    artifact=escape_cell(item.get("proof_artifact")),
                    proof_hash=escape_cell(short_hash(item.get("proof_item_hash"))),
                    command=escape_cell(item.get("proof_command")),
                )
            )
        lines.append("")
    for group in runbook.get("groups") or []:
        lines.extend(
            [
                f"## {group.get('resolution_scope')} ({group.get('count')})",
                "",
                "| Area | Priority | Mode | Evidence | Proof artifact | Proof command | Proof hint | Action |",
                "| --- | --- | --- | --- | --- | --- | --- | --- |",
            ]
        )
        for item in group.get("items") or []:
            lines.append(
                "| {area} | {priority} | {mode} | {evidence} | {proof_artifact} | {proof_command} | {proof_hint} | {action} |".format(
                    area=escape_cell(item.get("area")),
                    priority=escape_cell(item.get("priority")),
                    mode=escape_cell(item.get("resolution_mode")),
                    evidence=escape_cell(item.get("evidence_key")),
                    proof_artifact=escape_cell(item.get("proof_artifact")),
                    proof_command=escape_cell(item.get("proof_command")),
                    proof_hint=escape_cell(item.get("proof_hint")),
                    action=escape_cell(item.get("action")),
                )
            )
        lines.append("")
    errors = (runbook.get("verification") or {}).get("errors") or []
    warnings = (runbook.get("verification") or {}).get("warnings") or []
    if warnings:
        lines.extend(["## Verification Warnings", ""])
        lines.extend(f"- {warning}" for warning in warnings)
        lines.append("")
    if errors:
        lines.extend(["## Verification Errors", ""])
        lines.extend(f"- {error}" for error in errors)
        lines.append("")
    return "\n".join(lines)


def count_by(items: list[Any], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in items:
        if not isinstance(item, dict):
            continue
        value = str(item.get(key) or "unknown")
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def format_counts(counts: dict[str, Any]) -> str:
    if not counts:
        return "none"
    return ", ".join(f"{key}={value}" for key, value in sorted(counts.items()))


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def hash_payload(payload: Any) -> str:
    return sha256(json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def short_hash(value: Any) -> str:
    text = str(value or "")
    return text[:12] if len(text) >= 12 else "n/a"


def escape_cell(value: Any) -> str:
    return str(value or "n/a").replace("|", "\\|").replace("\n", " ")


if __name__ == "__main__":
    raise SystemExit(main())
