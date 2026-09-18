# Agentic RAG Codon Optimization Platform

뇌 영역·세포 유형별 근거를 활용해 유전자치료용 CDS 설계를 지원하는 프로덕션 지향 연구 프로토타입입니다.

## 3분 프로젝트 요약

### 해결하려는 문제

유전자치료용 CDS 설계에는 서로 다른 생물학 문헌과 데이터, 제약조건을 고려한 서열 최적화, 품질 검증이 함께 필요합니다. 단순히 후보 서열을 생성하는 데 그치지 않고, 검색한 근거와 최적화 과정, 산출 결과를 추적하고 동일한 조건에서 재현할 수 있도록 만드는 것이 핵심 문제였습니다.

### 접근 방식

본 저장소는 다음 네 계층을 하나의 감사 가능한 워크플로로 연결합니다.

1. **근거 수집·검색** — PDF, Markdown, TXT, JSON 문서와 GTEx·Allen·사용자 정의 데이터를 정규화하고, hybrid retrieval, alias expansion, facet-aware reranking, source-priority scoring으로 필요한 근거를 검색합니다.
2. **에이전트 오케스트레이션** — planner, retriever, optimizer, QC writer 역할을 분리하고 실행 이력을 보존해, 결과만 반환하는 불투명한 단일 호출을 피했습니다.
3. **제약 기반 최적화** — seed가 고정된 NSGA-II로 동의 코돈 CDS 후보를 탐색하며 GC/CpG, motif, splice·polyA signal, restriction site, codon pair, 5-prime GC, hairpin, complexity 제약을 평가합니다.
4. **검증·출처 추적** — golden test, API contract, UI/API smoke test, 서명된 산출물 bundle, provenance hash, preflight evidence를 통해 결과를 검토하고 재현할 수 있도록 구성했습니다.

### 시스템 흐름

```mermaid
flowchart LR
    A[문서·구조화 데이터] --> B[수집·정규화]
    B --> C[Hybrid retrieval·reranking]
    C --> D[Agentic workflow]
    D --> E[Seeded NSGA-II 최적화]
    E --> F[QC·후보 비교]
    F --> G[보고서·manifested ZIP bundle]
    G --> H[보관·hash·audit·CI 검증]
```

### 주요 검토 항목

| 검토 항목 | 저장소에서 확인할 수 있는 근거 |
| --- | --- |
| 문제 해결 범위 | 근거 기반 CDS 설계부터 검색·최적화·QC·감사까지 하나의 워크플로로 구현 |
| 백엔드·UI | FastAPI API와 Next.js dashboard |
| 재현성 | seed 고정 최적화, canonical configuration hash, lockfile, snapshot, manifest |
| 품질 검증 | API contract, response/value golden test, backend regression, frontend build, browser smoke test |
| 운영 준비도 | Docker Compose, Postgres adapter, 인증·rate limit, readiness gate, operational audit, deployment runbook |
| 산출물 | Markdown·HTML·JSON·PDF 보고서와 hash 검증 가능한 ZIP evidence bundle |

### 빠른 검증 방법

저장소 루트에서 다음 명령으로 전체 검증 체인을 실행할 수 있습니다.

```powershell
python scripts/preflight.py
python scripts/portfolio_readiness_matrix.py --strict
python scripts/final_portfolio_check.py --workflow .github/workflows/ci.yml --output-dir backend/app/data/runtime --output-json backend/app/data/runtime/final_portfolio_check.json
python scripts/verify_final_portfolio_check.py --path backend/app/data/runtime/final_portfolio_check.json
```

CI는 backend compile, OpenAPI contract, golden·regression test, frontend build, browser UI 경로, 서명된 산출물, Docker Compose/Postgres 배포 구성을 함께 검증합니다.

### 기술 스택

- **Backend:** Python, FastAPI, Pydantic, SQLite/Postgres
- **Frontend:** Next.js, React, TypeScript, Playwright
- **Retrieval:** local JSON, 선택형 pgvector/Qdrant, hash 또는 OpenAI embedding
- **Operations:** Docker Compose, GitHub Actions, SHA-256 manifest, 선택형 HMAC signing

### 한계와 향후 보완점

- seed tRNA/codon-availability prior에는 신뢰도가 낮거나 placeholder인 값이 명시되어 있으며, release에 고정된 정량 데이터로 해석해서는 안 됩니다.
- 검색 품질은 수집 문서와 구조화 reference snapshot의 범위·최신성에 영향을 받습니다.
- 프로덕션 적용 전에는 검증된 source-data refresh, 배포 secret 설정, readiness·audit gate 통과가 필요합니다.

## 상세 기능

- FastAPI backend for scoring, gene-to-CDS resolution, optimization, evidence retrieval, and QC report export.
- Next.js dashboard for target entry, candidate comparison, evidence review, QC summary, data provenance, source-snapshot coverage, and operational readiness.
- MANE-aware Ensembl transcript selection.
- Local document ingestion for PDF, Markdown, text, and JSON.
- Hybrid local RAG retrieval over seed evidence, ingested documents, and structured GTEx/Allen/CUSTOM priors, including alias expansion, facet-aware reranking, source-priority scoring, and retrieval coverage evaluation.
- RAG diagnostics for corpus/source distribution, reproducible chunking policy, token and embedding health, regression macro metrics, source-provenance regression evidence, weak-case surfacing, and operational recommendations.
- Manifested RAG evaluation ZIP bundles with retrieval trace, source provenance hashes, score CSV, retrieved chunks JSONL, RAG status, structured manifest, artifact verification, and semantic cross-checks.
- Data catalog and refresh API for reproducible GTEx/Allen reference panel ingestion.
- Data provenance and release-bundle audit for source hashes, dataset/source-file summary hashes, release metadata, RAG/document coverage, and refresh history.
- Production caveat tracking for low-confidence or placeholder tRNA/codon-availability priors so seed matrices cannot be mistaken for release-pinned quantitative data.
- Baseline audit recording for local seed/ingested data so reproducibility checks can start before a live external refresh.
- Data lockfile support for pinning and drift-checking a reproducible data/RAG/document state.
- Data release lockfile support for pinning structured dataset releases, file hashes, and release drift independently from runtime refresh state.
- Data release bundle verification in deployment readiness and production audit evidence, including required file, row-count, source-byte, and promotion-status checks.
- Data snapshot ZIP export and archived semantic verification for structured data, document extracts, RAG index, refresh log, manifests, external source snapshots, and artifact-level SHA-256 audit manifests.
- External source snapshot coverage and backfill controls for reproducible bundled/manual structured records.
- Seeded NSGA-II synonymous optimizer with repair and manufacturing-policy scoring for GC/CpG/motif/polyA/splice/restriction-site/codon-pair/5-prime-GC/hairpin/secondary-structure-proxy/complexity constraints.
- Optimizer diagnostics for objective inventory, benchmark quality bands, weak-case surfacing, and tuning recommendations.
- Manifested optimizer benchmark ZIP bundles with benchmark metrics, case provenance fingerprints, recommendation summary trade-off/Pareto evidence, recommended-candidate folding evidence, diagnostics, case catalog, candidate diagnostics, config snapshots, artifact verification, archive storage, readiness/audit evidence, and semantic cross-checks.
- ORF validation and QC gate checks for CDS inputs and recommended candidates.
- QC report bundles preserve recommended-candidate RNA folding/proxy evidence with stable hashes for archive and audit review.
- Candidate diagnostics for feasible-set counts, Pareto-front representatives, score ranges, codon-level diversity, and recommendation-audit trade-off regret against best-by-metric alternatives.
- Optimizer reproducibility manifests with canonical config hashes, score-config hashes, objective inventory, seed, repair policy, target hash, and source CDS hash.
- CUSTOM-style tissue-aware codon multipliers from structured priors.
- Markdown, HTML, JSON, and PDF QC report export.
- Manifested QC report ZIP bundles with JSON, Markdown, HTML, PDF, report-format summary hashes, candidate ranking CSV, provenance files, artifact verification, and semantic cross-checks for report and optimizer reproducibility metadata.
- SQLite run persistence with request/design/report/trace storage.
- Persistent agent memory index for session, semantic-rule, and artifact summaries derived from saved design runs.
- SQLite background job store for long-running single-gene design, batch gene design, and data refresh work.
- SQLite operational audit log for report/export/verify/job/data-refresh actions with actor hashes and request IDs.
- ZIP audit bundle export for persisted runs and background jobs, including file-level SHA-256 artifact manifests.
- Artifact verification endpoints for run, job, and data snapshot export bundles.
- Governance attestation export for OpenAPI, data, RAG, storage, and artifact-ledger state with hash verification.
- Production audit export for deployment readiness, security, storage, data provenance, RAG/optimizer diagnostics, QC/data snapshot archive semantics, governance verification, artifact-ledger state, and hash-pinned required action evidence.
- Optional Postgres runtime adapter for run/job/audit stores, with redacted `DATABASE_URL` status and migration DDL.
- SQLite-to-Postgres migration bundle export, dry-run import planning, and row-level source/target parity hashing.
- Immutable local artifact archive for exported run/job/data snapshot/QC report bundles with SHA-256 indexing, archive re-verification endpoints, and QC bundle semantic re-checks.
- Append-only artifact archive ledger with hash-chain verification for stored export history.
- Dry-run-first artifact retention cleanup with ledger tombstones for auditable archive pruning.
- External timestamp provider configuration evidence for regulated artifact archive promotion review.
- JSON and Prometheus-style metrics for request, run, job, RAG, and structured-data observability.
- Local agentic workflow orchestration for planner, retriever, optimizer, and QC writer roles.
- Docker Compose configuration for backend/frontend deployment.
- API shape and value-level golden fixtures for deterministic scoring, validation, optimization, and QC regression checks.

## Local Development

Backend:

```powershell
cd backend
pip install -r requirements.txt
python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000
```

Optional environment:

```powershell
$env:APP_DATA_DIR = 'C:\path\to\persistent\data'
$env:CORS_ORIGINS = 'http://127.0.0.1:3000,http://127.0.0.1:3001'
$env:API_KEYS = 'replace-with-a-long-random-key'
$env:API_KEY_ROLES = 'replace-with-a-long-random-key=admin'
$env:RATE_LIMIT_PER_MINUTE = '120'
$env:ARTIFACT_SIGNING_KEY = 'replace-with-a-long-random-signing-secret'
$env:ARTIFACT_SIGNING_KEY_ID = 'local-hmac-sha256'
$env:ARTIFACT_RETENTION_DAYS = '0'
$env:ARTIFACT_RETENTION_KEEP_MIN = '100'
```

`API_KEYS` is optional and disabled by default for local development. When set, protected API routes accept either `X-API-Key: <key>` or `Authorization: Bearer <key>`. `API_KEY_ROLES` can map keys to `viewer`, `operator`, or `admin` roles using `key=admin;other=viewer,operator` or a JSON object. If roles are omitted, configured keys default to `admin` for backwards-compatible single-key deployments. For browser-based demos, set the matching `NEXT_PUBLIC_API_KEY` in the frontend environment.
`ARTIFACT_SIGNING_KEY` is optional. When set, exported run/job/data snapshot manifests are signed with HMAC-SHA256 and verify endpoints validate the signature.

For a production-style configuration, copy `.env.production.example`, replace every placeholder secret, set `STORAGE_BACKEND=postgres`, and confirm `/api/v1/deployment/readiness` has no `fail` gates before deployment.

Validate the production template or a filled production env file:

```powershell
python scripts/validate_production_env.py --template
python scripts/validate_production_env.py --path .env.production
python scripts/production_audit.py --template --skip-api
```

Frontend:

```powershell
cd frontend
npm install
npm run dev
```

Open `http://127.0.0.1:3001`.

## Docker

Docker is not installed in the current workspace environment, but deployment files are included:

```powershell
docker compose up --build
```

For a production-like Compose plan with Postgres, auth, retention, and artifact signing required:

```powershell
Copy-Item .env.production.example .env.production
python scripts/validate_production_env.py --path .env.production
python scripts/compose_preflight.py
docker compose --profile postgres --env-file .env.production -f docker-compose.yml -f docker-compose.production.yml up --build
```

Then open:

- Frontend: `http://127.0.0.1:3000`
- Backend docs: `http://127.0.0.1:8000/docs`

Runtime settings are available at `http://127.0.0.1:8000/api/v1/settings`.
Security status is available at `http://127.0.0.1:8000/api/v1/security/status`.
Operational audit summary is available at `http://127.0.0.1:8000/api/v1/audit/summary`.
Archived artifacts are available at `http://127.0.0.1:8000/api/v1/artifacts`.
Readiness is available at `http://127.0.0.1:8000/api/v1/health/ready`.
Reference data catalog is available at `http://127.0.0.1:8000/api/v1/data/catalog`, and the target-level live/seed coverage matrix is available at `http://127.0.0.1:8000/api/v1/data/coverage`.
Deployment notes are in `docs/DEPLOYMENT.md`.

## Verification

Run the full local preflight from the repository root:

```powershell
python scripts/preflight.py
python scripts/preflight.py --output-json backend/app/data/runtime/preflight_latest.json
```

The preflight compiles the backend, validates the production env template, checks Docker Compose deployment shape, checks the OpenAPI contract, verifies the structured import CLI preview path, verifies the GTEx/Allen data-refresh CLI planning and validation paths, runs response/value golden tests, runs the final portfolio evidence check and its verifier, runs the manual backend regression suite, builds the frontend, starts the backend and frontend if needed, runs the HTTP smoke test, runs a separate signed-artifact smoke profile with `ARTIFACT_SIGNING_KEY` enabled, and verifies the browser UI smoke path. The readiness smoke also verifies that tRNA/codon-availability seed caveats are exposed through the data provenance gate, and the artifact smoke checks indexed semantic summaries including workflow trace archives. Use `--skip-frontend`, `--skip-smoke`, `--skip-signing-smoke`, or `--skip-ui-smoke` for narrower diagnostics. `--output-json` writes a reproducible evidence file with timestamps, executed commands, durations, skipped checks, pass/fail status, compact stdout/stderr, parsed JSON details where commands emit JSON, `mode_hash`, `skipped_hash`, `checks_hash`, and `preflight_hash`; the production audit repeats those hashes in Markdown for promotion review.

For a deployment audit artifact, run:

```powershell
python scripts/production_audit.py --path .env.production --require-api
python scripts/production_audit.py --path .env.production --require-api --preflight-evidence backend/app/data/runtime/preflight_latest.json
python scripts/production_promotion_runbook.py --write-result backend/app/data/runtime/production_audit_write_result.json --output backend/app/data/runtime/production_promotion_runbook.md
```

This writes JSON and Markdown reports under `backend/app/data/runtime/production_audits`, combining production env validation, Compose checks, optional preflight evidence validation, and live operational API status when the backend is reachable. When `--preflight-evidence` is provided, the audit requires a passing, recent preflight payload with the expected backend, contract, data-refresh, golden, and regression checks.
The promotion runbook command accepts either the audit write-result or a `production_audit.json`, verifies the gap hashes, and groups remaining actions by operator environment, source-data promotion, artifact refresh, and optimizer validation scope.

Individual contract and golden checks are also available from the repository root:

```powershell
python scripts/portfolio_readiness_matrix.py --strict
python scripts/final_portfolio_check.py --workflow .github/workflows/ci.yml --output-dir backend/app/data/runtime --output-json backend/app/data/runtime/final_portfolio_check.json
python scripts/verify_final_portfolio_check.py --path backend/app/data/runtime/final_portfolio_check.json
python scripts/api_contract_test.py
python scripts/golden_response_test.py
python scripts/golden_value_test.py
```

`portfolio_readiness_matrix.py` emits a hash-pinned JSON evidence matrix for the PDF scope: real data ingestion, Agentic RAG retrieval/evaluation, multi-objective codon optimization, QC export, operator UI, deployment operations, and automated audit/CI coverage.
`final_portfolio_check.py` regenerates the matrix, CI production evidence chain, and portfolio submission summary artifacts, then records their SHA-256 values and payload hashes in `final_portfolio_check.json`; `verify_final_portfolio_check.py` recomputes those hashes before portfolio submission.

Backend smoke and manual regression checks can be run from `backend`:

```powershell
cd backend
python -c "from tests import test_optimizer as t; [getattr(t, name)() for name in dir(t) if name.startswith('test_')]; print('manual tests passed')"
python scripts/smoke_test.py
cd ..\frontend
npm.cmd run smoke:ui
```

Install `pytest` if you want standard test discovery.
