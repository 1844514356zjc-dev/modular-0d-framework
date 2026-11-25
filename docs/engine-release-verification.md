# Engine release verification

Local checks performed on 2026-11-15, before publishing this repository.

| Check | Environment / method | Result |
| --- | --- | --- |
| Runtime lock installation | Python 3.10.22; `uv pip install --require-hashes --only-binary=:all: -r requirements.txt` | Passed |
| Dependency consistency | `uv pip check --python .venv-py310/bin/python` | Passed |
| PEP 8 | `python -m flake8 src tests benchmarks examples --count --show-source --statistics` | 0 violations |
| Complete tests | Python 3.10.22; `python -m pytest -ra --strict-config --strict-markers` | 68 passed |
| Companion-cycle references | `python benchmarks/validate_literature.py` | All comparison rows passed their declared tolerances |
| Linux x86-64 wheels | Python 3.10, manylinux_2_28, hash-checked binary-only dry run | Resolved 40 runtime packages |
| Deployment YAML | Parse workflow and Compose files; check trigger and service port declarations | Passed |
| Base image availability | Docker Registry manifest for `python:3.10-slim-bookworm` | Available |
| README local links | Resolve relative Markdown targets | Passed |
| Container build and runtime | Docker CLI unavailable on this workstation | Not executed locally |
| GitHub Actions execution | Workflow prepared; repository not published by this task | Not executed remotely |

The workflow includes a container build, Compose health check, engine/scavenging
smoke calculation, and companion-cycle smoke calculations. Those configured
checks are not reported here as completed runs. Dependency resolution and YAML
parsing do not substitute for a Docker build.

These checks establish numerical/software behavior for the tested conditions;
they do not constitute experimental engine validation. See the root README for
model scope and the distinction between engine ideal-gas properties and the
companion laboratory's CoolProp calculations.

## English dashboard refactor

Subsequent checks on 2026-11-15:

- Python 3.10: full suite **71 passed**; Flake8 across source, tests, examples
  and benchmarks: **0 violations**.
- Focused dashboard suite rerun after the final download behavior update:
  **6 passed**.
- Headless Chrome/Playwright: English-only rendered text, simulation button,
  two-case RPM study, all four CSV/JSON downloads, log P–V toggle, and mobile
  width at 390 px passed.
- Downloaded simulation JSON retained both sweep cases and all 1,441 samples.
- Desktop thermodynamic/gas-exchange/sweep views and mobile header/KPIs were
  captured and visually inspected; sidebar expander contrast was corrected.

These checks concern the default engine application. The separately preserved
real-fluid dashboard retains its earlier interface.

## Modular 0D composition refactor

Checks performed on 2026-11-15 after extracting mechanism, exchange, property,
and heat-model components:

- Full Python 3.10 suite: **83 passed**, including legacy adapters, independent
  slider-crank derivatives, valve 720° timing, mixed port/valve composition,
  dictionary input preservation, CoolProp derivative/entropy checks, analytic
  Otto efficiency and dynamic UI topology controls.
- Flake8 across source/tests/benchmarks/examples: **0 violations**; Ruff passed.
- Headless composition example ran all three mechanism/exchange arrangements.
  Air-standard Otto efficiency: 0.6018928468, versus 0.6018928294 analytically
  (compression ratio 10, gamma 1.4).
- Native Chrome/Playwright: opposed-piston simulation, parameter sweep and four
  downloads; four-stroke valves and CoolProp Air; closed Otto energy audit;
  closed Nitrogen gas spring; log P–V axes; English text and mobile width.
- Rendered four-stroke flow areas and closed-cycle energy budgets were visually
  inspected after waiting for Plotly SVG rendering.
- Runtime component examples using Air, Nitrogen and CO2 ran for both open
  topologies at the tested 10 mg heat-source setting. This is not a guarantee
  for arbitrary state ranges; the EOS validity and gas-phase guards apply.

No new runtime dependency was introduced. Docker/GitHub remote execution was
not performed locally; the existing CI container smoke step now exercises the
opposed, four-stroke and closed presets.
