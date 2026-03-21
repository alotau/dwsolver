# Implementation Plan: Architecture Diagrams Directory

**Branch**: `007-architecture-diagrams` | **Date**: 2026-03-21 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/007-architecture-diagrams/spec.md`

## Summary

Create an `architecture/` directory at the project root containing three Mermaid diagrams (threading/control-flow, component map, input data flow) plus an index README, and update the root `README.md` to link there. All content is pure documentation — no source code changes, no build system changes.

## Technical Context

**Language/Version**: Markdown / Mermaid diagram syntax (GitHub-native rendering, no version dependency)
**Primary Dependencies**: None — Mermaid renders natively in GitHub Markdown as of 2022
**Storage**: Plain `.md` files committed to the repository
**Testing**: Visual inspection in GitHub Markdown preview; no programmatic test runner needed
**Target Platform**: GitHub.com Markdown renderer
**Project Type**: Documentation
**Performance Goals**: N/A
**Constraints**: Diagrams must parse and render without errors in GitHub's Mermaid renderer; no external images or CDN resources
**Scale/Scope**: 4 new files (`architecture/README.md` + 3 diagram files) + 1 updated file (`README.md`)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | ✅ PASS | No solver logic changed; test suite unaffected |
| II. Thread Safety | ✅ PASS | No code changes; diagrams document existing primitives |
| III. Cross-Platform | ✅ PASS | Markdown/Mermaid has no platform dependency |
| IV. Repair Before Extension | ✅ PASS | This is documentation of existing behavior, not a new solver feature |
| V. CLI-First, Library-Ready | ✅ PASS | No API or interface changes |

All gates pass. No violations to justify.

## Project Structure

### Documentation (this feature)

```text
specs/007-architecture-diagrams/
├── plan.md              # This file
├── research.md          # Phase 0 output — threading, module, and data-flow research
├── data-model.md        # Phase 1 output — diagram entity model
├── quickstart.md        # Phase 1 output — how to update diagrams
├── contracts/           # Phase 1 output — per-diagram content contracts
│   ├── threading-flow-contract.md
│   ├── component-map-contract.md
│   └── data-flow-contract.md
└── tasks.md             # Phase 2 output (/speckit.tasks — NOT created here)
```

### Deliverable files (repository root)

```text
architecture/
├── README.md                  # Index: one-sentence description per diagram
├── threading-flow.md          # FR-002, FR-003, FR-006, FR-007 — Mermaid stateDiagram-v2
├── component-map.md           # FR-004, FR-006, FR-007 — Mermaid graph LR
└── data-flow.md               # FR-005, FR-006, FR-007 — Mermaid flowchart TD

README.md                      # FR-009 — add link to architecture/ directory
```

**Structure Decision**: Single flat `architecture/` directory at the project root. No subdirectories needed at this scope (3 diagrams). The `.md` extension means each file is a self-contained Markdown document containing both the Mermaid block and its prose explanation (FR-007).

## Phase 0 Research

*See [research.md](./research.md) for full findings. Key decisions:*

- **Diagram format**: Mermaid — GitHub native, no plugin, version-stable syntax
- **Threading diagram type**: `stateDiagram-v2` — best fit for showing thread states and transitions gated by sync primitives
- **Component diagram type**: `graph LR` — left-to-right call graph between modules
- **Data flow diagram type**: `flowchart TD` — top-down input → transform → output flow

## Phase 1 Design

*See [data-model.md](./data-model.md), [quickstart.md](./quickstart.md), and [contracts/](./contracts/) for full detail.*

### Diagram content summary

**threading-flow.md** (`stateDiagram-v2`):
- Two swimlane groups: Master Thread, Subproblem Thread-N
- States: Init → InitialSolve → WaitMasterReady → [Phase 1 region] → Drain → [Phase 2 region] → Done
- All sync primitives labeled on transitions: `master_lp_ready_cv`, `customers` semaphore, `service_queue_mutex`, `next_iteration_cv`, `sub_data_mutex[i]`, `master_mutex`
- Phase 1 region annotated as conditional (only when auxiliary variables detected)

**component-map.md** (`graph LR`):
- 7 dwsolver module nodes + 1 GLPK external block
- Edges labeled with primary call relationship: `calls`, `spawns`, `reads globals`, `optional post-solve`

**data-flow.md** (`flowchart TD`):
- Input tier: guidefile → `num_clients` count + N subproblem filenames + 1 master filename
- Process tier: initial solves → Phase 1 (conditional) → Phase 2 iterations → final master solve → optional rounding
- Output tier: `relaxed_solution`, `zero_vars`, optional `integer_solution`
