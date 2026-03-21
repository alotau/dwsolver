---

description: "Task list for feature 007-architecture-diagrams implementation"
---

# Tasks: Architecture Diagrams Directory

**Input**: Design documents from `/specs/007-architecture-diagrams/`
**Prerequisites**: plan.md ✅, spec.md ✅, research.md ✅, data-model.md ✅, contracts/ ✅, quickstart.md ✅

**Tests**: Not applicable — this is a documentation feature. Acceptance is verified by visual inspection in GitHub Markdown preview (SC-002) and by PR review against contract acceptance tests.

---

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies on incomplete tasks)
- **[Story]**: Which user story this task belongs to (US1, US2, US3)

---

## Phase 1: Setup

**Purpose**: Create the `architecture/` directory at the project root (FR-001)

- [ ] T001 Create `architecture/` directory at project root (satisfies FR-001)

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Skeleton index and root README link — must exist before user story phases complete

**⚠️ CRITICAL**: Both files must be created before Phase 3 begins. Diagram entries will be added to `architecture/README.md` incrementally in each user story phase.

- [ ] T002 [P] Create `architecture/README.md` with title, one-line intro, and three placeholder diagram entries (threading-flow, component-map, data-flow) per FR-008
- [ ] T003 [P] Add "Architecture" section with link to `architecture/` to root `README.md` per FR-009 — place it after the existing badges/intro block

**Checkpoint**: Foundation ready — user story phases can now proceed

---

## Phase 3: User Story 1 — Threading and control-flow diagram (Priority: P1) 🎯 MVP

**Goal**: `architecture/threading-flow.md` containing a valid `stateDiagram-v2` Mermaid diagram and all required prose sections, covering the full master + subproblem threading lifecycle including the Phase 1→2 drain step

**Independent Test**: A reviewer unfamiliar with the source can correctly identify all 9 sync primitives and describe the master-subproblem handshake by reading only this file (SC-001); diagram renders without errors in GitHub preview (SC-002)

### Implementation for User Story 1

- [ ] T004 [P] [US1] Write Mermaid `stateDiagram-v2` block in `architecture/threading-flow.md` — include all 18 required states (12 master + 6 subproblem), the `state "Phase 1 (conditional)" as Phase1` region, the standalone Drain + PushPhase2Duals transition states, the `state "Phase 2 — column generation" as Phase2` region, and all 9 sync primitive labels (master_lp_ready_cv, customers semaphore, service_queue_mutex, next_iteration_cv, master_mutex, sub_data_mutex[i]) per `specs/007-architecture-diagrams/contracts/threading-flow-contract.md`
- [ ] T005 [US1] Add all required prose sections below the diagram in `architecture/threading-flow.md`: (1) one-paragraph introduction, (2) "How to read this diagram" bullet walkthrough of one full DW iteration naming each sync primitive in fire order, (3) "Phase 1→2 transition" paragraph explaining the drain step and referencing the race it fixed, (4) "Single-subproblem variant (N=1)" note per `specs/007-architecture-diagrams/contracts/threading-flow-contract.md`
- [ ] T006 [US1] Update `architecture/README.md` — replace the threading-flow placeholder entry with a complete one-sentence description linking to `threading-flow.md`

**Checkpoint**: US1 complete — threading-flow.md is fully readable and passes all 5 contract acceptance tests

---

## Phase 4: User Story 2 — Module component map (Priority: P2)

**Goal**: `architecture/component-map.md` containing a valid `graph LR` Mermaid diagram and all required prose sections, covering all 7 dwsolver modules plus GLPK as an external block with labeled call/dependency edges

**Independent Test**: A developer can identify which source file handles building the master LP column from a subproblem solution by reading only this file (spec US2 independent test); diagram renders without errors (SC-002)

### Implementation for User Story 2

- [ ] T007 [P] [US2] Write Mermaid `graph LR` block in `architecture/component-map.md` — include all 8 nodes (dw_main, dw_phases, dw_subprob, dw_support, dw_globals with distinct shape, dw_blas, dw_rounding, GLPK 4.44), `subgraph External` wrapping GLPK, `subgraph Modules` wrapping the 7 dwsolver nodes, all required labeled edges, and a shared-state annotation for dw_globals per `specs/007-architecture-diagrams/contracts/component-map-contract.md`
- [ ] T008 [US2] Add all required prose sections below the diagram in `architecture/component-map.md`: (1) one-paragraph introduction, (2) "Module responsibilities" with 1–2 sentences per node, (3) "GLPK boundary" explanation of vendored location and why shown as single block, (4) "dw_rounding: optional post-processing" note per `specs/007-architecture-diagrams/contracts/component-map-contract.md`
- [ ] T009 [US2] Update `architecture/README.md` — replace the component-map placeholder entry with a complete one-sentence description linking to `component-map.md`

**Checkpoint**: US2 complete — component-map.md passes all 4 contract acceptance tests

---

## Phase 5: User Story 3 — Input data flow diagram (Priority: P3)

**Goal**: `architecture/data-flow.md` containing a valid `flowchart TD` Mermaid diagram and all required prose sections, covering guidefile → LP inputs → all process stages → output files, with conditional Phase 1 and rounding paths clearly annotated

**Independent Test**: A first-time user can prepare a valid guidefile and correctly ordered LP files for a two-subproblem problem by reading this file alongside the README (spec US3 independent test); diagram renders without errors (SC-002)

### Implementation for User Story 3

- [ ] T010 [P] [US3] Write Mermaid `flowchart TD` block in `architecture/data-flow.md` — include all 15 required nodes in 6 subgraphs (Input files, Startup, Phase 1 conditional, Phase 2, Finalization, Outputs), the `BUILD_MASTER -.->|no auxiliary variables| PH2` dashed bypass edge, `[if aux vars needed]` label on the Phase 1 edge, `[if rounding_flag]` label on the rounding path, and a guidefile format annotation on the GF node per `specs/007-architecture-diagrams/contracts/data-flow-contract.md`
- [ ] T011 [US3] Add all required prose sections below the diagram in `architecture/data-flow.md`: (1) one-paragraph introduction, (2) "Guidefile format" section with verbatim line format and a two-subproblem example, (3) "Phase 1: when does it run?" short explanation, (4) "Output files" table mapping each output filename to its condition per `specs/007-architecture-diagrams/contracts/data-flow-contract.md`
- [ ] T012 [US3] Update `architecture/README.md` — replace the data-flow placeholder entry with a complete one-sentence description linking to `data-flow.md`

**Checkpoint**: US3 complete — data-flow.md passes all 5 contract acceptance tests; architecture/README.md now has all three final entries

---

## Phase 6: Polish and Cross-Cutting Concerns

**Purpose**: Render verification, final index review, FR/SC cross-check

- [ ] T013 Verify all three Mermaid diagrams render without errors — open each file in GitHub Markdown preview or paste each diagram block into mermaid.live; fix any syntax errors (SC-002)
- [ ] T014 Final review of `architecture/README.md`: confirm all three diagram entries are present with accurate one-sentence descriptions and working relative links (FR-008, SC-003)
- [ ] T015 Cross-check all requirements against delivered files: FR-001 (dir exists), FR-002 (threading states), FR-003 (Phase 1 region + Drain state), FR-004 (7 modules + GLPK), FR-005 (data flow), FR-006 (Mermaid syntax), FR-007 (prose in each file), FR-008 (architecture/README.md), FR-009 (root README link), SC-001–SC-004

---

## Dependency Graph

```
T001 (create dir)
  └── T002 (architecture/README.md skeleton)   ─┐ parallel
  └── T003 (root README.md link)               ─┘

T002, T003 complete → user story phases can begin

T004 (threading Mermaid) ─┐ all parallel (different files)
T007 (component Mermaid)  ─┤
T010 (data-flow Mermaid)  ─┘

T004 → T005 (threading prose, same file)
T005 → T006 (README.md threading entry)

T007 → T008 (component prose, same file)
T008 → T009 (README.md component entry)

T010 → T011 (data-flow prose, same file)
T011 → T012 (README.md data-flow entry)

T006, T009, T012 → T013, T014, T015 (polish — all artifacts complete)
```

## Parallel Execution Examples

**Foundation (can run together)**:
- T002 + T003 — different files, no inter-dependency

**Diagram blocks (can all run together after T001 is done)**:
- T004 + T007 + T010 — three different files, no inter-dependency

**Polish (can run together after T012)**:
- T013 + T014 + T015 — independent review tasks

## Implementation Strategy

**MVP scope (just US1)**: Complete T001 → T002 → T003 → T004 → T005 → T006 and the repo has a usable threading-flow diagram linked from root README. This alone satisfies SC-001 and makes the hardest-to-understand part of the codebase accessible.

**Full delivery**: All 15 tasks complete → 3 diagrams + index + root link → all FRs and SCs satisfied.

## Task Count Summary

| Phase | Tasks | User story |
|-------|-------|-----------|
| Setup | 1 (T001) | — |
| Foundational | 2 (T002–T003) | — |
| Phase 3 | 3 (T004–T006) | US1 (P1) |
| Phase 4 | 3 (T007–T009) | US2 (P2) |
| Phase 5 | 3 (T010–T012) | US3 (P3) |
| Polish | 3 (T013–T015) | — |
| **Total** | **15** | |

Parallel opportunities identified: 3 (foundation pair, diagram block triple, polish triple)
