# Research: Architecture Diagrams Directory

**Feature**: 007-architecture-diagrams
**Date**: 2026-03-21

---

## Decision 1: Diagram format

- **Decision**: Mermaid in `.md` files
- **Rationale**: GitHub has rendered Mermaid natively since 2022 (no viewer action, no extension, no external hosting). Plain text — diffs well, versioned with the source, renders on the repo landing page.
- **Alternatives considered**: PlantUML (requires server or plugin), DOT/Graphviz (not natively rendered on GitHub), SVG exports (binary, doesn't diff), draw.io XML (binary, requires plugin)

---

## Decision 2: Threading diagram type

- **Decision**: `stateDiagram-v2`
- **Rationale**: The threading model is fundamentally a state machine — each thread moves through discrete states (Idle, Solving, WaitingForMaster, InQueue, WaitingForNextIteration). `stateDiagram-v2` supports `note` annotations for sync primitives and can be divided into concurrent regions to separate master and subproblem activity lanes.
- **Alternatives considered**: `sequenceDiagram` (too narrow — shows message sequence but not state ownership or the N-way fan-out); `flowchart` (doesn't capture state ownership per thread clearly)

---

## Decision 3: Component diagram type

- **Decision**: `graph LR` (left-to-right directed graph)
- **Rationale**: Module dependencies are a DAG with a clear entry point (`dw_main`) at left. LR orientation makes the call hierarchy readable without crossing edges.
- **Alternatives considered**: `classDiagram` (implies OOP relationships that don't apply in C); `graph TD` (top-down works but the chain is wider than it is tall)

---

## Decision 4: Data flow diagram type

- **Decision**: `flowchart TD` (top-down flowchart)
- **Rationale**: Input → process → output is a natural vertical flow. Top-down reads as "data enters at top, results fall to bottom." Supports subgraphs for the conditional Phase 1 region.
- **Alternatives considered**: `sequenceDiagram` (appropriate only for time-ordered message passing, not data transformation)

---

## Research: Threading model (from source)

### Synchronization primitives (all defined in `dw_globals.c`)

| Primitive | Type | Purpose |
|-----------|------|---------|
| `master_lp_ready_mutex` + `master_lp_ready_cv` | mutex + condition variable | Gates subproblem threads until master LP is set up |
| `customers` | semaphore (named macOS / unnamed Linux) | Subproblem posts when it has a result ready; master waits N times per iteration |
| `service_queue_mutex` | mutex | Protects the circular `service_queue` array (tail write by subproblem, head read by master) |
| `next_iteration_mutex` + `next_iteration_cv` | mutex + condition variable | Master broadcasts to wake all subproblems for next iteration |
| `master_mutex` | mutex | Protects master LP dual-value reads from subproblems |
| `sub_data_mutex[i]` | per-subproblem mutex array | Protects `sub_data[i].r` convexity-constraint dual update |
| `reduced_cost_mutex` | mutex | Protects reduced-cost check shared state |
| `glpk_mutex` | mutex | Protects GLPK global state (GLPK is not fully thread-safe internally) |
| `fputs_mutex` | mutex | Protects logging/output (serializes `dw_printf` calls) |

### Threading flow (one full DW iteration)

**Startup:**
1. Master thread: reads guidefile, spawns N `subproblem_thread` threads via `pthread_create`
2. Each subproblem thread: reads its LP file, runs initial GLPK simplex solve
3. Each subproblem thread: checks `master_lp_ready` flag — if 0, waits on `master_lp_ready_cv`
4. Master thread: builds master LP structure, sets `master_lp_ready = 1`, broadcasts `master_lp_ready_cv`
5. Subproblem threads wake, complete column-structure setup

**Per-iteration (Phase 1 or Phase 2) — master side:**
1. Master calls `phase_1_iteration()` or `phase_2_iteration()`
2. Master waits: `sem_wait(&customers)` — blocks until a subproblem posts
3. Master reads: locks `service_queue_mutex`, reads `service_queue[head]`, advances head, unlocks
4. Master processes the subproblem result (checks reduced cost, builds master column)
5. Repeats steps 2–4 for all N subproblems
6. Solves master LP, extracts dual values
7. Pushes updated duals into `sub_data[i]` (under `sub_data_mutex[i]`) and `row_duals` (under `master_mutex`)
8. Increments `current_iteration`, broadcasts `next_iteration_cv`

**Per-iteration — subproblem side (`signal_availability()`):**
1. Subproblem solves its LP (GLPK simplex)
2. Subproblem: locks `service_queue_mutex`, enqueues `my_id` in `service_queue[tail]`, advances tail, unlocks
3. Subproblem: posts `sem_post(&customers)` — signals master a result is ready
4. Subproblem: locks `next_iteration_mutex`, waits `pthread_cond_wait(&next_iteration_cv, ...)` until `current_iteration` advances, unlocks

**Phase 1 → Phase 2 transition drain (introduced in fix-007):**
1. Master sets `sub_data[j].phase_one = 0` for all j (marking Phase 2 active)
2. Master drains N in-flight results without using them: `for i in 0..N { sem_wait; dequeue from service_queue; }`
3. Master pushes fresh Phase 2 dual values from the just-solved master LP into all `sub_data[i]`
4. Master increments `current_iteration`, rebroadcasts `next_iteration_cv`
5. Subproblems wake with `phase_one=0` and correct Phase 2 duals, begin first real Phase 2 iteration

---

## Research: Module call graph (from `#include` and call-site analysis)

| Caller | Callee | Relationship |
|--------|--------|-------------|
| `dw_main` | `dw_support` | `init_globals`, `init_signals`, `process_cmdline`, `prepare_D`, `prepare_md`, `get_solution`, `dw_printf` |
| `dw_main` | `dw_subprob` | `pthread_create` with `subproblem_thread` as entry point |
| `dw_main` | `dw_phases` | `phase_1_iteration`, `phase_2_iteration` |
| `dw_main` | `dw_rounding` | `check_col_integrality` (optional; only when `rounding_flag` set) |
| `dw_main` | GLPK | Master LP creation, column/row manipulation, simplex solve |
| `dw_phases` | `dw_subprob` | Reads subproblem results from shared `sub_data` |
| `dw_phases` | `dw_support` | `dw_printf` |
| `dw_phases` | GLPK | Reads dual values from master LP |
| `dw_subprob` | `dw_blas` | Sparse matrix ops (`dw_daxpy`, etc.) for column construction |
| `dw_subprob` | `dw_support` | `organize_solution`, `prepare_column` |
| `dw_subprob` | GLPK | Subproblem LP solve (simplex / MIP) |
| `dw_rounding` | `dw_blas` | Sparse matrix ops |
| `dw_rounding` | `dw_support` | Utility functions |
| `dw_rounding` | GLPK | LP manipulation for rounded solution |
| `dw_support` | GLPK | Parses LP files, utility calls |
| `dw_globals` | — | Defines all global variables; no calls (pure data) |
| `dw_blas` | — | Pure math utilities; no external deps beyond C stdlib |

**Note on `dw_globals`**: It is a definitions-only translation unit — no functions, purely variable definitions. All other modules access globals implicitly via `dw.h` extern declarations. It is shown in the component diagram as a shared state node that all modules read/write, not as a caller/callee.

---

## Research: Input file format and data flow

### Guidefile format
```
<N>                          ← line 1: number of subproblems
<subproblem_1.lp>            ← lines 2..N+1: subproblem LP filenames (one per line)
...
<subproblem_N.lp>
<master.lp>                  ← line N+2: master LP filename
[optional extra lines ignored]
```

### Data flow summary

```
Inputs:
  guidefile ──► num_clients (N), subproblem filenames, master filename
  subproblem_i.lp ──► each subproblem LP (GLPK format)
  master.lp ──► master LP structure (GLPK format)

Process:
  1. Parse guidefile (dw_support/process_cmdline)
  2. Load master LP & N subproblem LPs via GLPK (dw_main, dw_subprob)
  3. Each subproblem: initial GLPK simplex solve
  4. Master: build initial reduced master LP (add convexity constraints)
  [Phase 1 — conditional, only if auxiliary variables detected]
  5. Phase 1 iterations: column generation eliminating auxiliary vars
     (master: phase_1_iteration; subproblems: solve + signal_availability)
  6. Phase 1→2 drain: discard stale results, push Phase 2 duals
  [Phase 2 — always]
  7. Phase 2 iterations: main DW column generation
     (master: phase_2_iteration; subproblems: solve + signal_availability)
  8. Termination check: reduced cost ≤ tolerance or max iterations reached
  9. Final master LP solve

Outputs:
  relaxed_solution   ← continuous optimal (always)
  zero_vars          ← variables at zero in optimal basis (always)
  zeros_rounded      ← rounded solution from zeros (only if -r/--round set)
  integerized_zeros  ← integerized solution from zeros (only if -i/--integerize set)
  basis snapshots    ← per-iteration LP bases (only if --write-bases set)
  int LP dumps       ← intermediate integer LPs (only if --write-int-probs set)
```

---

## NEEDS CLARIFICATION: none

All research items resolved. No blockers for Phase 1 design.
