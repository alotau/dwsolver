# Data Model: Architecture Diagrams Directory

**Feature**: 007-architecture-diagrams
**Date**: 2026-03-21

This feature is pure documentation. There is no persistent data store or runtime schema. The "entities" here are the conceptual elements that each diagram must accurately represent. This document defines those entities and their required attributes so that diagram authors have a canonical reference.

---

## Entity 1: Thread

Represents one concurrent execution unit in the dwsolver runtime.

| Attribute | Description |
|-----------|-------------|
| `id` | `master` or `subproblem-{0..N-1}` |
| `role` | `coordinator` (master) or `worker` (subproblem) |
| `entry_function` | `main()` for master; `subproblem_thread()` for workers |
| `states` | Set of discrete states the thread passes through (see State entity) |

**Validation**: Every box in the threading flow diagram must be attributable to exactly one thread.

---

## Entity 2: ThreadState

A discrete phase of a thread's lifecycle. Used to populate the `stateDiagram-v2` nodes.

| State name | Owner thread | Description |
|------------|-------------|-------------|
| `Init` | master | Reading guidefile, initialising globals |
| `SpawnSubproblems` | master | `pthread_create` × N |
| `BuildMasterLP` | master | Setting up initial reduced master problem |
| `[Phase1Region]` | master | Outer loop of `phase_1_iteration` calls (conditional) |
| `WaitN_Ph1` | master | `sem_wait(&customers)` × N per Phase 1 iteration |
| `Drain` | master | Discarding N stale results at Phase 1→2 transition |
| `PushPhase2Duals` | master | Writing fresh duals after drain |
| `[Phase2Region]` | master | Outer loop of `phase_2_iteration` calls |
| `WaitN_Ph2` | master | `sem_wait(&customers)` × N per Phase 2 iteration |
| `FinalSolve` | master | Last master LP solve after column generation |
| `OptionalRounding` | master | `dw_rounding` post-processing (if `rounding_flag`) |
| `WriteOutputs` | master | Writing solution files |
| `InitialSolve` | subproblem-N | GLPK simplex on initial subproblem LP |
| `WaitMasterReady` | subproblem-N | `pthread_cond_wait(&master_lp_ready_cv, ...)` |
| `SetupColumns` | subproblem-N | Building column structures from master, computing initial duals |
| `SolveLP` | subproblem-N | GLPK simplex (or `glp_intopt` if integer mode) per iteration |
| `EnqueueSelf` | subproblem-N | Lock `service_queue_mutex`, enqueue `my_id`, `sem_post(&customers)` |
| `WaitNextIteration` | subproblem-N | `pthread_cond_wait(&next_iteration_cv, ...)` |
| `Done` | all | Thread exits |

---

## Entity 3: SynchronizationPrimitive

A concurrency control object. Every primitive must appear as a labeled edge or annotation in the threading flow diagram.

| Name | Type | Direction | Used by |
|------|------|-----------|---------|
| `master_lp_ready_mutex` | mutex | — | master (lock/unlock), subproblem-N (lock/unlock) |
| `master_lp_ready_cv` | condition variable | master→subproblems | broadcast by master; waited by subproblem-N |
| `customers` | semaphore | subproblem-N→master | posted by subproblem-N; waited by master |
| `service_queue_mutex` | mutex | — | master (read head), subproblem-N (write tail) |
| `next_iteration_mutex` | mutex | — | master (lock/unlock), subproblem-N (lock/unlock) |
| `next_iteration_cv` | condition variable | master→subproblems | broadcast by master; waited by subproblem-N |
| `master_mutex` | mutex | — | master (writes duals), subproblem-N (reads duals) |
| `sub_data_mutex[i]` | per-subproblem mutex | — | master (writes `r`), subproblem-N (reads `r`) |

---

## Entity 4: Module

A dwsolver source module. Appears as a node in the component diagram.

| Module | Source file | Primary responsibility |
|--------|-------------|----------------------|
| `dw_main` | `src/dw_main.c` | Entry point; orchestrates startup, Phase 1/2 loops, output |
| `dw_phases` | `src/dw_phases.c` | `phase_1_iteration()` and `phase_2_iteration()` logic |
| `dw_subprob` | `src/dw_subprob.c` | `subproblem_thread()` — worker thread entry; LP solve + signal |
| `dw_support` | `src/dw_support.c` | Init, CLI parsing, D-matrix setup, solution output, logging |
| `dw_globals` | `src/dw_globals.c` | Defines all global variables (sync primitives, LP objects) |
| `dw_blas` | `src/dw_blas.c` | Portable BLAS-like sparse matrix operations |
| `dw_rounding` | `src/dw_rounding.c` | Optional post-solve integer rounding |
| GLPK | `src/glp*.c` (vendored) | LP/MIP solver; treated as a single external block |

---

## Entity 5: InputFile

A file consumed at startup. Drives the data-flow diagram.

| File | Format | Consumed by | Order in guidefile |
|------|--------|-------------|-------------------|
| guidefile | Plain text (line-oriented) | `process_cmdline` in `dw_support` | CLI arg `-g` |
| subproblem LP | GLPK LP format | `subproblem_thread` in `dw_subprob` | Lines 2..N+1 of guidefile |
| master LP | GLPK LP format | `main()` in `dw_main` | Line N+2 of guidefile |

---

## Entity 6: OutputFile

A file produced at completion. Terminates the data-flow diagram.

| File | Condition | Content |
|------|-----------|---------|
| `relaxed_solution` | Always | Continuous optimal variable values |
| `zeros_rounded` | Only if `-r`/`--round` is set | Variable values after rounding heuristic applied to the relaxed solution |
| `integerized_zeros` | Only if `-i`/`--integerize` is set | Variables fixed to zero after integerization |
| `intprob_iter_*.cpxlp` | Only if `--write-int-probs` is set | Intermediate LP snapshots for integer subproblems |
| `basis_iteration_*` | Only if `--write-bases` is set | Per-iteration LP basis snapshots |

---

## State transitions (threading flow)

```
master_lp_ready_cv broadcast:
  SpawnSubproblems → BuildMasterLP (master)
  InitialSolve → WaitMasterReady → SetupColumns (subproblem-N)

customers sem_post:
  SolveLP → EnqueueSelf (subproblem-N posts)

customers sem_wait:
  WaitN_Ph1 or WaitN_Ph2 (master waits N times)

next_iteration_cv broadcast:
  Phase1 or Phase2 iteration end → (master broadcasts)
  WaitNextIteration → SolveLP (subproblem-N wakes)

Drain transition (Phase 1→2):
  master drains N results from service_queue (no dual push)
  then pushes Phase 2 duals, increments current_iteration, rebroadcasts
```
