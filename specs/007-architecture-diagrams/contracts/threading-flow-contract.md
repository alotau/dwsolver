# Contract: threading-flow.md

**Diagram type**: `stateDiagram-v2`
**File**: `architecture/threading-flow.md`
**Requirements covered**: FR-002, FR-003, FR-006, FR-007

---

## Required nodes (states)

Every state in this list MUST appear in the diagram. Label text may be condensed but must be unambiguous.

### Master thread states
- `Init` — reading guidefile, initializing globals
- `SpawnSubproblems` — `pthread_create` × N
- `BuildMasterLP` — GLPK master LP construction
- `Phase1_WaitResults` — `sem_wait(&customers)` × N per Phase 1 iteration *(inside Phase 1 region)*
- `Phase1_SolveMaster` — GLPK simplex on master LP *(inside Phase 1 region)*
- `Drain` — discarding N stale transition results (service_queue dequeue × N) *(Phase 1→2 transition)*
- `PushPhase2Duals` — writing Phase 2 duals to `sub_data[i]` *(Phase 1→2 transition)*
- `Phase2_WaitResults` — `sem_wait(&customers)` × N per Phase 2 iteration *(inside Phase 2 region)*
- `Phase2_SolveMaster` — GLPK simplex on master LP *(inside Phase 2 region)*
- `FinalSolve` — last master LP solve
- `WriteOutputs` — writing `relaxed_solution`, `zeros_rounded`, `integerized_zeros` *(basis/intermediate LPs only when `--write-bases` / `--write-int-probs` are set)*
- `OptionalRounding` — `dw_rounding` (shown with `[if -r/--round]` guard)

### Subproblem thread states (applies to all N workers)
- `InitialSolve` — GLPK simplex on subproblem LP (initial)
- `WaitMasterReady` — `pthread_cond_wait(&master_lp_ready_cv, ...)`
- `SetupColumns` — build column structures from master
- `SolveLP` — GLPK simplex (or `glp_intopt` if integer) per iteration
- `EnqueueSelf` — lock `service_queue_mutex`, append `my_id`, `sem_post(&customers)`
- `WaitNextIteration` — `pthread_cond_wait(&next_iteration_cv, ...)`

---

## Required regions / annotations

- **Phase 1 region**: All Phase 1 master states MUST be enclosed in a labeled `state "Phase 1 (conditional)" as Phase1` block. A `note` or annotation MUST indicate it runs only when the initial master LP has auxiliary variables.
- **Phase 1→2 transition**: The `Drain` and `PushPhase2Duals` states MUST appear as distinct states between the Phase 1 and Phase 2 regions — not inside either region.
- **Phase 2 region**: All Phase 2 master states MUST be enclosed in a labeled `state "Phase 2 — column generation" as Phase2` block.

---

## Required edge labels (sync primitives)

Every sync primitive must appear at least once as a label on a transition edge or as a `note` adjacent to the state it governs.

| Primitive | On which edge / note |
|-----------|---------------------|
| `master_lp_ready_cv` broadcast | `BuildMasterLP → (master broadcasts master_lp_ready_cv)` |
| `master_lp_ready_cv` wait | `InitialSolve → WaitMasterReady` |
| `customers` sem_post | `EnqueueSelf` note or edge label |
| `customers` sem_wait | `Phase1_WaitResults` or `Phase2_WaitResults` note or edge label |
| `service_queue_mutex` | `EnqueueSelf` note AND `Phase1_WaitResults`/`Phase2_WaitResults` note |
| `next_iteration_cv` broadcast | edge from `Phase1_SolveMaster` → end of Phase 1 iter, and `Phase2_SolveMaster` → end of Phase 2 iter |
| `next_iteration_cv` wait | `WaitNextIteration` note or edge label |
| `sub_data_mutex[i]` | `PushPhase2Duals` note |
| `master_mutex` | `Phase1_SolveMaster` or `Phase2_SolveMaster` note (master writes duals) |

---

## Required prose sections in `threading-flow.md`

The file must contain the following sections in order:

1. A one-paragraph introduction explaining what the diagram shows and who the intended reader is.
2. The Mermaid code block (the diagram itself).
3. **"How to read this diagram"** — a bulleted walkthrough of one complete DW iteration, naming each sync primitive in the order it fires.
4. **"Phase 1→2 transition"** — a short paragraph specifically explaining the drain step and why it is needed (references the race condition fixed in commit `7576722`).
5. **"Single-subproblem variant (N=1)"** — a note explaining what differs when `num_clients = 1` (no fan-out, semaphore posts once).

---

## Acceptance test

A reviewer can correctly answer all of the following by reading the diagram and its prose:

1. What does the master thread wait on while subproblems are solving? (`customers` semaphore)
2. What does a subproblem do after it finishes solving? (enqueues itself, posts `customers`, waits on `next_iteration_cv`)
3. Why is there a Drain step between Phase 1 and Phase 2? (to discard stale results carrying Phase 1 duals)
4. Are Phase 1 states always executed? (No — only when auxiliary variables are detected)
5. Which mutex protects the `service_queue`? (`service_queue_mutex`)
