# Threading and Control-Flow Diagram

This diagram shows the complete lifecycle of the master thread and the N subproblem
threads in a single DWSOLVER run — from startup through Phase 1 (if needed), the
Phase 1→2 drain transition, Phase 2 column generation, and final output.
It is the primary reference for anyone modifying concurrency logic: every mutex,
semaphore, and condition variable that gates a state transition is labeled on the
relevant edge or state.

```mermaid
stateDiagram-v2

    [*] --> fork_start
    state fork_start <<fork>>
    fork_start --> Master
    fork_start --> SubN

    %% ── MASTER THREAD ─────────────────────────────────────────────────────────
    state "Master Thread" as Master {
        [*] --> Init
        Init --> SpawnSubproblems
        SpawnSubproblems --> BuildMasterLP

        note right of SpawnSubproblems
            pthread_create(&threads[i], subproblem_thread) × N
        end note

        note right of BuildMasterLP
            Broadcasts master_lp_ready_cv when master LP is ready.
            Subproblem threads waiting on master_lp_ready_cv are released.
        end note

        BuildMasterLP --> Phase1         : aux vars detected
        BuildMasterLP --> Phase2         : no auxiliary variables

        %% ── PHASE 1 (conditional) ────────────────────────────────────────────
        state "Phase 1 (conditional — only when auxiliary variables detected)" as Phase1 {
            [*] --> Ph1_WaitN

            Ph1_WaitN --> Ph1_Solve      : N results collected
            Ph1_Solve --> Ph1_WaitN      : not converged — broadcast next_iteration_cv
            Ph1_Solve --> [*]            : Phase 1 converged

            note right of Ph1_WaitN
                sem_wait(customers) × N
                Lock service_queue_mutex → read service_queue[head] → advance head → unlock
            end note

            note right of Ph1_Solve
                Solve master LP (glp_simplex)
                Extract duals (master_mutex)
                Push duals to sub_data[i] (sub_data_mutex[i])
                Increment current_iteration
                Broadcast next_iteration_cv
            end note
        }

        Phase1 --> Drain

        %% ── PHASE 1 → 2 TRANSITION ───────────────────────────────────────────
        Drain --> PushPhase2Duals

        note right of Drain
            Discard N in-flight stale results from Phase 1:
            sem_wait(customers) × N
            Lock service_queue_mutex → advance head × N → unlock
            (results NOT used — they carry Phase 1 duals)
        end note

        note right of PushPhase2Duals
            Lock master_mutex → extract Phase 2 row_duals → unlock
            For each i: lock sub_data_mutex[i] → write sub_data[i].r → unlock
            Increment current_iteration
            Broadcast next_iteration_cv
            (subproblems wake with phase_one=0 and Phase 2 duals)
        end note

        PushPhase2Duals --> Phase2

        %% ── PHASE 2 ──────────────────────────────────────────────────────────
        state "Phase 2 — Column Generation" as Phase2 {
            [*] --> Ph2_WaitN

            Ph2_WaitN --> Ph2_Solve      : N results collected
            Ph2_Solve --> Ph2_WaitN      : continue — broadcast next_iteration_cv
            Ph2_Solve --> [*]            : converged or max iterations reached

            note right of Ph2_WaitN
                sem_wait(customers) × N
                Lock service_queue_mutex → read service_queue[head] → advance head → unlock
            end note

            note right of Ph2_Solve
                Solve master LP (glp_simplex)
                Extract duals (master_mutex)
                Push duals to sub_data[i] (sub_data_mutex[i])
                Increment current_iteration
                Broadcast next_iteration_cv
            end note
        }

        Phase2 --> FinalSolve
        FinalSolve --> OptRounding       : rounding_flag set (--rounding_flag)
        FinalSolve --> WriteOutputs
        OptRounding --> WriteOutputs
        WriteOutputs --> [*]
    }

    %% ── SUBPROBLEM THREAD × N ────────────────────────────────────────────────
    state "Subproblem Thread (× N, one per subproblem LP)" as SubN {
        [*] --> InitialSolve
        InitialSolve --> WaitMasterReady

        note right of WaitMasterReady
            Lock master_lp_ready_mutex
            If master_lp_ready == 0: pthread_cond_wait(master_lp_ready_cv)
            Unlock master_lp_ready_mutex
        end note

        WaitMasterReady --> SetupColumns
        SetupColumns --> SolveLP

        SolveLP --> EnqueueSelf

        note right of EnqueueSelf
            Lock service_queue_mutex
            service_queue[tail] = my_id; advance tail
            sem_post(customers)
            Unlock service_queue_mutex
        end note

        EnqueueSelf --> WaitNextIter

        note right of WaitNextIter
            Lock next_iteration_mutex
            While current_iteration unchanged:
              pthread_cond_wait(next_iteration_cv, next_iteration_mutex)
            Unlock next_iteration_mutex
        end note

        WaitNextIter --> SolveLP         : next iteration (phase_one flag checked inside SolveLP)
        WaitNextIter --> [*]             : solver signals done
    }

    state join_end <<join>>
    Master --> join_end
    SubN --> join_end
    join_end --> [*]
```

---

## How to read this diagram

Follow one complete DW iteration from the master's perspective, noting each sync primitive in the order it fires:

1. **Master waits for N results** — enters `Ph1_WaitN` or `Ph2_WaitN` and calls `sem_wait(customers)` once per subproblem. The `customers` semaphore is the signal that a subproblem has finished solving and enqueued itself.

2. **Subproblem posts a result** — after `SolveLP` completes, the subproblem enters `EnqueueSelf`: it locks `service_queue_mutex`, appends its own ID to `service_queue[tail]`, advances the tail pointer, calls `sem_post(customers)` (decrementing the master's wait count by one), and unlocks.

3. **Master dequeues and processes** — on each `sem_wait` return, the master locks `service_queue_mutex`, reads `service_queue[head]`, advances the head pointer, and unlocks. This gives it the index of the subproblem that reported. Master pushes a new column to the reduced master LP for that result.

4. **Master solves and pushes duals** — after collecting all N results, master calls `glp_simplex` on the master LP (protected by `master_mutex` for dual reads), then iterates over `sub_data[i]`, writing updated `r` values under `sub_data_mutex[i]`.

5. **Master wakes subproblems** — master increments `current_iteration` and calls `pthread_cond_broadcast(next_iteration_cv)` under `next_iteration_mutex`. All subproblems waiting in `WaitNextIter` unblock, re-lock, check the iteration counter, and return to `SolveLP`.

6. **Repeat** — steps 1–5 repeat until reduced costs satisfy the optimality tolerance or the maximum iteration count is reached.

---

## Phase 1→2 transition

When Phase 1 converges, the last `phase_1_iteration()` call has already broadcast `next_iteration_cv` to release subproblems. Those subproblems wake up, solve with `phase_one=1` and Phase 1 duals, and enqueue their results — but master is no longer in Phase 1's wait loop. It has moved on to clearing auxiliary variables and setting `sub_data[j].phase_one = 0` for all j.

**The race**: if master immediately entered `phase_2_iteration()`, subproblem results still carrying Phase 1 duals and the Phase 1 convexity objective would be read as Phase 2 results. The reduced-cost check would use the Phase 1 convexity-constraint dual (`r_phase1`), which can be ≤ 0 at Phase 1 convergence, causing premature termination at a suboptimal point.

**The drain fix** (introduced in commit `7576722`): master explicitly discards all N in-flight transition results by calling `sem_wait(customers)` × N and dequeuing from `service_queue` without using the results. It then pushes fresh Phase 2 duals from the just-solved master LP into every `sub_data[i]`, increments `current_iteration`, and rebroadcasts `next_iteration_cv`. Subproblems wake with `phase_one=0` and correct Phase 2 duals, and their next iteration is the first genuine Phase 2 solve.

---

## Single-subproblem variant (N=1)

When `num_clients = 1`, the fan-out is trivial: there is exactly one subproblem thread, and the master's `sem_wait(customers)` loop runs only once per iteration. The `service_queue` still exists but always holds at most one entry. All synchronization primitives are still used and the diagram remains accurate; the only behavioral difference is that "N results collected" means "1 result collected".
