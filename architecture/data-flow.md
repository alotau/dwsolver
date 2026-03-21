# Input Data Flow Diagram

This diagram shows what input files DWSOLVER requires, how they are consumed, every major
transformation stage (startup, optional Phase 1, Phase 2 column generation, optional rounding),
and what output files are produced.
It is the primary reference for first-time users setting up a problem and for developers
tracing where a particular file is read or written.

```mermaid
flowchart TD
    subgraph Inputs["Input Files"]
        GF["guidefile\n(-g flag)"]
        SLP["subproblem_i.lp\n(GLPK LP format, × N)"]
        MLP["master.lp\n(GLPK LP format)"]
    end

    subgraph Startup["Startup"]
        PARSE["process_cmdline — dw_support\nReads guidefile:\n  line 1 → num_clients = N\n  lines 2..N+1 → subproblem filenames\n  line N+2 → master filename"]
        LOAD["Load LPs via GLPK\ndw_main + dw_subprob"]
        INIT_SOLVE["Initial subproblem solves\nGLPK simplex, × N parallel threads"]
        BUILD_MASTER["Build reduced master LP\ndw_main + GLPK\nAdd convexity constraints (one per subproblem)"]
    end

    subgraph Ph1Region["Phase 1 — Auxiliary Variable Elimination (conditional)"]
        PH1["phase_1_iteration()\nColumn generation to drive\nauxiliary variables to zero"]
        DRAIN["Phase 1→2 Drain\nDiscard N stale results\nPush Phase 2 duals\nRebroadcast next_iteration_cv"]
    end

    subgraph Ph2Region["Phase 2 — Column Generation"]
        PH2["phase_2_iteration()\nMain Dantzig-Wolfe column generation loop\nIterates until reduced cost ≤ tolerance\nor max_phase2_iterations reached"]
    end

    subgraph Final["Finalization"]
        FINAL["Final master LP solve\nGLPK simplex on fully-columned master LP"]
        ROUND["Rounding — dw_rounding\n(only if --round / -r)"]
    end

    subgraph Outputs["Output Files"]
        OUT_RS["relaxed_solution\n(always)"]
        OUT_ZV["zero_vars\n(always)"]
        OUT_INT["integer_solution\n(only if --round / -r)"]
    end

    GF --> PARSE
    SLP --> LOAD
    MLP --> LOAD
    PARSE --> LOAD
    LOAD --> INIT_SOLVE
    INIT_SOLVE --> BUILD_MASTER

    BUILD_MASTER -->|"auxiliary variables detected"| PH1
    BUILD_MASTER -.->|"no auxiliary variables — skip Phase 1"| PH2

    PH1 --> DRAIN
    DRAIN --> PH2

    PH2 --> FINAL
    FINAL --> ROUND
    FINAL --> OUT_RS
    FINAL --> OUT_ZV
    ROUND --> OUT_INT
```

> **Solid arrows** = always-executed path.
> **Dashed arrow** = bypass path when no auxiliary variables are detected (Phase 1 skipped entirely).

---

## Guidefile format

The guidefile is a plain-text file, one token per line:

```
line 1:     <N>                ← integer: number of subproblems
line 2:     <subproblem_1.lp>  ← filename of the first subproblem LP
line 3:     <subproblem_2.lp>  ← filename of the second subproblem LP
...
line N+1:   <subproblem_N.lp>  ← filename of the Nth subproblem LP
line N+2:   <master.lp>        ← filename of the master LP
```

**Important**: subproblem filenames come *before* the master filename. Getting this order wrong will cause DWSOLVER to load the wrong LP as the master problem.

**Example** — a two-subproblem problem:
```
2
sub1.lp
sub2.lp
master.lp
```

All filenames are resolved relative to the working directory at the time `dwsolver` is invoked. All LP files must be in [GLPK LP format](https://en.wikibooks.org/wiki/GLPK/LP_Format).

---

## Phase 1: when does it run?

Phase 1 runs only when the initial master LP is infeasible with the starting basis — which occurs when the LP requires auxiliary (artificial) variables to reach a feasible point. DWSOLVER detects this by checking the master LP column structure after construction: if any auxiliary variables (`y_i`) have been introduced to handle constraints that cannot be satisfied by the initial subproblem extreme points alone, `need_phase_one` is set and the Phase 1 loop runs.

If the initial LP is already feasible (no auxiliary variables), DWSOLVER skips Phase 1 entirely and proceeds directly to Phase 2. The bypass dashed arrow in the diagram above shows this path.

---

## Output files

| Output file | Produced | Contents |
|-------------|----------|----------|
| `relaxed_solution` | Always | Continuous optimal variable values from the final master LP solve |
| `zeros_rounded` | Only if `-r` / `--round` passed | Solution with selected zero-valued variables rounded according to the rounding heuristic |
| `integerized_zeros` | Only if `-i` / `--integerize` passed | Integer-adjusted variant of `zeros_rounded`, enforcing integrality on chosen variables |
| `basis_iteration_*` | Only if `--write-bases` passed | Per-iteration LP basis snapshots from the master solve |
| `phase1_step_*.cpxlp`, `master_step_*.cpxlp` | Only if `--write-int-probs` passed | CPLEX LP snapshots of the master problem during Phase 1 (`phase1_step_*.cpxlp`) and Phase 2 (`master_step_*.cpxlp`) |

All output files are written to the working directory.
