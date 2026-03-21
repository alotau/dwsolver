# Contract: data-flow.md

**Diagram type**: `flowchart TD`
**File**: `architecture/data-flow.md`
**Requirements covered**: FR-005, FR-006, FR-007

---

## Required nodes

| Node ID | Label | Tier |
|---------|-------|------|
| `GF` | `guidefile\n(-g flag)` | Input |
| `MLP` | `master.lp\n(GLPK format)` | Input |
| `SLP` | `subproblem_i.lp\n(GLPK format, × N)` | Input |
| `PARSE` | `process_cmdline\n(dw_support)\n→ num_clients, filenames` | Transform |
| `LOAD` | `Load LPs via GLPK\n(dw_main + dw_subprob)` | Transform |
| `INIT_SOLVE` | `Initial subproblem solves\n(GLPK simplex, × N parallel)` | Transform |
| `BUILD_MASTER` | `Build reduced master LP\n(dw_main + GLPK)` | Transform |
| `PH1` | `Phase 1 iterations\nphase_1_iteration()\n[conditional]` | Transform (conditional) |
| `DRAIN` | `Phase 1→2 drain\nDiscard stale results\nPush Phase 2 duals` | Transform |
| `PH2` | `Phase 2 iterations\nphase_2_iteration()\ncolumn generation loop` | Transform |
| `FINAL` | `Final master LP solve\n(GLPK simplex)` | Transform |
| `ROUND` | `Rounding\n(dw_rounding)\n[if rounding_flag]` | Transform (optional) |
| `OUT_RS` | `relaxed_solution` | Output |
| `OUT_ZV` | `zero_vars` | Output |
| `OUT_INT` | `integer_solution\n[if rounding_flag]` | Output (conditional) |

---

## Required edges

```
GF --> PARSE
MLP --> LOAD
SLP --> LOAD
PARSE --> LOAD
LOAD --> INIT_SOLVE
INIT_SOLVE --> BUILD_MASTER
BUILD_MASTER --> PH1
PH1 --> DRAIN
DRAIN --> PH2
BUILD_MASTER --> PH2  (dashed / annotated: "if no aux vars — skip Phase 1")
PH2 --> FINAL
FINAL --> ROUND
FINAL --> OUT_RS
FINAL --> OUT_ZV
ROUND --> OUT_INT
```

The `BUILD_MASTER → PH2` bypass edge MUST be shown as a dashed line (`-.->`) with the annotation `no auxiliary variables`.

The `PH1` node MUST carry a visual indicator that it is conditional (use a note or an `[if aux vars needed]` label on the edge from `BUILD_MASTER`).

The `ROUND` node MUST carry a `[if rounding_flag]` label.

---

## Required subgraph groupings

| Subgraph | Nodes inside |
|----------|-------------|
| `Input files` | `GF`, `MLP`, `SLP` |
| `Startup` | `PARSE`, `LOAD`, `INIT_SOLVE`, `BUILD_MASTER` |
| `Phase 1 (conditional)` | `PH1`, `DRAIN` |
| `Phase 2 — column generation` | `PH2` |
| `Finalization` | `FINAL`, `ROUND` |
| `Outputs` | `OUT_RS`, `OUT_ZV`, `OUT_INT` |

---

## Required guidefile format annotation

Somewhere in the diagram (as a `note` on the `GF` node, or as a separate text block below the diagram) the guidefile format MUST be shown:

```
line 1:   <N>              ← number of subproblems
lines 2..N+1: <sp_i.lp>   ← subproblem LP filenames
line N+2: <master.lp>      ← master LP filename
```

---

## Required prose sections in `data-flow.md`

1. A one-paragraph introduction explaining what inputs the solver needs and what it produces.
2. The Mermaid code block.
3. **"Guidefile format"** — the exact line-by-line format (verbatim from the annotation above), with an example for a two-subproblem case.
4. **"Phase 1: when does it run?"** — a short explanation that Phase 1 only activates when auxiliary variables are introduced by infeasibility in the initial LP.
5. **"Output files"** — a table mapping each output filename to the condition under which it is produced.

---

## Acceptance test

A first-time user can correctly answer all of the following by reading the diagram and its prose:

1. How many lines does the guidefile have for a problem with 3 subproblems? (5 lines: 1 count + 3 subproblem names + 1 master name)
2. Do subproblem filenames come before or after the master filename in the guidefile? (before)
3. Under what condition is `integer_solution` produced? (when `--rounding_flag` is passed)
4. Is Phase 1 always executed? (No — only if the initial master LP requires auxiliary variables)
5. What file format are the LP inputs? (GLPK LP format)
