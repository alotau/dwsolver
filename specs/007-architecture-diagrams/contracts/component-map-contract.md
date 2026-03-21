# Contract: component-map.md

**Diagram type**: `graph LR`
**File**: `architecture/component-map.md`
**Requirements covered**: FR-004, FR-006, FR-007

---

## Required nodes

Every module in this list MUST appear as a separate node.

| Node ID | Label | Type |
|---------|-------|------|
| `dw_main` | `dw_main` | dwsolver module |
| `dw_phases` | `dw_phases` | dwsolver module |
| `dw_subprob` | `dw_subprob` | dwsolver module |
| `dw_support` | `dw_support` | dwsolver module |
| `dw_globals` | `dw_globals\n(shared state)` | dwsolver module |
| `dw_blas` | `dw_blas` | dwsolver module |
| `dw_rounding` | `dw_rounding\n(optional)` | dwsolver module |
| `GLPK` | `GLPK 4.44\n(vendored)` | External block |

The GLPK node MUST be placed inside a `subgraph External` block to visually distinguish it from the seven dwsolver modules.

`dw_globals` MUST be rendered differently from caller/callee modules (e.g., using a different shape such as `[(dw_globals\nshared state)]` cylinder or `{dw_globals}` diamond) to indicate it is a definitions-only data hub, not a functional caller.

---

## Required edges

| From | To | Label |
|------|----|-------|
| `dw_main` | `dw_support` | `init / CLI / output` |
| `dw_main` | `dw_subprob` | `pthread_create` |
| `dw_main` | `dw_phases` | `phase_1_iteration\nphase_2_iteration` |
| `dw_main` | `dw_rounding` | `check_col_integrality\n[if rounding_flag]` |
| `dw_main` | `GLPK` | `master LP ops` |
| `dw_phases` | `dw_subprob` | `reads results` |
| `dw_phases` | `dw_support` | `dw_printf` |
| `dw_phases` | `GLPK` | `dual value reads` |
| `dw_subprob` | `dw_blas` | `sparse matrix ops` |
| `dw_subprob` | `dw_support` | `organize_solution\nprepare_column` |
| `dw_subprob` | `GLPK` | `subproblem LP solve` |
| `dw_rounding` | `dw_blas` | `sparse matrix ops` |
| `dw_rounding` | `dw_support` | `utility functions` |
| `dw_rounding` | `GLPK` | `LP manipulation` |
| `dw_support` | `GLPK` | `parse LP files` |
| all modules | `dw_globals` | `read/write globals` (shown as a single shared annotation, not N edges) |

**Note on `dw_globals` edges**: Rather than drawing 7 edges to `dw_globals` (which would clutter the diagram), use a single `note` or `subgraph` annotation indicating that all modules access globals via `dw.h` extern declarations.

---

## Required visual groupings

Use a `subgraph Modules` block enclosing the 7 dwsolver nodes, and a separate `subgraph External` block for the GLPK node. This makes the boundary between dwsolver code and vendored code immediately visible.

---

## Required prose sections in `component-map.md`

1. A one-paragraph introduction explaining what the diagram shows.
2. The Mermaid code block.
3. **"Module responsibilities"** — a brief (1–2 sentence) description of each node's role, matching the "Primary responsibility" column in `data-model.md`.
4. **"GLPK boundary"** — a sentence explaining what GLPK is, where it lives (`src/glp*.c`), and why it is shown as a single block.
5. **"dw_rounding: optional post-processing"** — a note that `dw_rounding` is not part of the core iterative DW loop.

---

## Acceptance test

A developer can correctly identify which source file to modify for each of these scenarios by reading only the component diagram:

1. "I need to change how reduced costs are computed and sent to subproblems." → `dw_phases`
2. "I need to change how the subproblem LP is loaded from disk." → `dw_subprob` (or `dw_support`)
3. "I need to add a new synchronization primitive." → `dw_globals`
4. "I want to improve the integer rounding heuristic." → `dw_rounding`
