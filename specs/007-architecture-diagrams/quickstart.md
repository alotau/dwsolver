# Quickstart: Maintaining Architecture Diagrams

**Feature**: 007-architecture-diagrams
**Date**: 2026-03-21

---

## Where the diagrams live

```
architecture/
├── README.md            — index with one-sentence descriptions
├── threading-flow.md    — stateDiagram-v2: master + subproblem thread states
├── component-map.md     — graph LR: module call/dependency graph
└── data-flow.md         — flowchart TD: input → process → output
```

## When to update a diagram

| If you change... | Update this diagram |
|------------------|-------------------|
| Any mutex, semaphore, or condition variable | `threading-flow.md` |
| The Phase 1→2 transition logic | `threading-flow.md` |
| Threading startup or teardown | `threading-flow.md` |
| Which module calls which function | `component-map.md` |
| Adding or renaming a source module | `component-map.md` |
| Guidefile format (field order, new fields) | `data-flow.md` |
| New input or output files | `data-flow.md` |
| New Phase 1 or Phase 2 termination path | `data-flow.md` + `threading-flow.md` |

**Convention**: A PR that touches threading or module structure MUST include a diagram update in the same commit. This is part of the definition of done (spec assumption 5).

---

## How to preview locally

GitHub renders Mermaid natively when you push. To preview before pushing:

1. **VS Code**: Install the "Markdown Preview Mermaid Support" extension, then open the `.md` file and press `Cmd+Shift+V`.
2. **mermaid.live**: Copy the fenced code block content to https://mermaid.live for instant rendering.
3. **GitHub web editor**: Edit the file directly in the GitHub UI — the preview tab renders Mermaid.

No build step, no local install required.

---

## Mermaid syntax reference (project-specific)

### stateDiagram-v2 (threading-flow.md)
```
stateDiagram-v2
    [*] --> Init
    Init --> SpawnSubproblems
    note right of SpawnSubproblems: pthread_create × N

    state "Phase 1 (conditional)" as Phase1 {
        WaitN_Ph1 --> SolveLP_Ph1 : next_iteration_cv broadcast
    }
```

### graph LR (component-map.md)
```
graph LR
    subgraph External
        GLPK["GLPK 4.44 (vendored)"]
    end
    dw_main -->|phase_1_iteration| dw_phases
    dw_main -->|pthread_create| dw_subprob
```

### flowchart TD (data-flow.md)
```
flowchart TD
    GF[guidefile] --> PC[process_cmdline\ndw_support]
    PC --> N["num_clients = N\nsubproblem filenames\nmaster filename"]
```

---

## Updating the index

After adding, removing, or renaming a diagram, update `architecture/README.md` to keep the index accurate. The README requires exactly one line per diagram: its filename and a single sentence describing what it shows.
