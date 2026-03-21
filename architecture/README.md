# Architecture

This directory contains architecture diagrams for the DWSOLVER repository.
All diagrams are authored in [Mermaid](https://mermaid.js.org/) syntax and render natively in GitHub Markdown — no plugins or external viewers required.

## Diagrams

| Diagram | Description |
|---------|-------------|
| [threading-flow.md](./threading-flow.md) | Threading and control-flow diagram: master thread, N subproblem threads, all synchronization primitives, Phase 1/2 regions, and the Phase 1→2 drain transition |
| [component-map.md](./component-map.md) | Component diagram: the seven dwsolver source modules, their call and dependency relationships, and the vendored GLPK library as a single external block |
| [data-flow.md](./data-flow.md) | Input data-flow diagram: guidefile format, LP input files, every major transformation stage (Phase 1, Phase 2, rounding), and all output files |

## Keeping diagrams in sync

See [specs/007-architecture-diagrams/quickstart.md](../specs/007-architecture-diagrams/quickstart.md) for a table of "if you change X, update diagram Y" maintenance rules.
Any PR that alters threading logic, module structure, or the input/output file format must include a corresponding diagram update.
