# Feature Specification: Architecture Diagrams Directory

**Feature Branch**: `007-architecture-diagrams`
**Created**: 2026-03-21
**Status**: Draft
**Input**: "Add architecture diagrams directory with flow and component diagrams"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Understand solver threading and control flow (Priority: P1)

A new contributor or reviewer opens the repository and wants to understand how the master thread and subproblem threads communicate — which synchronization primitives gate each transition, how the service queue works, and what constitutes one full Dantzig-Wolfe iteration.

**Why this priority**: This is the single hardest part of the codebase to grasp from source alone. The Phase 1→2 transition race that was recently fixed is a direct consequence of this flow being undocumented. Any future contributor touching threading must understand this first.

**Independent Test**: A reviewer who has not read the source can correctly answer "what happens after `pthread_cond_broadcast` fires at the end of a Phase 2 iteration?" solely from the flow diagram.

**Acceptance Scenarios**:

1. **Given** a contributor reading the flow diagram, **When** they trace the path from master broadcasting to a subproblem re-solving, **Then** every mutex, semaphore, and condition variable on that path is labeled.
2. **Given** the Phase 1→2 drain pattern now in `dw_main.c`, **When** a reviewer reads the diagram, **Then** the drain step is shown as a distinct state between Phase 1 completion and the first Phase 2 iteration.
3. **Given** a new contributor, **When** they look at any state box in the diagram, **Then** the owning thread (master or subproblem-N) is identifiable.

---

### User Story 2 - Understand module responsibilities and dependencies (Priority: P2)

A developer adding a new feature wants to know which source file to modify without reading every `.c` file. They need a component view showing the seven dwsolver modules and their relationships to each other and to the vendored GLPK library.

**Why this priority**: The component map is useful for onboarding and for scoping change impact, but less urgent than the flow diagram since module names are already visible in the directory listing.

**Independent Test**: A developer can determine "which file handles building the master LP column from a subproblem solution" by reading only the component diagram.

**Acceptance Scenarios**:

1. **Given** the component diagram, **When** a developer inspects `dw_phases`, **Then** they can see it depends on shared synchronization state from `dw_globals` and on GLPK for LP operations.
2. **Given** the component diagram, **When** a reviewer inspects the GLPK dependency, **Then** it is shown as a single external block distinct from the seven dwsolver modules.
3. **Given** the component diagram, **When** a developer reads `dw_rounding`, **Then** they can see it is an optional post-processing step, not part of the core iterative loop.

---

### User Story 3 - Understand the input file format and data flow (Priority: P3)

A user setting up a new Dantzig-Wolfe problem wants to understand what the guidefile, master LP file, and subproblem LP files are, the order they must appear, and how they are consumed to produce the output files.

**Why this priority**: This serves a wider audience (end-users as well as developers) but is lower priority than the code-facing diagrams because the README partially covers the input format.

**Independent Test**: A first-time user can prepare a valid guidefile and correctly ordered LP files for a two-subproblem problem by reading only this diagram alongside the README.

**Acceptance Scenarios**:

1. **Given** the data-flow diagram, **When** a user reads it, **Then** they understand that the guidefile line count controls the number of subproblems and that subproblem filenames are listed before the master filename.
2. **Given** the data-flow diagram, **When** a user traces the path from input files to `relaxed_solution`, **Then** every major transformation (initial solve, Phase 1 column generation, Phase 2 column generation, final master solve) is represented.

---

### Edge Cases

- Diagrams must remain accurate when `num_clients = 1` — the single-subproblem path omits some synchronization steps; the diagram notes this variant or an annotation explains it.
- The Phase 1 section of the flow diagram must be clearly marked as conditionally executed (only when auxiliary variables are needed), not always present.
- The diagram format must render in GitHub Markdown without any browser extensions or third-party viewers.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The repository MUST contain an `architecture/` directory at the project root.
- **FR-002**: The directory MUST contain a threading and control-flow diagram covering the master thread, N subproblem threads, and all synchronization primitives (mutexes, semaphores, condition variables) across one full DW iteration.
- **FR-003**: The flow diagram MUST distinguish Phase 1 (auxiliary variable elimination) and Phase 2 (column generation) as separately labeled regions, with the Phase 1→2 drain step explicitly shown as a distinct state.
- **FR-004**: The directory MUST contain a component diagram showing the seven dwsolver source modules (`dw_main`, `dw_phases`, `dw_subprob`, `dw_support`, `dw_globals`, `dw_blas`, `dw_rounding`), their call/dependency relationships, and GLPK as a single external block.
- **FR-005**: The directory MUST contain a data-flow diagram showing how the guidefile, master LP, and subproblem LP files are consumed and how results flow to the output files.
- **FR-006**: All diagrams MUST be authored in Mermaid syntax so they render natively in GitHub Markdown without plugins or external hosting.
- **FR-007**: Each diagram MUST be accompanied by prose explaining what it shows and what a reader should take away from it.
- **FR-008**: An `architecture/README.md` MUST list all diagrams with a one-sentence description of each, serving as an index.
- **FR-009**: The project root `README.md` MUST be updated to link to the `architecture/` directory so it is discoverable from the repository landing page.

### Key Entities

- **Diagram**: A Mermaid-sourced visual artifact stored as a `.md` file and renderable on GitHub without plugins.
- **`architecture/` directory**: The single versioned location for all architecture documentation, kept at the project root.
- **Synchronization primitive**: Any mutex, semaphore, or condition variable that gates a thread-state transition; each must appear labeled in the flow diagram.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A reviewer unfamiliar with the codebase can correctly describe the master-subproblem handshake (`sem_post` → `sem_wait` → service queue dequeue → broadcast → `pthread_cond_wait`) after reading only the flow diagram — verifiable during PR review.
- **SC-002**: All diagrams render without errors or broken blocks in a GitHub Markdown preview.
- **SC-003**: The `architecture/` directory contains at minimum three diagrams: threading flow, component map, and input data flow.
- **SC-004**: The root `README.md` links to the architecture directory; a user on the GitHub repository landing page can reach the architecture index in one click.

## Assumptions

- Mermaid is the correct format: GitHub has natively rendered Mermaid in Markdown since 2022 with no viewer action required.
- The `architecture/` directory lives at the project root (not inside `docs/` or `specs/`) so it is immediately visible when browsing the repository.
- The vendored GLPK source files in `src/` (`glp*.c`) are treated as a single external block in the component diagram rather than individual nodes.
- Diagrams are plain text files versioned with the source code; no binary image files are committed.
- Diagrams will be kept in sync with code changes as part of the definition of done for any future PR that alters the threading or module structure.
