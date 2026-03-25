# Security Policy

## Supported Versions

Only the current major release series receives security updates.

| Version | Supported          |
|---------|--------------------|
| 2.x     | :white_check_mark: |
| 1.x and earlier | :x:      |

## Reporting a Vulnerability

**Please do NOT open a public GitHub issue for security vulnerabilities.**

Disclosing a vulnerability publicly before a fix is available puts all users at risk. Instead, report it privately through one of the following channels:

1. **GitHub Private Security Advisory** (preferred):  
   [https://github.com/alotau/dwsolver/security/advisories/new](https://github.com/alotau/dwsolver/security/advisories/new)  
   This requires a GitHub account. The advisory is visible only to you and the maintainers until a fix is released.

2. **Email fallback** (if you cannot use GitHub):  
   Send a detailed report to **joey@littlepancake.com** with the subject line `[SECURITY] dwsolver vulnerability report`.

Please include in your report:
- A clear description of the vulnerability and its potential impact
- Steps to reproduce, or a minimal proof-of-concept input file
- The version of dwsolver and GLPK you are using
- Your operating system and compiler version

## Response Timeline

| Milestone | Target |
|-----------|--------|
| Acknowledgment | Within **14 business days** of report receipt |
| Resolution or coordinated disclosure date agreed | Within **60 days** of acknowledgment |

If a fix requires more than 60 days, the maintainer will notify you with a status update and agree on an extended timeline. We will not disclose the vulnerability publicly before a fix is available, except in cases where the issue is already being actively exploited.

## Scope

### In scope

The following categories constitute security vulnerabilities for this project:

- **Memory safety** — buffer overflows, use-after-free, out-of-bounds reads/writes in solver or library code
- **Integer overflow** — arithmetic overflow in loop indices or matrix dimension calculations that could be triggered by a crafted LP input file
- **Data races** — thread-safety violations in multi-threaded subproblem solving paths
- **GLPK dependency** — vulnerabilities in the GLPK shared library that materially affect dwsolver's security posture (e.g., a crafted `.lp` file exploiting a GLPK parser bug)

### Out of scope

The following are **not** treated as security vulnerabilities:

- **Denial of service via pathological LP inputs** — solving NP-hard linear programs can require exponential time; resource exhaustion from large or degenerate inputs is an inherent mathematical property, not a vulnerability.
- **Build-time tool vulnerabilities** — vulnerabilities in `automake`, `libtool`, or other tools used only to build the software from source.
- **Theoretical issues with no known exploit path** — academic findings without a demonstrated practical impact on dwsolver users.

## Security Design Notes

This section summarises the relevant security characteristics of the codebase to assist security researchers.

- **No cryptographic operations.** dwsolver is a numerical linear programming solver. It performs no encryption, hashing, key derivation, or any other cryptographic function. All OpenSSF Best Practices criteria relating to cryptography are marked "Not applicable" for this project.

- **Input surface.** The solver reads LP problem data from files in CPLEX LP format (`.cpxlp`) or GLPK's native MPS format. Input parsing is delegated entirely to the GLPK library. There is no shell command construction, no SQL, and no network I/O — the attack surface is limited to local file input and the GLPK shared library's parser.

- **Integer overflow risk.** Loop indices and matrix dimension arithmetic use C `int`. Very large LP instances (approaching `INT_MAX` rows or columns) could theoretically trigger signed integer overflow. This is documented as a known limitation; users working with extremely large models should validate their input dimensions.

- **Thread safety model.** Subproblem solving uses POSIX threads (`pthreads`). All shared solver state is protected by `pthread_mutex_t` mutexes. The public library API (`dw_solve()`) is **not** safe for concurrent calls on the same context object from multiple threads; each concurrent solver invocation must use a separate context.

- **No network access.** dwsolver makes no network connections at any point during initialisation, solving, or output. There is no telemetry, no update check, and no remote data retrieval.

## Disclosure Policy

1. Reporter submits a private Security Advisory or email.
2. Maintainer acknowledges the report within 14 business days.
3. A fix is developed on a private branch. The reporter is kept informed of progress.
4. If the severity warrants, a CVE is requested via the GitHub Security Advisory workflow.
5. Once the fix is ready and tagged, the Security Advisory is published, the `CHANGELOG.md` `### Security` section is updated, and the reporter is credited (with their permission).
6. We follow coordinated disclosure — we will not publish details before a fix is available.
