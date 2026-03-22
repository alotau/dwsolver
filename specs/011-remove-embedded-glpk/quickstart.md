# Developer Quickstart: Building dwsolver After Feature 011

**Feature**: 011-remove-embedded-glpk  
**Audience**: Contributors and integrators building from source  
**Date**: 2026-03-22

After this feature lands, dwsolver no longer bundles GLPK source code.
You must install GLPK ≥ 4.65 on your system before running `./configure`.

---

## 1. Install GLPK

### macOS (Homebrew)

```bash
brew install glpk
```

Verify: `pkg-config --modversion glpk`

### Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y libglpk-dev
```

Verify: `pkg-config --modversion glpk`

### Fedora / RHEL / Rocky

```bash
sudo dnf install glpk-devel
```

### Build GLPK from source (any platform)

Use this if your distro ships a version older than 4.65, or if you need
a specific version:

```bash
wget https://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz
tar xzf glpk-5.0.tar.gz
cd glpk-5.0
./configure --prefix=/usr/local
make -j$(nproc)
sudo make install
sudo ldconfig           # Linux only
```

Then set `PKG_CONFIG_PATH` if you installed to a non-standard prefix:

```bash
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH
```

---

## 2. Build dwsolver

```bash
# From the project root
autoreconf -fi          # only needed after a fresh clone or configure.ac change
./configure
make -j$(nproc)
```

### macOS with named semaphores (required for multi-threaded use on macOS)

```bash
./configure --enable-named-semaphores --enable-recursive-mutex
make -j$(nproc)
```

### Confirm GLPK was found

The configure output should contain a line like:

```
checking for glpk >= 4.65... yes
```

If it prints `no` or an error, re-check that `pkg-config` can find GLPK:

```bash
pkg-config --exists 'glpk >= 4.65' && echo OK || echo MISSING
pkg-config --modversion glpk
```

---

## 3. Run the test suite

```bash
make check
```

A successful run produces `PASS` for all tests listed in `tests/`.

---

## 4. Docker build

The Dockerfile handles all dependency installation automatically.  No
system GLPK is required on the Docker host:

```bash
docker build -t dwsolver .
docker run --rm dwsolver --help
```

---

## 5. Linking against `libdwsolver`

Consumers of `libdwsolver` link only against `libdwsolver.so`; they do **not**
need to add `-lglpk` to their own linker flags.  The linker records the
transitive dependency.  Use pkg-config:

```bash
# Compile and link
cc myapp.c $(pkg-config --cflags --libs dwsolver) -o myapp
```

If you use `--static`, pkg-config will include the GLPK static archive
transitively:

```bash
cc myapp.c $(pkg-config --cflags --libs --static dwsolver) -o myapp
```

---

## 6. Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `configure: error: GLPK >= 4.65 is required` | GLPK not installed or too old | Install/upgrade GLPK (§1) |
| `pkg-config: command not found` | pkg-config missing | `sudo apt-get install pkg-config` or `brew install pkg-config` |
| `error while loading shared libraries: libglpk.so.40` | Runtime lib missing | `sudo apt-get install libglpk40` (Ubuntu) or `sudo ldconfig` after installing from source |
| Build fails with `lpx_read_cpxlp undeclared` | Old object files present | `make clean && make` |
| macOS: runtime crash in multi-threaded mode | Named semaphores not enabled | Rebuild with `--enable-named-semaphores --enable-recursive-mutex` |
