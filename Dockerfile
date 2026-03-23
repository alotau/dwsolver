# syntax=docker/dockerfile:1

# ---------------------------------------------------------------------------
# Stage 1 – build
# ---------------------------------------------------------------------------
FROM ubuntu:24.04 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        automake \
        pkg-config \
        libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY . .

# Touch pre-generated Autotools files so make doesn't try to re-run automake
# (automake-1.17 required by configure is not available on Ubuntu 24.04).
# Also refresh config.guess/config.sub from the installed automake so the
# build works on aarch64 (Apple Silicon) as well as x86_64.
RUN touch aclocal.m4 configure config.h.in build-aux/ltmain.sh \
         Makefile.in src/Makefile.in \
    && find /usr/share/automake-* -name config.guess | head -1 | xargs -I{} cp {} build-aux/config.guess \
    && find /usr/share/automake-* -name config.sub   | head -1 | xargs -I{} cp {} build-aux/config.sub

# --disable-shared produces a real linked binary at src/dwsolver rather than a
# libtool wrapper script, so COPY in the runner stage works correctly.
# libglpk-dev on Ubuntu does not ship a glpk.pc file, so bypass pkg-config by
# exporting the CFLAGS/LIBS variables that PKG_CHECK_MODULES checks first.
RUN GLPK_CFLAGS=-I/usr/include GLPK_LIBS=-lglpk ./configure --disable-shared && make

# ---------------------------------------------------------------------------
# Stage 2 – runtime (binary only)
# ---------------------------------------------------------------------------
FROM ubuntu:24.04 AS runner

RUN apt-get update && apt-get install -y --no-install-recommends libglpk40 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data
COPY --from=builder /build/src/dwsolver /usr/local/bin/dwsolver

ENTRYPOINT ["/usr/local/bin/dwsolver"]
