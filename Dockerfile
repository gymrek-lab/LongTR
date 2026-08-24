# syntax=docker/dockerfile:1
#
# Portable, fully-static LongTR, built against musl libc on Alpine Linux.
#
# The resulting binary has NO dynamic dependencies and runs on any x86-64 Linux
# (CPU with SSE4.1 or newer). Network / remote-file support (libcurl, OpenSSL, S3,
# GCS) is intentionally disabled via the Makefile's PORTABLE=1 path -- that is what
# makes a clean, universal static link possible. Read local files only.
#
# Build the distributable tarball into ./dist :
#     DOCKER_BUILDKIT=1 docker build \
#         --build-arg VERSION="$(git describe --tags --always)" \
#         --target artifact --output type=local,dest=./dist .
#
# Or build the full image and smoke-test the binary:
#     docker build --target build -t longtr:portable .
#     docker run --rm longtr:portable ./LongTR --help
#
# NOTE: htslib and spoa are cloned fresh from upstream at build time (their dirs are
# excluded via .dockerignore). If an upstream change ever breaks the build, pin them
# to known-good tags in the Makefile's HTSLIB/SPOA targets.
#
ARG ALPINE_VERSION=3.20

# ---------------------------------------------------------------------------
FROM alpine:${ALPINE_VERSION} AS build

# Build toolchain + STATIC system libraries. No libcurl/openssl here: htslib is
# configured with --disable-libcurl under PORTABLE=1, so no network stack is linked.
RUN apk add --no-cache \
        build-base cmake git autoconf automake libtool bash perl linux-headers file \
        zlib-dev zlib-static \
        bzip2-dev bzip2-static \
        xz-dev xz-static \
        libdeflate-dev libdeflate-static

WORKDIR /build
COPY . /build

# Build the vendored deps (htslib without network, spoa with portable SSE4.1) and
# then link LongTR + DenovoFinder fully statically. Two make invocations so the
# dependency archives are guaranteed built before the final link (the Makefile does
# not declare that edge explicitly).
RUN make PORTABLE=1 HTSLIB-docker SPOA-docker \
 && make -j"$(nproc)" PORTABLE=1 LongTR DenovoFinder

# Assert the binary really is static; fail the build otherwise.
RUN file ./LongTR \
 && if ldd ./LongTR 2>/dev/null | grep -q '=>'; then \
        echo "ERROR: LongTR is dynamically linked" >&2; ldd ./LongTR >&2; exit 1; \
    fi \
 && cp LongTR LongTR-debug && cp DenovoFinder DenovoFinder-debug \
 && strip LongTR DenovoFinder

# Package a versioned tarball under /dist.
ARG VERSION=dev
RUN DST="LongTR-${VERSION}-portable-linux-x86_64" \
 && mkdir -p "/dist/${DST}/scripts" \
 && cp LongTR DenovoFinder VizAln VizAlnPdf README.md "/dist/${DST}/" \
 && cp scripts/filter_haploid_vcf.py scripts/filter_vcf.py \
       scripts/generate_aln_html.py scripts/html_alns_to_pdf.py "/dist/${DST}/scripts/" \
 && tar -czf "/dist/${DST}.tar.gz" -C /dist "${DST}" \
 && rm -rf "/dist/${DST}" \
 && ls -l /dist

# ---------------------------------------------------------------------------
# Minimal stage holding only the artifacts, so `docker build --output` extracts
# just the tarball (+ the raw static binary) with no image layers.
FROM scratch AS artifact
COPY --from=build /dist/ /
COPY --from=build /build/LongTR /LongTR
COPY --from=build /build/LongTR-debug /LongTR-debug
