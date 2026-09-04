#!/bin/sh
# Build the standalone C benchmarks of the gsDesign integration routines.
#
# Baseline sources are the pristine gsDesign 3.11.0 files (deps-src), compiled
# with their symbols renamed to base_* so that they can be linked next to the
# patched versions (v2/ = patch 0004 sources, patch3/ = patch 0003 sources).
# Requires R (libR and headers) and clang/gcc. Usage: ./build.sh [--run]
set -e
cd "$(dirname "$0")"
SRC=${GSDESIGN_SRC:-../../deps-src/gsDesign/src}
RHOME=$(R RHOME)
RINC="$RHOME/include"
RLIB="$RHOME/lib"
CC=${CC:-clang}
CFLAGS=${CFLAGS:--O2 -falign-functions=64}

mkdir -p base_renamed
for f in gridpts h1 hupdate probpos rprobrej gsbound gsbound1; do
  perl -pe 's/\bprobrej\b/base_probrej/g; s/\bgsbound1\b/base_gsbound1/g; s/\bgsbound\b/base_gsbound/g;
            s/\bgridpts1\b/base_gridpts1/g; s/\bgridpts\b/base_gridpts/g; s/\bh1\b/base_h1/g;
            s/\bhupdate\b/base_hupdate/g; s/\bprobneg\b/base_probneg/g; s/\bprobpos\b/base_probpos/g;
            s/#include "gsDesign.h"/#include "..\/base_gsDesign.h"/' "$SRC/$f.c" > "base_renamed/$f.c"
done
perl -pe 's/\bgsbound1\b/base_gsbound1/g; s/\bgsbound\b/base_gsbound/g; s/\bprobrej\b/base_probrej/g;
          s/\bgsdensity\b/base_gsdensity/g; s/\bstdnorpts\b/base_stdnorpts/g' "$SRC/gsDesign.h" > base_gsDesign.h

LINK="-L$RLIB -lR -Wl,-rpath,$RLIB"
# pristine routines only (primitive costs, grid sizes, per-routine timings)
$CC $CFLAGS -I"$RINC" -I"$SRC" bench_baseline.c $SRC/gridpts.c $SRC/h1.c $SRC/hupdate.c $SRC/probpos.c \
  $SRC/rprobrej.c $SRC/gsbound.c $SRC/gsbound1.c $LINK -o bench_baseline
# Newton iteration counts (instrumented copies of the pristine search routines)
$CC $CFLAGS -I"$RINC" -I"$SRC" bench_iters.c $SRC/gridpts.c $SRC/h1.c $SRC/hupdate.c $SRC/probpos.c \
  instr/gsbound.c instr/gsbound1.c $LINK -o bench_iters
# patch 0003 (same grid) versus baseline
$CC $CFLAGS -I"$RINC" -Ipatch3 bench_p3.c patch3/hupdate.c patch3/probpos.c patch3/rprobrej.c \
  patch3/gsbound1.c patch3/gsbound.c $SRC/gridpts.c $SRC/h1.c base_renamed/*.c $LINK -o bench_p3
# patch 0004 (JT and Gauss-Legendre) versus baseline, with and without erfc
for variant in "" "-DGS_USE_ERFC"; do
  out=bench_v2; [ -n "$variant" ] && out=bench_v2_erfc
  $CC $CFLAGS $variant -I"$RINC" -Iv2 bench_v2.c v2/gsquad.c v2/hupdate.c v2/probpos.c v2/rprobrej.c \
    v2/gsbound1.c v2/gsbound.c $SRC/gridpts.c $SRC/h1.c base_renamed/*.c $LINK -o $out
done
echo "built: bench_baseline bench_iters bench_p3 bench_v2 bench_v2_erfc"
if [ "$1" = "--run" ]; then
  ./bench_baseline > ../results-harness-baseline.txt
  ./bench_iters > ../results-harness-iters.txt
  ./bench_p3 > ../results-harness-patch3.txt
  ./bench_v2 > ../results-harness-v2.txt
  ./bench_v2_erfc > ../results-harness-v2-erfc.txt
  echo "results written to analysis/results-harness-*.txt"
fi
