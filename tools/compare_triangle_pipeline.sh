#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if (( $# != 2 )) || [[ ! $1 =~ ^([3-9]|10)$ ]]; then
    echo "Usage: bash $0 MAX_LOOP NEW_OUTPUT_DIRECTORY (3–10)" >&2
    exit 2
fi
maximum=$1
out=$2
geng=$(command -v geng || command -v nauty-geng)
labelg=$(command -v labelg || command -v nauty-labelg)
for binary in "build/triangle_seeds_L$maximum"; do
    [[ -x $binary ]] || { echo "Build first: $binary" >&2; exit 2; }
done
for ((loop=3; loop<=maximum; ++loop)); do
    [[ -x build/triangle_splits_L$loop ]] || { echo "Build split loop $loop first" >&2; exit 2; }
done
if [[ -e $out ]]; then
    echo "Use a fresh run directory; already exists: $out" >&2
    exit 2
fi
mkdir -p "$out/logs" "$out/geng" "$out/comparison"
printf 'backend\tloop\tvertices\twall\tuser\tsystem\n' > "$out/performance.tsv"
printf 'loop\tvertices\tedges\tpipeline\tgeng\tmissing\textra\tduplicates\tstatus\n' > "$out/counts.tsv"
TIMEFORMAT=$'%R\t%U\t%S'

# Generator timings include output writing and live logging, but exclude
# compilation, canonicalization and set comparison.
timed() {
    local backend=$1 loop=$2 vertices=$3 log=$4
    shift 4
    { time "$@" 2> "$log.stderr" | tee "$log"; } 2> "$out/time.tmp"
    printf '%s\t%s\t%s\t%s\n' "$backend" "$loop" "$vertices" "$(cat "$out/time.tmp")" >> "$out/performance.tsv"
}

echo '=== TRIANGLE SEEDS ==='
timed seeds all all "$out/logs/seeds.log" "build/triangle_seeds_L$maximum" "$out/seeds"
echo '=== VERTEX SPLITTING ==='
for ((loop=3; loop<=maximum; ++loop)); do
    timed splits "$loop" all "$out/logs/splits_L$loop.log" "build/triangle_splits_L$loop" "$out/seeds" "$out/splits_L$loop"
done

echo '=== GENG: 2-connected, minimum valence 3 ==='
failed=0
for ((loop=3; loop<=maximum; ++loop)); do
    for ((vertices=4; vertices<=2*(loop-1); ++vertices)); do
        edges=$((vertices+loop-1))
        (( edges <= vertices*(vertices-1)/2 )) || continue
        key=L${loop}_V${vertices}
        echo "geng L=$loop V=$vertices E=$edges"
        { time "$geng" -q -C -d3 "$vertices" "$edges:$edges" > "$out/geng/$key.g6" 2> "$out/logs/geng_$key.stderr"; } 2> "$out/time.tmp"
        printf 'geng\t%s\t%s\t%s\n' "$loop" "$vertices" "$(cat "$out/time.tmp")" >> "$out/performance.tsv"
        generated="$out/splits_L$loop/graphs_$key.g6"
        "$labelg" -q -g "$generated" | sort > "$out/comparison/$key.all.g6"
        uniq "$out/comparison/$key.all.g6" > "$out/comparison/$key.pipeline.g6"
        "$labelg" -q -g "$out/geng/$key.g6" | sort -u > "$out/comparison/$key.geng.g6"
        comm -23 "$out/comparison/$key.geng.g6" "$out/comparison/$key.pipeline.g6" > "$out/comparison/$key.missing.g6"
        comm -13 "$out/comparison/$key.geng.g6" "$out/comparison/$key.pipeline.g6" > "$out/comparison/$key.extra.g6"
        actual=$(wc -l < "$out/comparison/$key.pipeline.g6")
        reference=$(wc -l < "$out/comparison/$key.geng.g6")
        missing=$(wc -l < "$out/comparison/$key.missing.g6")
        extra=$(wc -l < "$out/comparison/$key.extra.g6")
        duplicates=$(( $(wc -l < "$generated") - actual ))
        status=match
        if (( missing || extra || duplicates )); then status=MISMATCH; failed=1; fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$loop" "$vertices" "$edges" "$actual" "$reference" "$missing" "$extra" "$duplicates" "$status" >> "$out/counts.tsv"
    done
done
{
    echo 'FINAL COUNTS — exact isomorphism-class comparison'
    awk -F '\t' '{for(i=1;i<=NF;i++) printf "%-12s",$i; print ""}
        NR>1 {a+=$4;b+=$5;m+=$6;e+=$7;d+=$8}
        END {printf "TOTAL: pipeline=%d geng=%d missing=%d extra=%d duplicates=%d\n",a,b,m,e,d}' "$out/counts.tsv"
    echo
    echo 'PERFORMANCE — seconds, compilation and comparison excluded'
    awk -F '\t' 'NR>1 {w[$1]+=$4;u[$1]+=$5;s[$1]+=$6}
        END {
            printf "%-20s %12s %12s %12s\n","Backend","Wall","User_CPU","System_CPU";
            printf "%-20s %12.3f %12.3f %12.3f\n","Triangle seeds",w["seeds"],u["seeds"],s["seeds"];
            printf "%-20s %12.3f %12.3f %12.3f\n","Splitting",w["splits"],u["splits"],s["splits"];
            printf "%-20s %12.3f %12.3f %12.3f\n","Pipeline total",w["seeds"]+w["splits"],u["seeds"]+u["splits"],s["seeds"]+s["splits"];
            printf "%-20s %12.3f %12.3f %12.3f\n","geng",w["geng"],u["geng"],s["geng"];
        }' "$out/performance.tsv"
} | tee "$out/summary.txt"
exit "$failed"
