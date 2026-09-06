#!/usr/bin/env bash
set -euo pipefail

# Also usable on a saved run; reporting never launches a generator.
write_reports() {
    local directory=$1 backend
    local reports=$directory/reports
    mkdir -p "$reports/counts" "$reports/performance"
    for backend in transient2 geng hairless; do
        if ! head -n 1 "$directory/counts.tsv" | tr "\t" "\n" | grep -qx "$backend"; then
            continue
        fi
        awk -F '\t' -v backend="$backend" 'BEGIN {
            print "loop\tvertices\tedges\tgraphs\tgeng_reference\tstatus"
        } NR == 1 {for (i=1;i<=NF;i++) column[$i]=i; next}
        NR > 1 {
            count=$(column[backend])
            status_column = (backend == "hairless" && column["hairless_status"]) ? column["hairless_status"] : column["status"]
            status = backend == "geng" ? "reference" : $(status_column)
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,count,$(column["geng"]),status
        }' "$directory/counts.tsv" > "$reports/counts/$backend.tsv"
        awk -F '\t' -v backend="$backend" 'BEGIN {
            print "loop\tvertices\twall_seconds\tuser_seconds\tsystem_seconds\tstatus"
        } NR > 1 && ($1 == backend || $1 == backend "_total" ||
                     $1 == backend "_failed" || $1 == backend "_total_incomplete") {
            status = ($1 ~ /failed|incomplete/) ? "INCOMPLETE" : "complete"
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", $2,$3,$4,$5,$6,status
        }' "$directory/performance.tsv" > "$reports/performance/$backend.tsv"
    done
    {
        echo "GRAPH COUNTS — 2-connected, minimum valence 3 (unsigned)"
        echo "Totals cover all compared loop orders; intermediate transient graphs are excluded."
        printf '%-14s %12s %12s %12s %12s\n' Backend Graphs Matches Mismatches Unavailable
        for backend in transient2 hairless geng; do
            [[ -f $reports/counts/$backend.tsv ]] || continue
            awk -F '\t' -v backend="$backend" 'NR > 1 {
                if ($4 ~ /^[0-9]+$/) graphs += $4
                matches += ($6 == "match")
                mismatches += ($6 == "MISMATCH")
                unavailable += ($6 != "match" && $6 != "MISMATCH" && $6 != "reference")
            } END {
                printf "%-14s %12d %12s %12s %12s\n", backend, graphs,
                    (backend == "geng" ? "reference" : matches+0),
                    (backend == "geng" ? "-" : mismatches+0),
                    (backend == "geng" ? "-" : unavailable+0)
            }' "$reports/counts/$backend.tsv"
        done
        echo "Matches/mismatches count (loop, vertices, edges) rows. Unavailable means missing or failed results; graph totals may be partial."
        printf '\nFINAL COUNTS BY LOOP AND VERTEX COUNT\n'
        awk -F '\t' '{for(i=1;i<=NF;i++) printf "%-17s", $i; printf "\n"}' "$directory/counts.tsv"
        echo "Detailed count tables: $reports/counts/"
    } > "$reports/counts.txt"
    {
        echo "PERFORMANCE — seconds, compilation excluded"
        printf '%-14s %12s %12s %12s %12s\n' Backend Wall User_CPU System_CPU Status
        for backend in transient2 hairless geng; do
            [[ -f $reports/performance/$backend.tsv ]] || continue
            awk -F '\t' -v backend="$backend" '$1 == "all" && $2 == "all" {
                printf "%-14s %12s %12s %12s %12s\n", backend,$3,$4,$5,$6
            }' "$reports/performance/$backend.tsv"
        done
        echo "Wall = elapsed time; user/system = CPU time."
        echo "Transient2 runs once through the loop limit. Hairless runs separately per loop and also computes automorphism/sign metadata."
        echo "Detailed timing tables: $reports/performance/"
    } > "$reports/performance.txt"
    cat "$reports/counts.txt"
    printf '\n'
    cat "$reports/performance.txt"
}

if [[ ${1:-} == --report ]]; then
    if (( $# != 2 )) || [[ ! -f $2/counts.tsv || ! -f $2/performance.tsv ]]; then
        echo "usage: bash $0 --report EXISTING_RUN_DIRECTORY (requires counts.tsv and performance.tsv)" >&2
        exit 2
    fi
    write_reports "$2"
    exit 0
fi

if (( $# < 3 || $# > 4 )) || [[ ! $1 =~ ^[3-9]$ ]]; then
    echo "usage: bash $0 MAX_LOOP NEW_RUN_DIRECTORY PIPELINE_BINARY_OR_DASH [HAIRLESS_BINARY_PREFIX] (loops 3-9; dash skips Transient2)" >&2
    exit 2
fi
max_loop=$1
run_directory=$2
pipeline=$3
hairless_prefix=${4:-}
hairless_threads=${HAIRLESS_THREADS:-1}
if [[ $pipeline == - && -z $hairless_prefix ]]; then
    echo "Hairless-only comparison requires HAIRLESS_BINARY_PREFIX" >&2
    exit 2
fi
if [[ -n $hairless_prefix ]]; then
    if [[ ! $hairless_threads =~ ^[1-9][0-9]*$ ]]; then
        echo "HAIRLESS_THREADS must be a positive integer" >&2
        exit 2
    fi
    for ((loop=3; loop<=max_loop; ++loop)); do
        if [[ ! -x ${hairless_prefix}${loop} ]]; then
            echo "Build the hairless generator first: ${hairless_prefix}${loop}" >&2
            exit 2
        fi
    done
fi
geng=${GENG:-}
if [[ -z $geng ]]; then
    geng=$(command -v geng || command -v nauty-geng || true)
fi
if [[ $pipeline != - && ! -x $pipeline ]] || [[ -z $geng ]] || ! command -v "$geng" >/dev/null; then
    echo "Build the pipeline first and install nauty (geng), or set GENG." >&2
    exit 2
fi
if [[ -e $run_directory ]]; then
    echo "Use a fresh run directory; already exists: $run_directory" >&2
    exit 2
fi
mkdir -p "$run_directory/geng"
counts=$run_directory/counts.tsv
performance=$run_directory/performance.tsv
printf 'loop\tvertices\tedges\ttransient2\tgeng\tstatus\n' > "$counts"
if [[ -n $hairless_prefix ]]; then
    printf 'loop\tvertices\tedges\ttransient2\tgeng\tstatus\thairless\thairless_status\n' > "$counts"
fi
if [[ $pipeline == - ]]; then
    printf 'loop\tvertices\tedges\thairless\tgeng\tstatus\n' > "$counts"
fi
printf 'backend\tloop\tvertices\twall_seconds\tuser_seconds\tsystem_seconds\n' > "$performance"
TIMEFORMAT=$'%R\t%U\t%S'
export LC_ALL=C

if [[ $pipeline != - ]]; then
echo "[Transient2] Starting through loop $max_loop. Full log: $run_directory/transient.log"
if { time "$pipeline" "$run_directory/transient" 2>&1 | tee "$run_directory/transient.log" | \
    awk '/^TransientGraph2/ {print; fflush()} /^Starting V=/ || /^Expanding V=/ {
        print "[Transient2] " $0; fflush()
    }'; } \
        2> "$run_directory/transient.time"; then
    printf 'transient2\tall\tall\t%s\n' "$(cat "$run_directory/transient.time")" >> "$performance"
else
    echo "TransientGraph2 failed. Partial output and timing remain in $run_directory; comparison stopped." >&2
    exit 1
fi

if ! grep -qx 'COUNT_SCOPE biconnected_min_degree_3' "$run_directory/transient.log"; then
    echo "Rebuild the pipeline: the comparison requires biconnected counts." >&2
    exit 2
fi
fi
failures=0
if [[ -n $hairless_prefix ]]; then
    mkdir -p "$run_directory/hairless"
    for ((loop=3; loop<=max_loop; ++loop)); do
        prefix=$run_directory/hairless/L${loop}
        echo "[Hairless] Starting L=$loop threads=$hairless_threads. Full log: $prefix.log"
        if { time OMP_NUM_THREADS="$hairless_threads" OMP_DYNAMIC=FALSE \
            "${hairless_prefix}${loop}" generate "$prefix" 2>&1 | tee "$prefix.log" | \
            awk '/^completed / && /kind=admissible/ {
                for(i=1;i<=NF;i++) if ($i ~ /^vertices=/) print "[Hairless] Finished " $i
                fflush()
            }'; } \
            2> "$prefix.time"; then
            printf 'hairless\t%s\tall\t%s\n' "$loop" "$(cat "$prefix.time")" >> "$performance"
        else
            printf 'hairless_failed\t%s\tall\t%s\n' "$loop" "$(cat "$prefix.time")" >> "$performance"
            touch "$prefix.failed"
            failures=$((failures + 1))
            echo "Hairless L=$loop failed; retaining its log and partial output and continuing the comparison."
        fi
    done
fi
for ((loop=3; loop<=max_loop; ++loop)); do
    for ((vertices=4; vertices<=2*(loop-1); ++vertices)); do
        edges=$((vertices + loop - 1))
        (( edges <= vertices*(vertices-1)/2 )) || continue
        actual=NA
        if [[ $pipeline != - ]]; then
        actual=$(awk -v l="$loop" -v v="$vertices" -v e="$edges" \
            '$1 == "COUNT" && $2 == l && $3 == v && $4 == e {print $5}' "$run_directory/transient.log")
        if [[ ! $actual =~ ^[0-9]+$ ]]; then
            echo "Missing or duplicate COUNT for L=$loop V=$vertices E=$edges; rebuild with the correct bounds." >&2
            exit 2
        fi
        fi
        prefix=$run_directory/geng/L${loop}_V${vertices}
        echo "[geng] L=$loop V=$vertices E=$edges"
        { time "$geng" -q -C -d3 "$vertices" "$edges:$edges" \
            > "$prefix.g6" 2> "$prefix.log"; } 2> "$prefix.time"
        reference=$(wc -l < "$prefix.g6")
        status=match
        if [[ $pipeline != - ]] && (( actual != reference )); then
            status=MISMATCH
            failures=$((failures + 1))
        fi
        if [[ -n $hairless_prefix ]]; then
            dimensions=$run_directory/hairless/L${loop}/generation_dimensions.tsv
            hairless_count=
            if [[ -f $dimensions ]]; then
                hairless_count=$(awk -F '\t' -v l="$loop" -v v="$vertices" -v e="$edges" \
                    '$1 == l && $2 == v && $3 == e {print $4}' "$dimensions")
            fi
            hairless_status=match
            if [[ ! $hairless_count =~ ^[0-9]+$ ]]; then
                hairless_count=NA
                hairless_status=UNAVAILABLE
            elif [[ -e $run_directory/hairless/L${loop}.failed ]]; then
                hairless_status=FAILED_RUN
            elif (( hairless_count != reference )); then
                hairless_status=MISMATCH
            fi
            if [[ $hairless_status != match ]]; then failures=$((failures + 1)); fi
            if [[ $pipeline == - ]]; then
                printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$loop" "$vertices" "$edges" \
                    "$hairless_count" "$reference" "$hairless_status" >> "$counts"
            else
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$loop" "$vertices" "$edges" "$actual" "$reference" "$status" \
                "$hairless_count" "$hairless_status" >> "$counts"
            fi
        else
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$loop" "$vertices" "$edges" "$actual" "$reference" "$status" >> "$counts"
        fi
        printf 'geng\t%s\t%s\t%s\n' "$loop" "$vertices" "$(cat "$prefix.time")" >> "$performance"
    done
done

geng_total=$(awk -F '\t' '$1 == "geng" {wall += $4; user += $5; system_time += $6}
    END {printf "geng_total\tall\tall\t%.3f\t%.3f\t%.3f\n", wall, user, system_time}' "$performance")
printf '%s\n' "$geng_total" >> "$performance"
if [[ -n $hairless_prefix ]]; then
    awk -F '\t' '$1 == "hairless" || $1 == "hairless_failed" {
        wall += $4; user += $5; system_time += $6; failed += ($1 == "hairless_failed")
    } END {
        printf "%s\tall\tall\t%.3f\t%.3f\t%.3f\n", \
            (failed ? "hairless_total_incomplete" : "hairless_total"), wall, user, system_time
    }' "$performance" >> "$performance"
fi
printf '\n'
write_reports "$run_directory"
echo "Count mismatches or run/report failures: $failures."
(( failures == 0 ))
