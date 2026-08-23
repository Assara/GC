#!/usr/bin/env bash
set -euo pipefail

if (( $# != 1 )); then
	echo "usage: $0 BENCHMARK_DIRECTORY" >&2
	exit 2
fi

benchmark_directory=$1
gc_stages=${benchmark_directory}/gc/stages.tsv
geng_stages=${benchmark_directory}/geng/stages.tsv

for required in \
	"${gc_stages}" \
	"${geng_stages}" \
	"${benchmark_directory}/gc/total.tsv" \
	"${benchmark_directory}/geng/total.tsv"; do
	if [[ ! -f ${required} ]]; then
		echo "missing benchmark result: ${required}" >&2
		exit 2
	fi
done

awk -F '\t' '
	NR == FNR {
		if (FNR > 1) gc[$3 FS $4]=$5
		next
	}
	FNR == 1 { next }
	{
		key=$3 FS $4
		if (!(key in gc) || gc[key] != $5) {
			printf "count mismatch V=%s E=%s: gc=%s geng=%s\n", $3, $4, (key in gc ? gc[key] : "missing"), $5 > "/dev/stderr"
			bad=1
		}
		seen[key]=1
	}
	END {
		for (key in gc) if (!(key in seen)) {
			printf "geng result missing for %s\n", key > "/dev/stderr"
			bad=1
		}
		exit bad
	}
' "${gc_stages}" "${geng_stages}"

gc_row=$(tail -n 1 "${benchmark_directory}/gc/total.tsv")
geng_row=$(tail -n 1 "${benchmark_directory}/geng/total.tsv")
gc_wall=$(awk -F '\t' '{print $4}' <<< "${gc_row}")
geng_wall=$(awk -F '\t' '{print $4}' <<< "${geng_row}")
ratio=$(awk -v gc="${gc_wall}" -v geng="${geng_wall}" 'BEGIN { if (geng == 0) print "inf"; else printf "%.6g", gc/geng }')

echo "Counts match at every (V,E) stage."
echo "GC total:   ${gc_row}"
echo "geng total: ${geng_row}"
echo "GC/geng wall-time ratio: ${ratio}"
