#!/usr/bin/env bash
set -euo pipefail

if (( $# != 1 )); then
	echo "usage: $0 BENCHMARK_DIRECTORY" >&2
	exit 2
fi

benchmark_directory=$1
gc_dimensions=${benchmark_directory}/gc/graphs/generation_dimensions.tsv
geng_dimensions=${benchmark_directory}/geng_finalized/generation_dimensions.tsv
gc_total=${benchmark_directory}/gc/total.tsv
geng_total=${benchmark_directory}/geng_finalized/total.tsv

for file in "${gc_dimensions}" "${geng_dimensions}" "${gc_total}" "${geng_total}"; do
	if [[ ! -f ${file} ]]; then
		echo "missing benchmark result: ${file}" >&2
		exit 1
	fi
done

if ! diff -u "${gc_dimensions}" "${geng_dimensions}"; then
	echo "final dimensions differ" >&2
	exit 1
fi

gc_row=$(tail -n 1 "${gc_total}")
geng_row=$(tail -n 1 "${geng_total}")
gc_wall=$(awk -F '\t' '{print $4}' <<< "${gc_row}")
geng_wall=$(awk -F '\t' '{print $4}' <<< "${geng_row}")
ratio=$(awk -v gc="${gc_wall}" -v geng="${geng_wall}" \
	'BEGIN { if (geng == 0) print "inf"; else printf "%.6g", gc/geng }')

echo "Unsigned and both signed dimensions match at every (V,E) stage."
echo "GC total:             ${gc_row}"
echo "geng+finalizer total: ${geng_row}"
echo "GC/(geng+finalizer) wall-time ratio: ${ratio}"
