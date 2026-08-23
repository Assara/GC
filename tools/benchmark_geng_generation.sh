#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 || $# > 2 )); then
	echo "usage: $0 LOOP [BENCHMARK_DIRECTORY]" >&2
	exit 2
fi

loop_number=$1
benchmark_directory=${2:-output/generation_benchmarks/loop_${loop_number}}
run_directory=${benchmark_directory}/geng

if [[ -n ${GENG:-} ]]; then
	geng=${GENG}
elif command -v nauty-geng >/dev/null 2>&1; then
	geng=$(command -v nauty-geng)
elif command -v geng >/dev/null 2>&1; then
	geng=$(command -v geng)
else
	echo "geng was not found; set GENG=/path/to/geng" >&2
	exit 2
fi
if [[ -e ${run_directory} ]]; then
	echo "refusing to overwrite existing benchmark directory: ${run_directory}" >&2
	exit 2
fi
if [[ ! -x /usr/bin/time ]]; then
	echo "GNU /usr/bin/time is required" >&2
	exit 2
fi

mkdir -p "${run_directory}/logs"
stages=${run_directory}/stages.tsv
echo -e "backend\tloop\tvertices\tedges\tgraph_count\tstage_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb" > "${stages}"

max_vertices=$((2 * (loop_number - 1)))
for ((vertices=2; vertices<=max_vertices; ++vertices)); do
	edges=$((loop_number + vertices - 1))
	log=${run_directory}/logs/V${vertices}.log
	metrics=${run_directory}/logs/V${vertices}.metrics
	maximum_edges=$((vertices * (vertices - 1) / 2))
	if ((vertices < 4 || edges > maximum_edges || 2 * edges < 3 * vertices)); then
		echo -e "geng\t${loop_number}\t${vertices}\t${edges}\t0\t0\t0\t0\t0" >> "${stages}"
		continue
	fi
	/usr/bin/time \
		-f "%e\t%U\t%S\t%M" \
		-o "${metrics}" \
		"${geng}" -C -d3 -u "${vertices}" "${edges}:${edges}" \
		2> "${log}"
	graph_count=$(awk '$1 == ">Z" { print $2 }' "${log}")
	if [[ -z ${graph_count} ]]; then
		echo "could not parse geng count from ${log}" >&2
		exit 1
	fi
	read -r wall_seconds user_seconds system_seconds max_rss_kb < "${metrics}"
	echo -e "geng\t${loop_number}\t${vertices}\t${edges}\t${graph_count}\t${wall_seconds}\t${user_seconds}\t${system_seconds}\t${max_rss_kb}" >> "${stages}"
done

{
	echo -e "backend\tloop\tworkers\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb"
	awk -F '\t' -v loop="${loop_number}" '
		NR > 1 {
			wall += $6; user_seconds += $7; system_seconds += $8
			if ($9 > rss) rss = $9
		}
		END { print "geng\t" loop "\t1\t" wall "\t" user_seconds "\t" system_seconds "\t" rss }
	' "${stages}"
} > "${run_directory}/total.tsv"

echo "geng benchmark complete: ${run_directory}"
