#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 || $# > 2 )); then
	echo "usage: $0 LOOP [BENCHMARK_DIRECTORY]" >&2
	exit 2
fi

loop_number=$1
benchmark_directory=${2:-output/generation_benchmarks/loop_${loop_number}}
run_directory=${benchmark_directory}/geng_finalized
finalizer=./geng_finalizer_${loop_number}

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
if [[ ! -x ${finalizer} ]]; then
	echo "missing ${finalizer}; build it separately with:" >&2
	echo "  make geng-finalizer GRAPH_STAGE_GENERATOR_LOOP=${loop_number}" >&2
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

mkdir -p "${run_directory}/logs" "${run_directory}/graphs"
stages=${run_directory}/stages.tsv
dimensions=${run_directory}/generation_dimensions.tsv
echo -e "backend\tloop\tvertices\tedges\tgraph_count\tstage_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb" > "${stages}"
echo -e "loop\tvertices\tedges\tunsigned_classes\todd_edge_classes\teven_edge_odd_vertex_classes" > "${dimensions}"

max_vertices=$((2 * (loop_number - 1)))
for ((vertices=2; vertices<=max_vertices; ++vertices)); do
	edges=$((loop_number + vertices - 1))
	maximum_edges=$((vertices * (vertices - 1) / 2))
	if ((vertices < 4 || edges > maximum_edges || 2 * edges < 3 * vertices)); then
		echo -e "geng_finalized\t${loop_number}\t${vertices}\t${edges}\t0\t0\t0\t0\t0" >> "${stages}"
		echo -e "${loop_number}\t${vertices}\t${edges}\t0\t0\t0" >> "${dimensions}"
		continue
	fi
	log=${run_directory}/logs/V${vertices}.log
	metrics=${run_directory}/logs/V${vertices}.metrics
	counts=${run_directory}/logs/V${vertices}.counts
	output=${run_directory}/graphs/loop_${loop_number}_vertices_${vertices}_admissible.gcg
	/usr/bin/time -f "%e\t%U\t%S\t%M" -o "${metrics}" \
		bash -c 'set -o pipefail; "$1" -C -d3 "$3" "$4:$4" 2>"$6" | "$2" "$3" "$5"' \
		_ "${geng}" "${finalizer}" "${vertices}" "${edges}" "${output}" "${log}" \
		> "${counts}"
	read -r graph_count odd_count even_count < "${counts}"
	read -r wall_seconds user_seconds system_seconds max_rss_kb < "${metrics}"
	echo -e "geng_finalized\t${loop_number}\t${vertices}\t${edges}\t${graph_count}\t${wall_seconds}\t${user_seconds}\t${system_seconds}\t${max_rss_kb}" >> "${stages}"
	echo -e "${loop_number}\t${vertices}\t${edges}\t${graph_count}\t${odd_count}\t${even_count}" >> "${dimensions}"
done

{
	echo -e "backend\tloop\tworkers\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb"
	awk -F '\t' -v loop="${loop_number}" '
		NR > 1 {
			wall += $6; user_seconds += $7; system_seconds += $8
			if ($9 > rss) rss = $9
		}
		END { print "geng_finalized\t" loop "\t1\t" wall "\t" user_seconds "\t" system_seconds "\t" rss }
	' "${stages}"
} > "${run_directory}/total.tsv"

echo "finalized geng benchmark complete: ${run_directory}"
