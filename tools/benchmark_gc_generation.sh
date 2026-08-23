#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
	echo "usage: $0 LOOP THREADS [BENCHMARK_DIRECTORY]" >&2
	exit 2
fi

loop_number=$1
threads=$2
benchmark_directory=${3:-output/generation_benchmarks/loop_${loop_number}}
generator=./graph_stage_generator_${loop_number}
run_directory=${benchmark_directory}/gc
graph_directory=${run_directory}/graphs

if [[ ! -x ${generator} ]]; then
	echo "missing ${generator}; build it separately with:" >&2
	echo "  make graph-stage-generator GRAPH_STAGE_GENERATOR_LOOP=${loop_number}" >&2
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

mkdir -p "${graph_directory}"

/usr/bin/time \
	-f "%e\t%U\t%S\t%M" \
	-o "${run_directory}/process_metrics.tsv" \
	env OMP_NUM_THREADS="${threads}" \
	"${generator}" generate "${graph_directory}" \
	> "${run_directory}/generator.log"

{
	echo -e "backend\tloop\tworkers\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb"
	read -r wall_seconds user_seconds system_seconds max_rss_kb \
		< "${run_directory}/process_metrics.tsv"
	echo -e "gc\t${loop_number}\t${threads}\t${wall_seconds}\t${user_seconds}\t${system_seconds}\t${max_rss_kb}"
} > "${run_directory}/total.tsv"

{
	echo -e "backend\tloop\tvertices\tedges\tgraph_count\tstage_seconds"
	awk -v loop="${loop_number}" '
		BEGIN { OFS="\t" }
		$0 ~ /kind=admissible/ {
			vertices=edges=graphs=seconds=""
			for (i=1; i<=NF; ++i) {
				split($i, part, "=")
				if (part[1] == "vertices") vertices=part[2]
				if (part[1] == "edges") edges=part[2]
				if (part[1] == "graphs") graphs=part[2]
				if (part[1] == "stage_elapsed_seconds") seconds=part[2]
			}
			if (vertices != "") print "gc", loop, vertices, edges, graphs, seconds
		}
	' "${run_directory}/generator.log"
} > "${run_directory}/stages.tsv"

echo "GC benchmark complete: ${run_directory}"
