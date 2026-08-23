#!/usr/bin/env bash

set -euo pipefail

script_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd -- "${script_directory}/.." && pwd)
cd -- "${repository_root}"

min_tadpoles=${MIN_TADPOLES:-3}
max_tadpoles=${1:-${MAX_TADPOLES:-5}}
dimension_directory=${2:-${DIMENSION_DIR:-output/graph_dimensions}}
dimension_table=${3:-${DIMENSION_TABLE:-${dimension_directory}/dimensions.tsv}}

case ${min_tadpoles} in
	''|*[!0-9]*)
		echo "error: MIN_TADPOLES must be an integer of at least 3" >&2
		exit 1
		;;
esac
case ${max_tadpoles} in
	''|*[!0-9]*)
		echo "error: MAX_TADPOLES must be an integer of at least 3" >&2
		exit 1
		;;
esac
if (( min_tadpoles < 3 )); then
	echo "error: MIN_TADPOLES must be an integer of at least 3" >&2
	exit 1
fi
if (( max_tadpoles < min_tadpoles )); then
	echo "error: MAX_TADPOLES must be at least MIN_TADPOLES" >&2
	exit 1
fi

mkdir -p -- "${dimension_directory}" "$(dirname -- "${dimension_table}")"
temporary_table=$(mktemp "${dimension_table}.tmp.XXXXXX")
cleanup() {
	rm -f -- "${temporary_table}"
}
trap cleanup EXIT

for (( loop_number = min_tadpoles; loop_number <= max_tadpoles; ++loop_number )); do
	echo "Building graph-stage generator for ${loop_number} tadpoles" >&2
	make graph-stage-generator GRAPH_STAGE_GENERATOR_LOOP="${loop_number}" >&2

	stage_directory="${dimension_directory}/loop_${loop_number}"
	mkdir -p -- "${stage_directory}"
	echo "Generating graph stages for ${loop_number} tadpoles in ${stage_directory}" >&2
	"./graph_stage_generator_${loop_number}" generate "${stage_directory}" >&2

	stage_dimension_log="${stage_directory}/generation_dimensions.tsv"
	if (( loop_number == min_tadpoles )); then
		cp -- "${stage_dimension_log}" "${temporary_table}"
	else
		sed '1d' "${stage_dimension_log}" >> "${temporary_table}"
	fi
done

mv -- "${temporary_table}" "${dimension_table}"
trap - EXIT

echo "Wrote graph dimension table to ${dimension_table}" >&2
cat -- "${dimension_table}"
