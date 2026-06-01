#!/usr/bin/env bash
set -euo pipefail

old_commit="${1:-d11ef3d}"
new_commit="${2:-HEAD}"
wheel="${3:-11}"
rounds="${4:-1}"
repeat="${5:-1}"
iterations="${6:-3}"
workload="${7:-split-contract}"
graph_file="${8:-}"

work_root="$(mktemp -d /tmp/gc-standardize3-compare-XXXXXX)"
old_dir="${work_root}/old"
new_dir="${work_root}/new"
mkdir -p "${old_dir}" "${new_dir}"

cleanup() {
	rm -rf "${work_root}"
}
trap cleanup EXIT

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${graph_file}" && "${graph_file}" != /* ]]; then
	graph_file="${repo_root}/${graph_file}"
fi
if [[ -n "${graph_file}" ]]; then
	echo "graph-file workloads are not supported by the current standardize3_profile.cpp" >&2
	exit 1
fi

extract_metric() {
	local label="$1"
	local file="$2"
	sed -n "s|^${label} = ||p" "${file}" | sed 's/ s$//' | tail -n 1
}

run_snapshot() {
	local commit="$1"
	local dest="$2"
	local output_file="$3"

	git -C "${repo_root}" archive "${commit}" | tar -x -C "${dest}"
	cp "${repo_root}/tools/standardize3_profile.cpp" "${dest}/tools/standardize3_profile.cpp"
	make -B -C "${dest}" standardize3-profile >/dev/null
	(
		cd "${dest}"
		./standardize3_profile "${wheel}" "${rounds}" "${repeat}" "${iterations}" "${workload}"
	) | tee "${output_file}"
}

old_output="${work_root}/old.txt"
new_output="${work_root}/new.txt"

echo "Comparing standardize3 commits"
echo "old = ${old_commit}"
echo "new = ${new_commit}"
echo "wheel = W${wheel}"
echo "rounds = ${rounds}"
echo "repeat = ${repeat}"
echo "iterations = ${iterations}"
echo "workload = ${workload}"
if [[ -n "${graph_file}" ]]; then
	echo "graph file = ${graph_file}"
fi
echo

run_snapshot "${old_commit}" "${old_dir}" "${old_output}"
echo
run_snapshot "${new_commit}" "${new_dir}" "${new_output}"
echo

old_total="$(extract_metric "standardize3 average" "${old_output}")"
new_total="$(extract_metric "standardize3 average" "${new_output}")"
old_attempts="$(extract_metric "create final attempts average" "${old_output}")"
new_attempts="$(extract_metric "create final attempts average" "${new_output}")"
old_sort="$(extract_metric "sort/filter average" "${old_output}")"
new_sort="$(extract_metric "sort/filter average" "${new_output}")"

python3 - <<'PY' "${old_total}" "${new_total}" "${old_attempts}" "${new_attempts}" "${old_sort}" "${new_sort}"
import sys

old_total, new_total, old_attempts, new_attempts, old_sort, new_sort = map(float, sys.argv[1:])

def pct_change(old, new):
    return ((new / old) - 1.0) * 100.0

print("Summary")
print(f"standardize3 average: {old_total:.6f} s -> {new_total:.6f} s ({pct_change(old_total, new_total):+.2f}%)")
print(f"create final attempts: {old_attempts:.6f} s -> {new_attempts:.6f} s ({pct_change(old_attempts, new_attempts):+.2f}%)")
print(f"sort/filter: {old_sort:.6f} s -> {new_sort:.6f} s ({pct_change(old_sort, new_sort):+.2f}%)")
PY
