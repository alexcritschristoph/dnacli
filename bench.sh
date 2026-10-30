#!/usr/bin/env bash
# Run the user-action benchmarks and keep a timestamped log in benchmarks/.
# Compare two logs with benchstat, if it is installed.
#
#   ./bench.sh                 every benchmark
#   ./bench.sh 'Draw|Pan'      only benchmarks whose name matches
#   COUNT=5 BENCHTIME=2s ./bench.sh
#
# The benchmarks read the fixtures in src/testdata/, which are not in this
# repository, so this is a script for a checkout that has them.
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -d src/testdata ]; then
	echo "bench.sh needs the fixtures in src/testdata/, which are not in this repository." >&2
	exit 1
fi

pattern="${1:-.}"
count="${COUNT:-3}"
benchtime="${BENCHTIME:-500ms}"
cores=$(getconf _NPROCESSORS_ONLN)
if [ "$(uname -s)" = Darwin ]; then
	cpu=$(sysctl -n machdep.cpu.brand_string)
	sha256() { shasum -a 256; }
else
	cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs || true)
	sha256() { sha256sum; }
fi

mkdir -p benchmarks
log="benchmarks/$(date +%Y-%m-%d_%H-%M-%S).txt"

# Build the test binary and warm the page cache before the measured run.
go test ./src -run '^$' -bench '^BenchmarkStartup' -benchtime=1x -cpu 1 >/dev/null

{
	echo "# dnacli interaction benchmarks"
	echo "# date:    $(date -Iseconds)"
	echo "# host:    $(hostname)"
	echo "# cpu:     $cpu"
	echo "# cores:   $cores"
	echo "# go:      $(go version)"
	echo "# code:    $(cat src/*.go | sha256 | cut -c1-12)"
	echo "# filter:  $pattern"
	echo "# samples: $count x $benchtime"
	echo
	go test ./src -run '^$' -bench "$pattern" -benchmem -count="$count" -benchtime="$benchtime" -cpu 1 -timeout 20m
} | tee "$log"

previous=$(ls -1t benchmarks/*.txt 2>/dev/null | sed -n 2p || true)
if [ -n "$previous" ] && command -v benchstat >/dev/null; then
	echo
	benchstat "$previous" "$log"
fi
