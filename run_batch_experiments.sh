#!/bin/bash

# Script to run batch experiments for various depths and configurations

set -e

DEPTHS=(4 8 12 16)
VARIANT="segmenttree9"
THREADS=8

echo "Starting batch experiments..."

for d in "${DEPTHS[@]}"; do
    echo "----------------------------------------------------------------"
    echo "Running Depth: $d, Updates: 5, Queries: 0"
    ./run_experiment_tmux.sh "$d" 5 0 "$VARIANT" "$THREADS"
    echo "Finished Depth: $d, Updates: 5, Queries: 0"
    sleep 2

    echo "----------------------------------------------------------------"
    echo "Running Depth: $d, Updates: 0, Queries: 5"
    ./run_experiment_tmux.sh "$d" 0 5 "$VARIANT" "$THREADS"
    echo "Finished Depth: $d, Updates: 0, Queries: 5"
    sleep 2
done

echo "All experiments completed."
