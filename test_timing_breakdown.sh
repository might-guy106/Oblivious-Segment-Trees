#!/bin/bash

# Quick test to see detailed timing breakdown
set -e

echo "Testing detailed timing breakdown for depth 8..."
./run_experiment_tmux.sh 8 0 1 segmenttree9 8

echo ""
echo "Showing the detailed log:"
LOG_FILE=$(ls -t logs/execution_online_st9_d8_u0_q1_*.log | head -1)
cat "$LOG_FILE"
