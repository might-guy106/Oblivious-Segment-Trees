#!/bin/bash

# Performance experiment helper for Segment Tree implementations
# This script generates the commands needed to run a 3-party MPC experiment
# Usage: ./run_experiment.sh <depth> <updates> <queries> [variant]

# Default parameters
DEPTH=${1:-20}
UPDATES=${2:-10}
QUERIES=${3:-10}
VARIANT=${4:-segmenttree8}

echo "========================================="
echo "3-Party MPC Experiment Commands"
echo "========================================="
echo "Variant: $VARIANT"
echo "Depth: $DEPTH"
echo "Updates: $UPDATES"
echo "Queries: $QUERIES"
echo "========================================="
echo ""
echo "Run these commands in THREE SEPARATE TERMINALS:"
echo ""
echo "Terminal 1 (Player 0 - Logger):"
echo "  ./prac -o -t 8 0 $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES"
echo ""
echo "Terminal 2 (Player 1):"
echo "  ./prac -o -t 8 1 localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES"
echo ""
echo "Terminal 3 (Player 2 - Server):"
echo "  ./prac -o -t 8 2 localhost localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES"
echo ""
echo "========================================="
echo "NOTE:"
echo "- Start Terminal 3 first (server)"
echo "- Then start Terminal 2"
echo "- Finally start Terminal 1"
echo "- Log files will be created in logs/ directory (from Player 0 only)"
echo "========================================="
echo ""
echo "To analyze results after completion, run:"
echo "  python3 analyze_performance.py logs/performance_*.csv"
echo ""
