#!/bin/bash

# Automated 3-party MPC experiment runner using tmux with 3-phase execution
# Phase 1: Online-only mode (-o) to determine resource requirements
# Phase 2: Preprocessing mode (-p) to precompute resources
# Phase 3: Online mode to run actual experiment with precomputed resources
#
# Usage: ./run_experiment_tmux.sh [--phase1-only|-O] <depth> <updates> <queries> [variant] [threads]
#
# Requirements: tmux must be installed
# Install: sudo apt-get install tmux  (Ubuntu/Debian)

set -e

# Check if tmux is installed
if ! command -v tmux &> /dev/null; then
    echo "Error: tmux is not installed!"
    echo "Please install it: sudo apt-get install tmux"
    exit 1
fi

# Option parsing
PHASE1_ONLY=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    -O|--phase1-only)
      PHASE1_ONLY=true
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [--phase1-only|-O] <depth> <updates> <queries> [variant] [threads]"
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "Unknown option: $1"
      exit 1
      ;;
    *)
      break
      ;;
  esac
done

# Default parameters
DEPTH=${1:-20}
UPDATES=${2:-10}
QUERIES=${3:-10}
VARIANT=${4:-segmenttree8}
THREADS=${5:-8}

# Create session name with timestamp
SESSION_NAME="segtree_exp_$(date +%Y%m%d_%H%M%S)"
RESOURCES_FILE="/tmp/resources_${SESSION_NAME}.txt"

echo "========================================="
if [ "$PHASE1_ONLY" = true ]; then
echo "Starting Phase 1 Only MPC Experiment"
else
echo "Starting 3-Phase MPC Experiment"
fi
echo "========================================="
echo "Variant: $VARIANT"
echo "Depth: $DEPTH"
echo "Updates: $UPDATES"
echo "Queries: $QUERIES"
echo "Threads: $THREADS"
echo "Session: $SESSION_NAME"
echo "========================================="
echo "Phase 1: Online-only mode (resource detection)"
if [ "$PHASE1_ONLY" != true ]; then
echo "Phase 2: Preprocessing mode (resource generation)"
echo "Phase 3: Online mode (actual experiment)"
fi
echo "========================================="
echo ""

# Ensure logs directory exists
mkdir -p logs logs/plots

# Create a new tmux session with first pane
tmux new-session -d -s "$SESSION_NAME" -n main

# Split window into three horizontal panes
tmux split-window -h -t "$SESSION_NAME:0.0"         # now 2 panes
tmux split-window -h -t "$SESSION_NAME:0.1"         # now 3 panes
tmux select-layout -t "$SESSION_NAME:0" even-horizontal

# Function to run only Phase 1 (online-only)
run_phase1_only() {
    echo "=== PHASE 1: Online-only mode (resource detection) ==="

    # Phase 1: Online-only mode (-o flag)
    # Pane 1 (middle): Player 0 - Logger (capture resources)
    tmux send-keys -t "$SESSION_NAME:0.0" "echo '=== PHASE 1: Player 0 (Logger) - Online-only mode ==='; sleep 1" C-m
    tmux send-keys -t "$SESSION_NAME:0.0" "./prac -o -t $THREADS 0 $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES 2>&1 | tee /tmp/phase1_output_${SESSION_NAME}.log" C-m

    # Pane 2 (right): Player 1
    tmux send-keys -t "$SESSION_NAME:0.1" "echo '=== PHASE 1: Player 1 - Online-only mode ==='; sleep 2" C-m
    tmux send-keys -t "$SESSION_NAME:0.1" "./prac -o -t $THREADS 1 localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

    # Pane 0 (left): Player 2 - Server
    tmux send-keys -t "$SESSION_NAME:0.2" "echo '=== PHASE 1: Player 2 (Server) - Online-only mode ==='; sleep 3" C-m
    tmux send-keys -t "$SESSION_NAME:0.2" "./prac -o -t $THREADS 2 localhost localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

    # Phase 1 only: extract resources and detach
    tmux send-keys -t "$SESSION_NAME:0.0" "echo 'Phase 1 only mode: Extracting resource requirements...'" C-m
    tmux send-keys -t "$SESSION_NAME:0.0" "grep 'Precomputed values used:' /tmp/phase1_output_${SESSION_NAME}.log | tail -1 | sed 's/.*T0 //' > $RESOURCES_FILE" C-m
    tmux send-keys -t "$SESSION_NAME:0.0" "echo 'Phase 1 only mode complete. Detaching in 2 seconds...'; sleep 2; tmux detach-client -s \"$SESSION_NAME\"" C-m
}



# Function to run all 3 phases
run_three_phases() {
    echo "=== PHASE 1: Online-only mode (resource detection) ==="

    # Phase 1: Online-only mode (-o flag)
    # Pane 1 (middle): Player 0 - Logger (capture resources)
    tmux send-keys -t "$SESSION_NAME:0.0" "echo '=== PHASE 1: Player 0 (Logger) - Online-only mode ==='; sleep 1" C-m
    tmux send-keys -t "$SESSION_NAME:0.0" "./prac -o -t $THREADS 0 $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES 2>&1 | tee /tmp/phase1_output_${SESSION_NAME}.log" C-m

    # Pane 2 (right): Player 1
    tmux send-keys -t "$SESSION_NAME:0.1" "echo '=== PHASE 1: Player 1 - Online-only mode ==='; sleep 2" C-m
    tmux send-keys -t "$SESSION_NAME:0.1" "./prac -o -t $THREADS 1 localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

    # Pane 0 (left): Player 2 - Server
    tmux send-keys -t "$SESSION_NAME:0.2" "echo '=== PHASE 1: Player 2 (Server) - Online-only mode ==='; sleep 3" C-m
    tmux send-keys -t "$SESSION_NAME:0.2" "./prac -o -t $THREADS 2 localhost localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

    # Wait for Phase 1 to complete and then run Phase 2
    tmux send-keys -t "$SESSION_NAME:0.0" "echo 'Phase 1 completed. Extracting resource requirements...'" C-m
    tmux send-keys -t "$SESSION_NAME:0.0" "grep 'Precomputed values used:' /tmp/phase1_output_${SESSION_NAME}.log | tail -1 | sed 's/.*T0 //' > $RESOURCES_FILE" C-m
    # tmux send-keys -t "$SESSION_NAME:0.0" "echo 'Waiting 3 seconds before Phase 2...'; sleep 3" C-m

    # Phase 2: Preprocessing mode (-p flag)
    tmux send-keys -t "$SESSION_NAME:0.0" "echo '=== PHASE 2: Preprocessing mode ==='" C-m
    tmux send-keys -t "$SESSION_NAME:0.0" "RESOURCES=\$(cat $RESOURCES_FILE); echo \"Captured resources: \$RESOURCES\"" C-m

    # Start preprocessing - Player 0 must start first to accept connections
    tmux send-keys -t "$SESSION_NAME:0.0" "echo '=== PHASE 2: Player 0 - Preprocessing ==='; sleep 1" C-m
    tmux send-keys -t "$SESSION_NAME:0.0" "./prac -t $THREADS -p 0" C-m

    tmux send-keys -t "$SESSION_NAME:0.1" "echo '=== PHASE 2: Player 1 - Preprocessing ==='; sleep 2" C-m
    tmux send-keys -t "$SESSION_NAME:0.1" "./prac -t $THREADS -p 1 localhost" C-m

    tmux send-keys -t "$SESSION_NAME:0.2" "echo '=== PHASE 2: Player 2 (Server) - Preprocessing ==='; sleep 3" C-m
    tmux send-keys -t "$SESSION_NAME:0.2" "bash -lc 'RESOURCES=\$(cat $RESOURCES_FILE); echo \"SHELL=\$SHELL\"; echo \"BASH_VERSION=\${BASH_VERSION:-none}\"; set -- \$RESOURCES; echo \"TokenCount=\$#\"; printf \"[arg:%s]\\n\" \"\$@\"; exec ./prac -t $THREADS -p 2 localhost localhost \"\$@\"'" C-m

    # Wait for Phase 2 to complete and then run Phase 3
    # tmux send-keys -t "$SESSION_NAME:0.0" "echo 'Phase 2 completed. Waiting 3 seconds before Phase 3...'; sleep 3" C-m

    # Phase 3: Online mode (no -o flag, using precomputed resources)
    tmux send-keys -t "$SESSION_NAME:0.0" "echo '=== PHASE 3: Online mode with precomputed resources ==='" C-m

    tmux send-keys -t "$SESSION_NAME:0.0" "echo '=== PHASE 3: Player 0 (Logger) - Online mode ==='; sleep 1" C-m
    tmux send-keys -t "$SESSION_NAME:0.0" "./prac -t $THREADS 0 $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

    tmux send-keys -t "$SESSION_NAME:0.1" "echo '=== PHASE 3: Player 1 - Online mode ==='; sleep 2" C-m
    tmux send-keys -t "$SESSION_NAME:0.1" "./prac -t $THREADS 1 localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

    tmux send-keys -t "$SESSION_NAME:0.2" "echo '=== PHASE 3: Player 2 (Server) - Online mode ==='; sleep 3" C-m
    tmux send-keys -t "$SESSION_NAME:0.2" "./prac -t $THREADS 2 localhost localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

    # Cleanup and detach after Phase 3
    tmux send-keys -t "$SESSION_NAME:0.0" "echo 'All 3 phases completed! Cleaning up...'; rm -f $RESOURCES_FILE /tmp/phase1_output_${SESSION_NAME}.log; rm -rf *t0; echo 'Detaching session in 2 seconds...'; sleep 2; tmux detach-client -s \"$SESSION_NAME\"" C-m
}

# Execute phases (optionally Phase 1 only)
if [ "$PHASE1_ONLY" = true ]; then
    echo "Phase 1 only mode enabled: Skipping Phases 2 and 3."
    run_phase1_only
else
    run_three_phases
fi

echo ""
if [ "$PHASE1_ONLY" = true ]; then
echo "Phase 1-only experiment started in tmux session: $SESSION_NAME"
echo ""
echo "Execution Flow:"
echo "  Phase 1: Online-only mode (-o) - Determines resource requirements"
else
echo "3-Phase experiment started in tmux session: $SESSION_NAME"
echo ""
echo "Execution Flow:"
echo "  Phase 1: Online-only mode (-o) - Determines resource requirements"
echo "  Phase 2: Preprocessing mode (-p) - Generates required resources"
echo "  Phase 3: Online mode - Runs experiment with precomputed resources"
fi
echo ""
echo "To attach to the session and view progress:"
echo "  tmux attach -t $SESSION_NAME"
echo ""
echo "Navigation tips:"
echo "  - Switch panes: Ctrl+b, then arrow keys"
echo "  - Detach: Ctrl+b, then d"
echo "  - Kill session: tmux kill-session -t $SESSION_NAME"
echo ""
echo "Log files will be created in logs/ directory (from Player 0)"
echo "Resources file: $RESOURCES_FILE"
echo ""
echo "Attaching to session in 2 seconds..."
sleep 2

# Attach to the session
tmux attach -t "$SESSION_NAME"

# After detaching or session ends, clean up
echo ""
if [ "$PHASE1_ONLY" = true ]; then
echo "Phase 1-only experiment completed. Cleaning up..."
else
echo "3-Phase experiment completed. Cleaning up..."
fi
tmux kill-session -t "$SESSION_NAME" 2>/dev/null || echo "Session already closed."
rm -f "$RESOURCES_FILE" "/tmp/phase1_output_${SESSION_NAME}.log" 2>/dev/null || true
echo "Session $SESSION_NAME has been deleted and temporary files cleaned up."
