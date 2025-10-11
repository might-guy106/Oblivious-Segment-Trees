#!/bin/bash

# Automated 3-party MPC experiment runner using tmux
# This script opens 3 tmux panes and runs the experiment automatically
# Usage: ./run_experiment_tmux.sh <depth> <updates> <queries> [variant]
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

# Default parameters
DEPTH=${1:-20}
UPDATES=${2:-10}
QUERIES=${3:-10}
VARIANT=${4:-segmenttree8}

# Create session name with timestamp
SESSION_NAME="segtree_exp_$(date +%Y%m%d_%H%M%S)"

echo "========================================="
echo "Starting 3-Party MPC Experiment"
echo "========================================="
echo "Variant: $VARIANT"
echo "Depth: $DEPTH"
echo "Updates: $UPDATES"
echo "Queries: $QUERIES"
echo "Session: $SESSION_NAME"
echo "========================================="
echo ""

# Ensure logs directory exists
mkdir -p logs logs/plots

# Create a new tmux session with first pane
tmux new-session -d -s "$SESSION_NAME" -n main

# Split window into three horizontal panes
tmux split-window -h -t "$SESSION:0.0"         # now 2 panes
tmux split-window -h -t "$SESSION:0.1"         # now 3 panes
tmux select-layout -t "$SESSION:0" even-horizontal

# Adjust pane sizes for better visibility
tmux select-layout -t "$SESSION_NAME" even-horizontal

# Send commands to each pane
# Pane 0 (left): Player 2 - Server (start first)
tmux send-keys -t "$SESSION_NAME:0.0" "echo 'Player 2 (Server) - Starting...'; sleep 1; ./prac -o -t 8 2 localhost localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

# Pane 1 (middle): Player 0 - Logger (start second)
tmux send-keys -t "$SESSION_NAME:0.1" "echo 'Player 0 (Logger) - Waiting for server...'; sleep 2; ./prac -o -t 8 0 $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES; echo 'Detaching the session in 2 seconds...'; sleep 2; tmux detach-client -s "$SESSION_NAME"" C-m

# Pane 2 (right): Player 1 (start last)
tmux send-keys -t "$SESSION_NAME:0.2" "echo 'Player 1 - Waiting for others...'; sleep 3; ./prac -o -t 8 1 localhost $VARIANT -d $DEPTH -u $UPDATES -q $QUERIES" C-m

echo ""
echo "Experiment started in tmux session: $SESSION_NAME"
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
echo ""
echo "Attaching to session in 2 seconds..."
sleep 2

# Attach to the session
tmux attach -t "$SESSION_NAME"

# After detaching or session ends, clean up
echo ""
echo "Experiment completed. Cleaning up tmux session..."
tmux kill-session -t "$SESSION_NAME" 2>/dev/null || echo "Session already closed."
echo "Session $SESSION_NAME has been deleted."
