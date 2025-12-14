#!/bin/bash

# Default Configuration (from README)
DEFAULT_DEST="segmentrees@172.27.21.172:/home/segmentrees/temp2"

# Use argument if provided, otherwise use default
if [ "$#" -ge 1 ]; then
    DEST="$1"
else
    DEST="$DEFAULT_DEST"
fi

if [ -z "$DEST" ]; then
    echo "Error: No destination provided and no default configured."
    echo "Usage: $0 [user@host:/path/to/destination]"
    exit 1
fi

echo "Syncing project to $DEST ..."
echo "Excluding: logs, .git, .vscode, build artifacts (*.o, *.d, prac binary)"

# rsync options:
# -a: archive mode (preserves permissions, times, symbolic links)
# -v: verbose
# -z: compress file data during the transfer
# -P: show progress during transfer
# --delete: delete extraneous files from dest dirs (optional, I'll allow the user to decide or omitted for safety)
# using ./ to sync contents of current dir

rsync -avzP --delete \
    --exclude 'logs' \
    --exclude '.venv' \
    --exclude '.git' \
    --exclude '.vscode' \
    --exclude '*.o' \
    --exclude '*.d' \
    --exclude 'prac' \
    ./ "$DEST"

echo "Sync complete."
