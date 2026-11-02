#!/bin/bash
# kill processes whose cmdline matches $1, except this script and its parent shell
for pid in $(pgrep -f "$1"); do
  if [ "$pid" != "$$" ] && [ "$pid" != "$PPID" ]; then kill "$pid" 2>/dev/null && echo "killed $pid"; fi
done
