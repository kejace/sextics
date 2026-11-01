#!/usr/bin/env bash
# Restore scratch work from the research branch onto this machine (run after cloning the fork).
set -eu
WT=${HODGE_RESEARCH_WT:-$HOME/devel/sextics-research}
R="$WT/research"
mkdir -p "$HOME/devel/sextics-work" /tmp/hodge/notes /tmp/hodge/results /tmp/hodge/discussions /tmp/hodge/lean
rsync -a "$R/engine/" "$HOME/devel/sextics-work/"
rsync -a "$R/notes/" /tmp/hodge/notes/
rsync -a "$R/results/" /tmp/hodge/results/
rsync -a "$R/discussions/" /tmp/hodge/discussions/
[ -d "$R/lean" ] && rsync -a "$R/lean/" /tmp/hodge/lean/
MEM="$HOME/.claude/projects/-home-$(whoami)-devel-sextics/memory"
mkdir -p "$MEM" && rsync -a "$R/claude-memory/" "$MEM/"
echo "restored engine -> ~/devel/sextics-work, notes/results/discussions -> /tmp/hodge, memory -> $MEM"
