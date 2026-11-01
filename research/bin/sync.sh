#!/usr/bin/env bash
# Snapshot all scratch work into the research branch of the fork and push.
# Safe to run repeatedly (cron); never touches the upstream checkout.
set -u
WT=${HODGE_RESEARCH_WT:-$HOME/devel/sextics-research}
cd "$WT" || exit 1
R=research
mkdir -p $R/engine $R/notes $R/results $R/discussions $R/workflows $R/claude-memory $R/lean
rsync -a --delete --exclude __pycache__ --exclude venv "$HOME/devel/sextics-work/" $R/engine/
[ -d /tmp/hodge/notes ] && rsync -a /tmp/hodge/notes/ $R/notes/
[ -d /tmp/hodge/results ] && rsync -a --exclude '*.npz' --exclude '*.olean' --exclude '*.zip' --max-size=20m /tmp/hodge/results/ $R/results/
[ -d /tmp/hodge/discussions ] && rsync -a /tmp/hodge/discussions/ $R/discussions/
[ -d /tmp/hodge/lean ] && rsync -a --exclude '*.olean' --max-size=20m /tmp/hodge/lean/ $R/lean/
for p in "$HOME"/.claude/projects/*sextics*/workflows/scripts; do [ -d "$p" ] && rsync -a "$p/" $R/workflows/; done
for m in "$HOME"/.claude/projects/*sextics*/memory; do [ -d "$m" ] && rsync -a "$m/" $R/claude-memory/; done
git add -A $R >/dev/null 2>&1
if ! git diff --cached --quiet; then
  git -c user.name="${GIT_USER:-kejace}" commit -qm "sync $(hostname) $(date -u +%Y-%m-%dT%H:%MZ)" && git pull -q --rebase fork research 2>/dev/null; git push -q fork research 2>&1 | tail -1
fi
