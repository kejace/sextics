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
[ -d /tmp/hodge/posts ] && mkdir -p $R/posts && rsync -a /tmp/hodge/posts/ $R/posts/
[ -d /tmp/hodge/remote ] && mkdir -p $R/remote && rsync -a --exclude '*.npz' --exclude '*.olean' --exclude '*.zip' --exclude '*.pkl' --max-size=20m /tmp/hodge/remote/ $R/remote/
[ -d /tmp/hodge/lean ] && rsync -a --exclude '*.olean' --max-size=20m /tmp/hodge/lean/ $R/lean/
for p in "$HOME"/.claude/projects/*sextics*/workflows/scripts; do [ -d "$p" ] && rsync -a "$p/" $R/workflows/; done
for m in "$HOME"/.claude/projects/*sextics*/memory; do [ -d "$m" ] && rsync -a "$m/" $R/claude-memory/; done
# per-host branch: brunos -> research (primary), anything else -> research-<shorthost>
case "$(hostname)" in brunos) BR=research;; *) BR=research-$(hostname -s | cut -c1-12);; esac
git rev-parse --verify -q "$BR" >/dev/null || git checkout -q -b "$BR"
[ "$(git branch --show-current)" = "$BR" ] || git checkout -q "$BR"
git add -A $R >/dev/null 2>&1
if ! git diff --cached --quiet; then
  git -c user.name="${GIT_USER:-kejace}" commit -qm "sync $(hostname -s) $(date -u +%Y-%m-%dT%H:%MZ)"
fi
git fetch -q fork "$BR" 2>/dev/null && git merge -q -X ours --no-edit "fork/$BR" 2>/dev/null
git push -q fork "$BR" 2>&1 | grep -v "^hint:" | tail -1
