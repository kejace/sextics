#!/usr/bin/env bash
# push local capacity dir (sources, scripts) to the remote box
rsync -az --exclude out --exclude logs --exclude '*.olean' -e "ssh -i $HOME/.ssh/google_compute_engine -o StrictHostKeyChecking=accept-new" /tmp/hodge/lean/capacity/ kejace@35.229.16.59:/tmp/hodge/lean/capacity/
