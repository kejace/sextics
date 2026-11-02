#!/bin/bash
# push single files into the remote packet dir (no directory sync). usage: rpush.sh file1 [file2 ...]  (paths relative to wave2-excl/)
cd /tmp/hodge/results/cycle-q2/wave2-excl || exit 1
for f in "$@"; do
  rsync -az -e "ssh -i /home/kejace/.ssh/google_compute_engine -o StrictHostKeyChecking=accept-new" "$f" "kejace@35.229.16.59:/tmp/hodge/results/cycle-q2/wave2-excl/$f" && echo "pushed $f"
done
