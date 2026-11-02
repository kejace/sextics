#!/bin/bash
# queue.sh NAME CPUS WAITSECS -- command : wait for capacity, then force-submit if still busy
NAME=$1; CPUS=$2; WAITS=$3; shift 3; [ "$1" == "--" ] && shift
cd /home/kejace/devel/sextics-work
python3 -m hodge.rjob submit --lane rank-one-biquartic --name $NAME --cpus $CPUS --wait $WAITS -- "$@"
rc=$?
if [ $rc -ne 0 ]; then
  echo "capacity wait expired (rc=$rc); forcing submit of $NAME at $(date)"
  python3 -m hodge.rjob submit --lane rank-one-biquartic --name $NAME --cpus $CPUS --force -- "$@"
fi
