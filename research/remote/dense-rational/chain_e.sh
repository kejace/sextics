#!/bin/bash
# resume the families process once the d=9 kernel run is done
cd /tmp/hodge/results/dense-rational
while pgrep -f "e1kernel.py direct 9" > /dev/null; do sleep 10; done
FP=$(pgrep -f "families.py deg6")
[ -n "$FP" ] && kill -CONT $FP
