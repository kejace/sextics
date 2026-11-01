#!/bin/bash
ulimit -v 14000000
cd /tmp/hodge/results/quadratic-cubic-pullback
exec python3 delta_served.py join012 join345 Fpart6 fermatX
