#!/bin/bash
# structured irreducible candidates (discriminant of the u-cubic, Fermat cones, transports): E1 screen + exact solve of survivors
set -x
python3 rcr2.py e1 --d 9 10 11 12 --sec all
python3 rcr2.py screen --file cands_w2.txt
echo SCREEN_DONE
