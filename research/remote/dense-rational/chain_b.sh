#!/bin/bash
cd /tmp/hodge/results/dense-rational
python3 families.py deg6 > fam_deg6.log 2>&1
python3 families.py deg5 > fam_deg5.log 2>&1
python3 families.py deg4 > fam_deg4.log 2>&1
python3 families.py deg7 > fam_deg7.log 2>&1
python3 families.py deg8 > fam_deg8.log 2>&1
