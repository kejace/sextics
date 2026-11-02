#!/bin/bash
# Hit-rate calibration vs family dimension on the positive control (family-hesse-web), Koszul formulation.
cd /tmp/hodge/results/four-plus-two-residual
C="python3 wave2-search/fpt_control.py"
$C kctl-lm-g0 family-hesse-web 400 8 lm 200 0 koszul
$C kctl-lbfgs-g0 family-hesse-web 400 8 lbfgs 400 0 koszul
$C kctl-lm-g2 family-hesse-web 800 8 lm 200 2 koszul
$C kctl-lm-g1 family-hesse-web 1600 8 lm 200 1 koszul
$C kctl-lbfgs-g1 family-hesse-web 1600 8 lbfgs 400 1 koszul
echo "[koszul-ctl] done"
