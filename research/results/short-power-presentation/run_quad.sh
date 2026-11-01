cd /tmp/hodge/results/short-power-presentation
python3 quad_polar.py build 6 m 2>&1 | grep -v "hodge engine"
echo "== Singular d6 mixed"; timeout 1800 Singular -q polar_d6_pm.sing 2>&1 | head -80
python3 quad_polar.py build 7 0 2>&1 | grep -v "hodge engine"
python3 quad_polar.py build 7 1 2>&1 | grep -v "hodge engine"
echo "== Singular d7 parD=0"; timeout 3600 Singular -q polar_d7_p0.sing 2>&1 | head -120
echo "== Singular d7 parD=1"; timeout 3600 Singular -q polar_d7_p1.sing 2>&1 | head -120
