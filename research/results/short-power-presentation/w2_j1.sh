cd /tmp/hodge/results/short-power-presentation
python3 -c "import flint; flint.ctx.threads=4"
echo "== e1vec 7 odd"; python3 e1vec.py 7 1 > w2_e1vec_7o.out 2>&1
echo "== testD7 even (full solve + served space)"; python3 test_D.py e1vec_d7_p0.json --image --family "K7-even D7 (cone over line X=M=0)" > w2_testD7e.out 2>&1
echo "== testD7 odd"; python3 test_D.py e1vec_d7_p1.json --image --family "K7-odd" > w2_testD7o.out 2>&1
echo "== e1vec 8 even"; python3 e1vec.py 8 0 > w2_e1vec_8e.out 2>&1
echo "== done"
