cd /tmp/hodge/results/short-power-presentation/verify
python3 e1kernel_v.py 1000000007 1e 1o 2e 2o 3e 3o 4e 4o 5e 5o 6e 6o 7e 7o 8e 8o > e1kernel_v_1e9p7.out 2>&1
python3 e1kernel_v.py 1000000007 9e 9o >> e1kernel_v_1e9p7.out 2>&1
