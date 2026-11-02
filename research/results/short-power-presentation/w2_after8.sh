cd /tmp/hodge/results/short-power-presentation
while pgrep -f "w2_kernel.py 8 0" > /dev/null; do sleep 30; done
python3 w2_kernel.py 8 0 --primes 8 --maxprimes 10 --threads 2 --numerators --tag more > w2_kernel_d8_p0_more.out 2>&1
echo "after8 done $(date)"
