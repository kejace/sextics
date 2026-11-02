cd /tmp/hodge/results/short-power-presentation
for spec in "7 0" "7 1" "8 0" "8 1"; do
  set -- $spec
  echo "== kernel d=$1 par=$2 $(date)"
  python3 w2_kernel.py $1 $2 --primes 3 --maxprimes 6 --threads 8 --numerators --tag kA > w2_kernel_d$1_p$2.out 2>&1
  tail -1 w2_kernel_d$1_p$2.out
done
echo "== done $(date)"
