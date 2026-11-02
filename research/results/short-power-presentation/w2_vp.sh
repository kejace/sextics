# remote varpro batch: usage bash w2_vp.sh <d> <par> <seed0> <nbatches> <restarts_per_batch> <procs>
cd /tmp/hodge/results/short-power-presentation
d=$1; par=$2; s0=$3; nb=$4; rs=$5; pr=$6
for ((b=0; b<nb; b++)); do
  seed=$((s0 + b*100000))
  echo "== varpro d=$d par=$par seed=$seed $(date)"
  python3 w2_varpro.py $d $par --restarts $rs --procs $pr --seed $seed --maxiter 500 --deflate --tag remote > w2_varpro_d${d}_p${par}_s${seed}.out 2>&1
  tail -1 w2_varpro_d${d}_p${par}_s${seed}.out | cut -c1-300
done
echo "== done $(date)"
