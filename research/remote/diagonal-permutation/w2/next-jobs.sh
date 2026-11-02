# Rotation plan (submit when a slot frees; 4 concurrent jobs max). All --force because the shared gate (load<=150) never opens.
cd /home/kejace/devel/sextics-work
S="python3 -m hodge.rjob submit --lane diagonal-permutation --cpus 8 --force --name search-diagonal-permutation"
# $S-j9x3  -- python3 w2pt_search.py --d 9 --tag j9x3 --extra kerfam --sector 1 --vanish X3 --workers 8 --restarts 150 --seed 51 --maxiter 300   (auto-submitted after k9exact)
# $S-j9y4  -- python3 w2pt_search.py --d 9 --tag j9y4 --extra kerfam --sector 1 --vanish Y4 --workers 8 --restarts 150 --seed 53 --maxiter 300
# $S-j9e   -- python3 w2pt_search.py --d 9 --tag j9e  --extra kerfam --workers 8 --restarts 60 --seed 61 --maxiter 300          (unrestricted, even)
# $S-j9o   -- python3 w2pt_search.py --d 9 --pD 1 --tag j9o --extra kerfam --workers 8 --restarts 60 --seed 63 --maxiter 300   (unrestricted, odd)
# $S-j9em  -- python3 w2pt_search.py --d 9 --tag j9em --sector -1 --workers 8 --restarts 150 --seed 71 --maxiter 300          (anti-invariant D: Jacobian only)
# kerfam10: python3 -m hodge.rjob submit --lane diagonal-permutation --cpus 8 --force --name search-diagonal-permutation-kerfam10 -- 'python3 w2_kerfam.py 10 0 && python3 w2_kerfam.py 10 1'
# $S-j10ep -- python3 w2pt_search.py --d 10 --tag j10ep --extra kerfam --sector 1 --workers 8 --restarts 40 --seed 81 --maxiter 300
