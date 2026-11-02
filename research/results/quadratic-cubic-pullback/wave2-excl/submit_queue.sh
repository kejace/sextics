#!/bin/bash
# wave-2 excl submission queue (waits for remote capacity); log in /tmp/hodge/results/quadratic-cubic-pullback/wave2-excl/submit_queue.log
cd /home/kejace/devel/sextics-work
python3 -m hodge.rjob submit --wait 14400 --lane quadratic-cubic-pullback --name excl-qcp-kernelA --cpus 10 -- python3 wave2_excl_kernel.py kernelA 10 "1:0:1,1:1:1,2:0:1,2:1:1,3:0:1,3:1:1,4:0:1,4:1:1,5:0:1,5:1:1,6:0:1,6:1:1,7:0:1,7:1:1,8:0:1,8:1:1,9:0:1,9:1:1,10:0:1:b,10:1:1:b,11:0:1:b,11:1:1:b,12:0:1:b,12:1:1:b,1:0:2,1:1:2,2:0:2,2:1:2,3:0:2,3:1:2,4:0:2,4:1:2,5:0:2,5:1:2,6:0:2,6:1:2,7:0:2,7:1:2,8:0:2,8:1:2,9:0:2,9:1:2"
python3 -m hodge.rjob submit --wait 14400 --lane quadratic-cubic-pullback --name excl-qcp-relaxC --cpus 8 -- "python3 wave2_excl_relax.py relaxC5e 5 0; python3 wave2_excl_relax.py relaxC5o 5 1; python3 wave2_excl_relax.py relaxC5a 5 all; python3 wave2_excl_relax.py relaxC6e 6 0; python3 wave2_excl_relax.py relaxC6o 6 1; python3 wave2_excl_relax.py relaxC6a 6 all exact=0 msolve=0"
python3 -m hodge.rjob submit --wait 14400 --lane quadratic-cubic-pullback --name excl-qcp-polarB --cpus 16 -- python3 wave2_excl_polar.py polarB 16 1 10
echo QUEUE DONE
