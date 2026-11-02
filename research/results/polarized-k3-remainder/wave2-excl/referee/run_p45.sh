#!/usr/bin/env bash
# sequential fixed-D replays (one process at a time), receipts -> referee/packet-replays.jsonl
ulimit -v 6000000
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl/
O=referee/packet-replays.jsonl
nice -n 10 python3 fixedD_koszul.py referee-fUU0-k3-d7 7 "[scale(mul(U[0],U[i]),Fr(1,42)) for i in range(3)] + [{},{},{}]" --target=k3 --out=$O
nice -n 10 python3 fixedD_koszul.py referee-fUU0sq-k3-d8 8 "[scale(mul(mul(U[0],U[0]),U[i]),Fr(1,42)) for i in range(3)] + [{},{},{}]" --target=k3 --out=$O
nice -n 10 python3 fixedD_koszul.py referee-fUU0sq-full-d8 8 "[scale(mul(mul(U[0],U[0]),U[i]),Fr(1,42)) for i in range(3)] + [{},{},{}]" --target=full --out=$O
echo P45-DONE
