#!/usr/bin/env bash
# kjob.sh <label> <cpuset> <gen args...>  -- generate Cap/<label>.lean with kron_gen.py (waits for the Kron lib olean), then time it
LABEL="$1"; CPUSET="$2"; shift 2
CAP=/tmp/hodge/lean/capacity
for i in $(seq 1 120); do [ -f $CAP/out/Reflect/Kron.olean ] && break; sleep 10; done
cd $CAP/gen && python3 kron_gen.py "$@" --name "$LABEL" && bash $CAP/run_module.sh "Cap.$LABEL" "$LABEL" "$CPUSET"
