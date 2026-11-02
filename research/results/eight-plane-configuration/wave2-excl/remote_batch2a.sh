#!/bin/bash
# wave-2 excl batch 2a (remote): exact Y-partials + Koszul-reduced E2 sampling at d = 7, 8 (random, Hamiltonian, gradient, arrangement fields)
cd /tmp/hodge/results/eight-plane-configuration
TAG=excl-epc-batch2a
python3 wave2-excl/w2x_koszul.py exact $TAG partials > wave2-excl/$TAG-partials-exact.out 2>&1
python3 wave2-excl/w2x_jsample.py 7 8 $TAG > wave2-excl/$TAG-jsample7.out 2>&1
python3 wave2-excl/w2x_jsample.py 8 8 $TAG > wave2-excl/$TAG-jsample8.out 2>&1
