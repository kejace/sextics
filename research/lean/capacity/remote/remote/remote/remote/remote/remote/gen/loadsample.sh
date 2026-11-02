#!/usr/bin/env bash
while true; do echo "{\"ts\": \"$(date +%Y-%m-%dT%H:%M:%S)\", \"load1\": $(cut -d' ' -f1 /proc/loadavg)}" >> /tmp/hodge/lean/capacity/results/load.jsonl; sleep 30; done
