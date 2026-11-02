#!/usr/bin/env bash
# Batch 1: template ceiling. Each line: label | gen args | cpuset
cd /tmp/hodge/lean/capacity
CPU=4
while IFS='|' read -r label args; do
  label=$(echo $label); args=$(echo $args)
  [ -z "$label" ] && continue
  cs="$CPU,$((CPU+1))"; CPU=$((CPU+2))
  ./rsubmit.sh "cap-$label" "cd /tmp/hodge/lean/capacity/gen && python3 gen_capacity.py $args --name $label && bash /tmp/hodge/lean/capacity/run_module.sh Cap.$label $label $cs"
done <<'LIST'
e2m3k10  | --F m3 --d 6 --kD 10 --kH 1 --parts e2
e2m3k20  | --F m3 --d 6 --kD 20 --kH 1 --parts e2
e2m3k40  | --F m3 --d 6 --kD 40 --kH 1 --parts e2
e2m3k60  | --F m3 --d 8 --kD 60 --kH 1 --parts e2
e2m3k80  | --F m3 --d 8 --kD 80 --kH 1 --parts e2
e2m3k100 | --F m3 --d 8 --kD 100 --kH 1 --parts e2
e2m3k150 | --F m3 --d 10 --kD 150 --kH 1 --parts e2
e1m3k10  | --F m3 --d 6 --kD 10 --kH 1 --parts e1
e1m3k40  | --F m3 --d 6 --kD 40 --kH 2 --parts e1
e1m3k150 | --F m3 --d 8 --kD 150 --kH 2 --parts e1
e1dnk5   | --F dense --d 4 --kD 5 --kH 1 --parts e1
e1dnk10  | --F dense --d 4 --kD 10 --kH 1 --parts e1
e1dnk20  | --F dense --d 4 --kD 20 --kH 1 --parts e1
e1dnk40  | --F dense --d 6 --kD 40 --kH 1 --parts e1
e1dnk80  | --F dense --d 6 --kD 80 --kH 1 --parts e1
fullm3k20 | --F m3 --d 6 --kD 20 --kH 1 --parts full
fulldnk20 | --F dense --d 6 --kD 20 --kH 1 --parts full
LIST
