cd /tmp/hodge/results/short-power-presentation
while [ ! -f relax_d7_m_w2_chart.ms ]; do sleep 20; done
while pgrep -f "[w]2_relax.py 7 m" > /dev/null; do sleep 10; done
echo "mixed d=7 relaxation built $(date); starting subfamily restrictions"
python3 w2_subfam2.py relax_d7_m_w2 --maxr 5 --timeout 300 --tag w2 > w2_subfam_d7_m.out 2>&1
echo "mixed d=7 subfamilies done $(date)"
