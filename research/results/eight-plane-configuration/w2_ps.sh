#!/bin/bash
# list my lane's local processes (excluding this script)
ps -eo pid,ppid,etime,pcpu,cmd | grep -E "w2_(jvarpro|localchain|hess|fam|k10|kernel|sing)" | grep -v "w2_ps.sh" | grep -v grep | cut -c1-160
