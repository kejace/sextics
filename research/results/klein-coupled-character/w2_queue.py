"""Local job queue for the lane: keeps at most MAXJ w2-* jobs running on the remote, submitting the next queued job
whenever the remote gate (load/memory) allows.  Jobs are read from w2_queue.jobs.json on every loop (so jobs can be
added/reordered without restarting).  Light: one ssh every ~90 s."""
import subprocess, time, json, sys, os
ENGINE = '/home/kejace/devel/sextics-work'
LANE = 'klein-coupled-character'
LANEDIR = '/tmp/hodge/results/klein-coupled-character'
MAXJ = 4
JOBS = LANEDIR + '/w2_queue.jobs.json'
STATE = LANEDIR + '/w2_queue.state.json'
LOG = LANEDIR + '/w2_queue.log'


def log(msg):
    with open(LOG, 'a') as f:
        f.write(time.strftime('%H:%M:%S ') + msg + '\n')


def remote_jobs():
    r = subprocess.run(['ssh', '-i', os.path.expanduser('~/.ssh/google_compute_engine'), '-o', 'ConnectTimeout=20',
                        'kejace@35.229.16.59',
                        'for p in /tmp/hodge/jobs/w2-*.pid; do n=$(basename $p .pid); pid=$(cat $p); if kill -0 $pid 2>/dev/null; then echo RUNNING $n; else echo DONE $n; fi; done'],
                       capture_output=True, text=True, timeout=60)
    running, done = set(), set()
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2:
            (running if parts[0] == 'RUNNING' else done).add(parts[1])
    return running, done


def main():
    submitted = set(json.load(open(STATE))) if os.path.exists(STATE) else set()
    idle_rounds = 0
    while True:
        try:
            jobs = json.load(open(JOBS))
        except Exception as ex:
            log(f'jobs file error {ex!r}'); time.sleep(60); continue
        try:
            running, done = remote_jobs()
        except Exception as ex:
            log(f'status error {ex!r}'); time.sleep(120); continue
        mine = {j['name'] for j in jobs}
        running = {n for n in running if n in mine}
        done = {n for n in done if n in mine}
        pending = [j for j in jobs if j['name'] not in submitted]
        if not pending and not running:
            idle_rounds += 1
            if idle_rounds > 20:
                log('all jobs submitted and finished; exiting'); break
        else:
            idle_rounds = 0
        if len(running) < MAXJ and pending:
            for j in pending:
                deps = j.get('deps', [])
                if any(dp not in done for dp in deps):
                    continue
                r = subprocess.run(['python3', '-m', 'hodge.rjob', 'submit', '--lane', LANE, '--name', j['name'], '--cpus', str(j.get('cpus', 8)), '--', j['cmd']],
                                   cwd=ENGINE, capture_output=True, text=True, timeout=600)
                if 'submitted ' in r.stdout:
                    submitted.add(j['name']); json.dump(sorted(submitted), open(STATE, 'w'))
                    log(f'submitted {j["name"]}; running={sorted(running)}')
                else:
                    log(f'gate closed for {j["name"]}: {(r.stderr or r.stdout).strip()[-160:]}')
                break
        time.sleep(90)


if __name__ == '__main__':
    main()
