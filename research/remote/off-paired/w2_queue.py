"""Load-aware job queue runner for lane off-paired.

Queue file w2_queue.txt: one job per line  "<cpus> <name> <command...>" (command runs in the lane dir with
PYTHONPATH set, OMP threads = 1; the command itself uses multiprocessing for cpus > 1).
Policy: at most BASE_CPUS (3) cores sustained at any load; bursts up to BURST_CPUS (16) only while the 1-min load
average is below LOAD_MAX (140).  Running jobs are tracked via /tmp/hodge/jobs/<name>.pid.  Poll every 60 s.
"""
import os, sys, time, shlex, subprocess, json
LANE_DIR = '/tmp/hodge/results/off-paired'; JOBS = '/tmp/hodge/jobs'
QUEUE = f'{LANE_DIR}/w2_queue.txt'; DONE = f'{LANE_DIR}/w2_queue.done'
BASE_CPUS, BURST_CPUS, LOAD_MAX = 3, 16, 140.0
running = {}   # name -> (pid, cpus, t0)


def load1():
    return float(open('/proc/loadavg').read().split()[0])


def alive(pid):
    try:
        os.kill(pid, 0); return True
    except OSError:
        return False


def log(msg):
    with open(f'{JOBS}/search-off-paired-queue.log', 'a') as f:
        f.write(time.strftime('%H:%M:%S ') + msg + '\n')


def launch(cpus, name, cmd):
    env = dict(os.environ, PYTHONPATH='/home/kejace/devel/sextics-work', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               FLINT_NUM_THREADS='1')
    logf = open(f'{JOBS}/{name}.log', 'a')
    p = subprocess.Popen(['bash', '-c', cmd], cwd=LANE_DIR, env=env, stdout=logf, stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, start_new_session=True)
    open(f'{JOBS}/{name}.pid', 'w').write(str(p.pid))
    running[name] = (p.pid, cpus, time.time())
    log(f'LAUNCH {name} cpus={cpus} pid={p.pid} load={load1():.0f} :: {cmd[:160]}')


while True:
    for n in list(running):
        pid, cpus, t0 = running[n]
        if not alive(pid):
            log(f'FINISHED {n} after {time.time()-t0:.0f}s'); del running[n]
    used = sum(c for _, c, _ in running.values())
    ld = load1()
    cap = BURST_CPUS if ld < LOAD_MAX else BASE_CPUS
    if os.path.exists(QUEUE):
        lines = [l for l in open(QUEUE).read().splitlines() if l.strip() and not l.startswith('#')]
        launched = []
        for l in lines:
            parts = l.split(None, 2)
            if len(parts) < 3:
                launched.append(l); continue
            cpus, name, cmd = int(parts[0]), parts[1], parts[2]
            if name in running or os.path.exists(f'{JOBS}/{name}.pid'):
                if name not in running and os.path.exists(f'{JOBS}/{name}.pid'):
                    try:
                        pid = int(open(f'{JOBS}/{name}.pid').read())
                        if alive(pid):
                            running[name] = (pid, cpus, time.time()); used += cpus
                            launched.append(l); continue
                    except Exception:
                        pass
                else:
                    launched.append(l); continue
            if used + cpus <= cap:
                launch(cpus, name, cmd); used += cpus; launched.append(l)
        if launched:
            rest = [l for l in lines if l not in launched]
            with open(QUEUE, 'w') as f:
                f.write('\n'.join(rest) + ('\n' if rest else ''))
            with open(DONE, 'a') as f:
                f.write('\n'.join(launched) + '\n')
    if os.path.exists(f'{LANE_DIR}/w2_queue.stop'):
        log('STOP file seen; exiting'); break
    time.sleep(60)
