"""Idempotent: add local-execution mode to hodge/rjob.py (runs on the compute box itself)."""
import sys
P = '/home/kejace/devel/sextics-work/hodge/rjob.py'
src = open(P).read()
if 'def is_local' in src:
    print('rjob: already patched'); sys.exit(0)
old_ssh = '''def ssh(cmd: str, check=False, capture=True, timeout=None):
    full = ['ssh', *SSH_OPTS, f'{USER}@{HOST}', cmd]
    return subprocess.run(full, check=check, capture_output=capture, text=True, timeout=timeout)


def rsync(src, dst, extra=()):
    cmd = ['rsync', '-az', '--exclude', '__pycache__', '--exclude', 'venv', '--exclude', '*.olean', '-e',
           'ssh ' + ' '.join(shlex.quote(o) for o in SSH_OPTS), *extra, src, dst]
    return subprocess.run(cmd, capture_output=True, text=True)
'''
new_ssh = '''def is_local():
    """True when this process already runs on the compute box (then everything executes locally, no ssh).

    Detection: HODGE_REMOTE_LOCAL=1 forces local mode, HODGE_REMOTE_LOCAL=0 forces ssh; otherwise the
    hostname is compared with HODGE_REMOTE_HOSTNAME (default 'hodge-sextics')."""
    v = os.environ.get('HODGE_REMOTE_LOCAL')
    if v is not None:
        return v not in ('', '0', 'false', 'no')
    import socket
    return socket.gethostname().split('.')[0] == os.environ.get('HODGE_REMOTE_HOSTNAME', 'hodge-sextics')


def ssh(cmd: str, check=False, capture=True, timeout=None):
    if is_local():
        return subprocess.run(['bash', '-c', cmd], check=check, capture_output=capture, text=True, timeout=timeout)
    full = ['ssh', *SSH_OPTS, f'{USER}@{HOST}', cmd]
    return subprocess.run(full, check=check, capture_output=capture, text=True, timeout=timeout)


def _strip_host(path):
    return path.split(':', 1)[1] if ':' in path and not path.startswith('/') else path


def rsync(src, dst, extra=()):
    if is_local():
        s, d = _strip_host(src), _strip_host(dst)
        if os.path.abspath(s.rstrip('/')) == os.path.abspath(d.rstrip('/')):
            return subprocess.CompletedProcess(['rsync-local-noop'], 0, '', '')
        cmd = ['rsync', '-a', '--exclude', '__pycache__', '--exclude', 'venv', '--exclude', '*.olean', *extra, s, d]
        return subprocess.run(cmd, capture_output=True, text=True)
    cmd = ['rsync', '-az', '--exclude', '__pycache__', '--exclude', 'venv', '--exclude', '*.olean', '-e',
           'ssh ' + ' '.join(shlex.quote(o) for o in SSH_OPTS), *extra, src, dst]
    return subprocess.run(cmd, capture_output=True, text=True)
'''
assert old_ssh in src, 'ssh/rsync block not found'
src = src.replace(old_ssh, new_ssh)
old_up = '''def up(lane=None):
    ssh('mkdir -p /tmp/hodge/results /tmp/hodge/jobs /tmp/hodge/notes ~/devel/sextics-work')
'''
new_up = old_up + '''    if is_local():
        print('up: running on the compute box itself; nothing to push')
        return
'''
assert old_up in src; src = src.replace(old_up, new_up)
old_pull = '''def pull(lane):
    dst = MIRROR / lane
'''
new_pull = '''def pull(lane):
    if is_local():
        print(f'pull {lane}: running on the compute box itself; lane dir {RESULTS / lane} is already local '
              f'(job logs in /tmp/hodge/jobs/)')
        return
    dst = MIRROR / lane
'''
assert old_pull in src; src = src.replace(old_pull, new_pull)
old_doc = "Conventions for agents: keep local processes light"
new_doc = ("When invoked ON the compute box (hostname hodge-sextics, or HODGE_REMOTE_LOCAL=1) every subcommand runs locally:\n"
           "no ssh, no rsync; `submit` still detaches the job with nohup/setsid and records its pid under /tmp/hodge/jobs.\n\n"
           + old_doc)
assert old_doc in src; src = src.replace(old_doc, new_doc, 1)
open(P, 'w').write(src); print('rjob: patched')
