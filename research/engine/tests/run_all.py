"""Run all engine tests:  python3 tests/run_all.py   (or: pytest tests/)"""
import os, sys, time, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
mods = ['test_jacobian', 'test_symmetry', 'test_varpro', 'test_sweep']
fails = 0
for m in mods:
    t0 = time.time()
    print(f'== {m}', flush=True)
    mod = importlib.import_module(m)
    for name in dir(mod):
        if name.startswith('test_'):
            try:
                getattr(mod, name)()
            except Exception as ex:
                fails += 1
                import traceback; traceback.print_exc()
                print(f'  FAIL {name}: {ex}', flush=True)
    print(f'== {m} done ({time.time()-t0:.1f}s)', flush=True)
print('ALL OK' if not fails else f'{fails} FAILURES')
sys.exit(1 if fails else 0)
