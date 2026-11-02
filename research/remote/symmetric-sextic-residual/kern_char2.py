import sys, time, json, pickle
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from charsolve import *
F, A, q = load()
groups = {'S3S2_++': Group([[0,1,2],[4,5]], [1,1]), 'S3S2_-+': Group([[0,1,2],[4,5]], [-1,1]), 'S3S2_+-': Group([[0,1,2],[4,5]], [1,-1]), 'S3S2_--': Group([[0,1,2],[4,5]], [-1,-1]),
          'S2S3_++': Group([[1,2],[3,4,5]], [1,1]), 'S2S3_-+': Group([[1,2],[3,4,5]], [-1,1]), 'S2S3_+-': Group([[1,2],[3,4,5]], [1,-1]), 'S2S3_--': Group([[1,2],[3,4,5]], [-1,-1])}
for d in [int(a) for a in sys.argv[1:]] or [8, 9]:
    for lab, G_ in groups.items():
        for pD in (0, 1):
            t0 = time.time(); info, polys = e1_kernel_char(F, A, q, d, pD, G_)
            info['seconds'] = round(time.time() - t0, 1); info['type'] = lab
            print(json.dumps(info), flush=True)
            pickle.dump((info, polys), open(f'{RES}/kchar_{lab}_d{d}_p{pD}.pkl', 'wb'))
