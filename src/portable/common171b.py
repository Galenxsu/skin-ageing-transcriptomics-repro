from pathlib import Path
import sys, os, csv, json, hashlib, datetime

# Clean-room paths are explicit environment inputs. No fallback to the historical
# project tree is permitted because that would invalidate isolation.
required = ['CLEANROOM_PHASE171A', 'CLEANROOM_PHASE171B', 'CLEANROOM_GCTX']
missing = [k for k in required if not os.environ.get(k)]
if missing:
    raise RuntimeError('Missing clean-room environment variables: ' + ','.join(missing))

A = Path(os.environ['CLEANROOM_PHASE171A']).resolve()
O = Path(os.environ['CLEANROOM_PHASE171B']).resolve()
M = Path(os.environ['CLEANROOM_GCTX']).resolve()
R = O.parent
H = O  # dependency overlay is stored inside the isolated Phase171B directory

for p in (A, O, M):
    if not p.exists():
        raise FileNotFoundError(p)

sys.dont_write_bytecode = True

def sha(p):
    with open(p, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def csvread(p):
    with open(p, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def csvwrite(p, rows, fields=None):
    with open(p, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader(); w.writerows(rows)

def js(p, x):
    p.write_text(json.dumps(x, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')

def jr(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def event(kind, **data):
    with open(O/'EXECUTION_EVENTS.jsonl', 'a', encoding='utf-8') as f:
        f.write(json.dumps({'utc':now(), 'event':kind, **data}, ensure_ascii=False, allow_nan=False)+'\n')

def dirs():
    for d in ('QA','RESULTS','CHECKPOINTS','NULL','MATRIX','SCRIPTS'):
        (O/d).mkdir(exist_ok=True)

