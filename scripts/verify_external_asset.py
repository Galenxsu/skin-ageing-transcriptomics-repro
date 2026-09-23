from pathlib import Path
import hashlib,sys
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
if len(sys.argv)!=3: raise SystemExit('usage: verify_external_asset.py FILE EXPECTED_SHA256')
got=sha(sys.argv[1]);print(got)
if got.lower()!=sys.argv[2].lower():raise SystemExit('SHA256_MISMATCH')
