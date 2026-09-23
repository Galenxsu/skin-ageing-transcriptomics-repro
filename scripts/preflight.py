from pathlib import Path
import shutil,sys,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
free=shutil.disk_usage(ROOT).free/1024**3
status={'free_disk_gib':free,'required_free_disk_gib':60,'full_null_authorized_by_resources':free>=60}
print(json.dumps(status,indent=2))
if free<60:
    raise SystemExit('RESOURCE_GATE_BLOCKED: full null requires >=60 GiB free disk; do not reduce B.')
