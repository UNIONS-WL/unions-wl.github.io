#!/usr/bin/env python3
"""Append one fair-share sample for this cluster to production/fairshare/<cluster>.json.

    fairshare.py fir rrg-kyi def-sfabbro def-waerbeke def-mjhudson
    fairshare.py nibi def-sfabbro def-waerbeke def-mjhudson

Each sample is [UTC time, account LevelFS, our running jobs, our pending jobs].
LevelFS is Slurm's fair-share level of the account's CPU allocation: above 1 the
account has used less than its share and its jobs gain priority, below 1 it has
used more. Samples older than KEEP_DAYS are dropped. Standard library only.
"""
import json
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

KEEP_DAYS = 45
cluster, accounts = sys.argv[1], sys.argv[2:]
path = Path(__file__).resolve().parents[1] / "fairshare" / f"{cluster}.json"
doc = json.loads(path.read_text()) if path.exists() else {"cluster": cluster, "accounts": {}}


def lines(*command):
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout.splitlines()


now = datetime.now(timezone.utc)
stamp = now.strftime("%Y-%m-%dT%H:%MZ")
cutoff = (now - timedelta(days=KEEP_DAYS)).strftime("%Y-%m-%dT%H:%MZ")
level = {}
for row in lines("sshare", "-n", "-P", "-A", ",".join(a + "_cpu" for a in accounts), "-o", "Account,User,LevelFS"):
    account, user, levelfs = row.split("|")
    if not user:
        level[account.strip().removesuffix("_cpu")] = float(levelfs)
jobs = [row.split("|") for row in lines("squeue", "--me", "-h", "-o", "%a|%T")]
for account in accounts:
    count = lambda state: sum(a == account + "_cpu" and s == state for a, s in jobs)
    samples = doc["accounts"].setdefault(account, [])
    samples.append([stamp, level[account], count("RUNNING"), count("PENDING")])
    doc["accounts"][account] = [s for s in samples if s[0] >= cutoff]
doc["updated"] = stamp
path.parent.mkdir(exist_ok=True)
path.write_text(json.dumps(doc, separators=(",", ":")) + "\n")
