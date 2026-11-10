#!/bin/bash
# Sample this cluster's fair share and push it to the pages repo. Run hourly from cron:
#   7 * * * * /path/to/unions-wl.github.io/production/tools/fairshare-push.sh fir rrg-kyi def-sfabbro def-waerbeke def-mjhudson
# The checkout must be a clone used for nothing else, on main, able to push.
set -euo pipefail
export PATH=/opt/software/slurm/bin:/usr/bin:/bin:$PATH
cd "$(dirname "$(realpath "$0")")/../.."
cluster=$1
git pull -q --rebase
python3 production/tools/fairshare.py "$@"
git add "production/fairshare/$cluster.json"
git commit -q -m "production: fair share, $cluster [auto]"
for attempt in 1 2 3; do
    git push -q && exit 0
    git pull -q --rebase
done
exit 1
