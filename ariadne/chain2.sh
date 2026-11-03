#!/bin/bash
cd /home/claude/ariadne
echo "[$(date -u +%H:%M)] Sherbrooke response surface (resumable) ..."
python3 response_surface.py --result results/sherbrooke_3shell_onsetP58.json --confirm-seeds 7 > /tmp/surf_sh.log 2>&1
echo "[$(date -u +%H:%M)] exit=$?"
if [ -f results/surface_sherbrooke_3shell_onsetP58.json ]; then
  python3 uncertainty.py --surface results/surface_sherbrooke_3shell_onsetP58.json --draws 2000 > /tmp/unc_sh.log 2>&1
  echo "[$(date -u +%H:%M)] sherbrooke uncertainty exit=$?"
  tail -4 /tmp/unc_sh.log
fi
echo "[$(date -u +%H:%M)] IDEAS sub-1 surface ..."
python3 response_surface.py --result results/ideas_sub-1_onset_rh-entorhinal.json --confirm-seeds 7 > /tmp/surf_id.log 2>&1
echo "[$(date -u +%H:%M)] exit=$? ; ALL DONE"
