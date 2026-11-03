#!/bin/bash
# Full overnight chain. Each step is resumable, so a container restart costs minutes, not hours.
cd /home/claude/ariadne
SURF=results/surface_stanford_hardi_onsetP39.json

echo "[$(date -u +%H:%M)] response surface (resumable) ..."
python3 response_surface.py --result results/stanford_hardi_onsetP39.json > /tmp/surf_s.log 2>&1
echo "[$(date -u +%H:%M)] surface exit=$? ; $(tail -2 /tmp/surf_s.log | head -1)"

if [ -f "$SURF" ]; then
  echo "[$(date -u +%H:%M)] uncertainty analysis ..."
  python3 uncertainty.py --surface "$SURF" --draws 2000 > /tmp/unc.log 2>&1
  echo "[$(date -u +%H:%M)] uncertainty exit=$?"
  python3 make_report.py --uncertainty results/uncertainty_stanford_hardi_onsetP39.json > /tmp/report.log 2>&1
  echo "[$(date -u +%H:%M)] report exit=$?"
  tail -6 /tmp/unc.log
else
  echo "[$(date -u +%H:%M)] FAILED: no surface produced"
fi

echo "[$(date -u +%H:%M)] virtual surgery, pre-registered (resumes per patient) ..."
python3 virtual_surgery2.py --workers 2 > /tmp/vs2.log 2>&1
echo "[$(date -u +%H:%M)] virtual surgery exit=$?"
tail -12 /tmp/vs2.log

echo "[$(date -u +%H:%M)] ALL DONE"
