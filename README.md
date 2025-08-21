# E-dna-sightings-trustscore

A tiny (WIP) tool that compares novel eDNA-based species sightings with external sightings databases to compute a **trust score** for each sighting’s reliability.

## Status
Early planning. Repo scaffold first; implementation to follow.

## MVP Goals
- Ingest eDNA results (CSV)
- Match taxa, location, and time against external sightings sources
- Compute a reproducible trust score
- Simple CLI + example notebook

## Quickstart
```bash
# clone and enter
git clone https://github.com/<your-org>/E-dna-sightings-trustscore.git
cd E-dna-sightings-trustscore

# (optional) create env
python -m venv .venv
# mac/linux
source .venv/bin/activate
# windows (PowerShell)
.venv\Scripts\Activate.ps1

# install deps (placeholder)
pip install -r requirements.txt
