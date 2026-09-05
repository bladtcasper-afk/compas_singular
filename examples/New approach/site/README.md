# Overview page

A visual summary of what the frame-field front end can mesh today: every domain in the
suite drawn from its real outline, with its route, patch count and element quality.

## Published at

<https://claude.ai/code/artifact/2384fac7-eeea-40bc-ba84-ef95b01b895b>

**To update it, republish to that same URL** — publishing `overview.html` without it mints
a new link and the old one goes stale.

## Regenerate

```
cd "compas_singular/examples/New approach/site"
C:/Users/Casper/anaconda3/envs/singular312/python.exe build_page.py
```

`build_page.py` calls `build_data.py` first, so one command does both stages.

| File | What it does |
|---|---|
| `build_data.py` | reads `../15_baseline.py` for the outlines and `../baseline.json` for the numbers, writes `domains.json` |
| `build_page.py` | reads `domains.json`, writes `overview.html` |
| `domains.json` | generated — SVG path data plus the per-domain figures |
| `overview.html` | generated — the page to publish |

## The rule this exists to enforce

**Nothing on the page is typed by hand.** Outlines are imported from `15_baseline.py`;
figures are read from `baseline.json`. Retyping an outline from memory rather than
importing it cost two rounds of false corrections: the suite samples the ellipse with 60
points, and the obvious guess — matching `DISC` — is 48. It measures a different polygon,
deterministically and plausibly.

So after changing anything that moves the baseline, run `15_baseline.py` first, then
rebuild here. The page follows the baseline; it never restates it.

## Design

Colours are the viewer's own legend, so the page and the 3D scenes read as one project:
separatrix orange as the single accent, cable purple for guide curves, patch-corner blue
for focus. Both light and dark themes are tokenised at `:root`.
