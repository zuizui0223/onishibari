#!/usr/bin/env python3
"""Create contact pages for human evidence screening; do NOT auto-assign sex.

Uses published iNaturalist photo URLs from the source audit. Retain per-photo
attribution and license in TSV. Contact pages are for internal visual screening
and should not be interpreted as a random sample of individual plants.
"""
import argparse
import csv
import io
import json
from pathlib import Path
import textwrap
import time
from urllib.parse import urlsplit, urlunsplit
from urllib.request import Request, urlopen

from PIL import Image, ImageDraw, ImageFont, ImageOps

BG = "#fafafa"
HEADERS = {"User-Agent": "onishibari-evidence-photo-screen/0.1 (academic research)"}

def medium_url(url):
    # iNaturalist photo API returns square thumbnail; medium is a documented size.
    # Keep host and photo path untouched; only request a larger size.
    if not url:
        return None
    p = urlsplit(url)
    x = p.path.replace("/square.", "/medium.").replace("/thumb.", "/medium.")
    return urlunsplit((p.scheme, p.netloc, x, p.query, p.fragment))

def fetch_photo(url):
    with urlopen(Request(url, headers=HEADERS), timeout=25) as response:
        content = response.read(3000000)
    with Image.open(io.BytesIO(content)) as img:
        return img.convert("RGB")

def page(rows, i, dest, size=320, cols=4, page_n=12):
    block = rows[i:i + page_n]
    nr = (len(block)+cols-1)//cols
    tile_w, tile_h = size, size + 55
    canvas = Image.new("RGB", (cols*tile_w, nr*tile_h+50), BG)
    d = ImageDraw.Draw(canvas)
    d.text((10,10), f"Daphne pseudomezereum - image screening page {i//page_n+1} | provisional", fill="black")
    for j, row in enumerate(block):
        x = (j%cols)*tile_w
        y = (j//cols)*tile_h+48
        if row["img"] is not None:
            shot = ImageOps.contain(row["img"], (size-12, size-8))
            canvas.paste(shot, (x+(size-shot.width)//2, y+(size-shot.height)//2))
        else:
            d.rectangle((x+5, y+5, x+size-5, y+size-5), outline="grey", width=1)
            d.text((x+15,y+60), "unavailable or unlicensed", fill="black")
        d.text((x+10,y+size+4), f"ID:{row['id']}  {row['date']}  {row['license']}", fill="black")
        d.text((x+10,y+size+22), "sex: UNKNOWN (requires manual scoring)", fill="black")
    dest.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(dest, quality=86)

def main():
    arg = argparse.ArgumentParser()
    arg.add_argument("--source", default="data/derived/inat_onishibari_opportunity_snapshot.json")
    arg.add_argument("--output-dir", default="data/screening_gallery")
    opts = arg.parse_args()
    result = json.loads(Path(opts.source).read_text(encoding="utf-8"))
    obs = result.get("sample") or []
    out = Path(opts.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    license_ok = {"cc0","cc-by","cc-by-nc","cc-by-sa","cc-by-nc-sa"}
    for item in obs:
        pics = item.get("photo_refs") or []
        p = pics[0] if pics else {}
        license_code = p.get("license_code") or "unlicensed"
        raw_url = p.get("url")
        url = medium_url(raw_url)
        img = None
        error = None
        if license_code in license_ok and url:
            for attempt in range(2):
                try:
                    img = fetch_photo(url)
                    break
                except Exception as exc:
                    error = f"{type(exc).__name__}: {exc}"
                    time.sleep(attempt+1)
        else:
            error = "photo license missing/disallowed or URL absent"
        rows.append({
            "id": item.get("observation_id"),
            "date": item.get("observed_on") or "NA",
            "license": license_code,
            "observation_url": item.get("observation_url"),
            "url": raw_url,
            "attribution": p.get("attribution") or "",
            "error": error or "",
            "img": img,
        })
        time.sleep(0.15)
    for i in range(0, len(rows), 12):
        page(rows, i, out / f"page_{i//12+1:02d}.jpg")
    with (out / "screening_manifest.tsv").open("w", encoding="utf-8", newline="") as f:
        columns = ["id","date","license","observation_url","url","attribution","error",
                   "flowers_visible","anthers_visible","stigma_visible","sex_morph","confidence","annotator_notes"]
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for r in rows:
            clean = {k:r.get(k, "") for k in columns}
            for col in ["flowers_visible","anthers_visible","stigma_visible","sex_morph","confidence"]:
                clean[col] = "UNSCORED"
            writer.writerow(clean)
    print("PHOTO_SCREEN_AUDIT=" + json.dumps({
        "total":len(rows),
        "retrieved":sum(r["img"] is not None for r in rows),
        "not_retrieved_or_unlicensed":sum(r["img"] is None for r in rows),
        "pages": (len(rows)+11)//12,
        "warning":"No sex annotations have been made"
    }),flush=True)

if __name__ == "__main__":
    main()
