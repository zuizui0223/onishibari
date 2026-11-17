#!/usr/bin/env python3
"""Make attributed contact sheets of every *licensed* iNat image.

No sex or floral anatomy labels are inferred automatically. Photo files are
not committed to git; workflow artifacts expire. Sheet captions are linked to
full license and attribution in the TSV manifest.
"""
import argparse
import csv
import io
import json
from pathlib import Path
import time
from urllib.error import URLError, HTTPError
from urllib.parse import urlsplit, urlunsplit
from urllib.request import Request, urlopen

from PIL import Image, ImageDraw, ImageOps

LICENSES = {"cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa"}
HEADERS = {"User-Agent": "onishibari-research-photo-screen/0.2"}
TILE_W = 420
TILE_H = 350
COLUMNS = 3
ROWS = 3

def medium_url(url):
    parsed = urlsplit(url or "")
    if parsed.scheme != "https" or not parsed.netloc:
        return None
    if not (parsed.netloc.endswith("inaturalist.org") or
            parsed.netloc.endswith("amazonaws.com")):
        return None
    path = parsed.path.replace("/square.", "/medium.").replace("/thumb.", "/medium.")
    return urlunsplit((parsed.scheme, parsed.netloc, path, parsed.query, parsed.fragment))

def load_image(url):
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers=HEADERS), timeout=35) as resp:
                raw = resp.read(3500000)
            with Image.open(io.BytesIO(raw)) as p:
                return p.convert("RGB"), ""
        except (HTTPError, URLError, TimeoutError, OSError) as err:
            error = f"{type(err).__name__}: {str(err)[:120]}"
            time.sleep(0.6*(attempt+1))
    return None, error

def render(rows, dest, page):
    count = len(rows)
    height = ROWS*TILE_H+54
    sheet = Image.new("RGB", (COLUMNS*TILE_W, height), "white")
    draw = ImageDraw.Draw(sheet)
    draw.text((15, 14), f"Daphne pseudomezereum: multi-photo screen | Page {page} | SEX UNKNOWN", fill="#202020")
    for j, rec in enumerate(rows):
        x = (j % COLUMNS)*TILE_W
        y = (j // COLUMNS)*TILE_H + 51
        photo = rec["image"]
        if photo is None:
            draw.rectangle((x+8,y+8,x+TILE_W-8,y+TILE_H-60), outline="#aaaaaa")
            draw.text((x+20,y+30), "SKIPPED: license/download", fill="#555555")
        else:
            resized = ImageOps.contain(photo, (TILE_W-16, TILE_H-70))
            sheet.paste(resized, (x+(TILE_W-resized.width)//2, y+(TILE_H-70-resized.height)//2))
        draw.text((x+12,y+TILE_H-63),
                  f"obs {rec['observation_id']} photo {rec['photo_index']}/{rec['photo_total']}", fill="black")
        draw.text((x+12,y+TILE_H-44),
                  f"{rec['observed_on']} | {rec['quality_grade']} | {rec['license'] or 'unlicensed'}",
                  fill="black")
        draw.text((x+12,y+TILE_H-25), "anthers / stigma / sex = NOT YET SCORED", fill="#666666")
    dest.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(dest, format="JPEG", quality=88, optimize=True)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", default="data/derived/inat_all_photo_snapshot.json")
    p.add_argument("--output-dir", default="data/all_photo_gallery")
    p.add_argument("--max-images", type=int, default=200)
    args = p.parse_args()
    source = json.loads(Path(args.source).read_text(encoding="utf-8"))
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for obs in source.get("observations", []):
        for idx, photo in enumerate(obs.get("photos", []),1):
            records.append({
                "observation_id": obs["observation_id"],
                "observed_on": obs.get("observed_on") or "unknown",
                "quality_grade": obs.get("quality_grade") or "unknown",
                "photo_id": photo.get("photo_id"),
                "photo_index": idx,
                "photo_total": len(obs["photos"]),
                "license": (photo.get("license_code") or "").lower(),
                "attribution": photo.get("attribution") or "",
                "url": photo.get("url") or "",
                "observation_url": obs.get("url"),
                "image": None,
                "error": "",
            })
    if len(records) > args.max_images:
        raise ValueError(f"Too many public photos ({len(records)}) for max-images")
    columns = ["observation_id","observed_on","quality_grade","photo_id",
               "photo_index","photo_total","license","attribution","url",
               "observation_url","error","open_flower","anthers_visible",
               "stigma_visible","functional_sex","manual_confidence"]
    for i, r in enumerate(records):
        url = medium_url(r["url"])
        if r["license"] not in LICENSES or not url:
            r["error"] = "license absent or blocked, or invalid media URL"
        else:
            r["image"], r["error"] = load_image(url)
        if (i+1) % 15 == 0:
            print(f"IMAGE_PROGRESS={i+1}/{len(records)}",flush=True)
        time.sleep(.04)
    for start in range(0, len(records), COLUMNS*ROWS):
        render(records[start:start+COLUMNS*ROWS],
               out / f"all_photos_page_{start//(COLUMNS*ROWS)+1:02d}.jpg",
               start//(COLUMNS*ROWS)+1)
    with (out / "all_photo_screen_manifest.tsv").open("w",encoding="utf-8",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=columns,delimiter="\t")
        writer.writeheader()
        for r in records:
            clean={k:r.get(k,"") for k in columns}
            for k in ["open_flower","anthers_visible","stigma_visible","functional_sex",
                      "manual_confidence"]:
                clean[k]="UNSCORED"
            writer.writerow(clean)
    available=sum(r["image"] is not None for r in records)
    print("INAT_MULTI_PHOTO=" + json.dumps({
        "observations":len(source.get("observations",[])),
        "images_total":len(records),
        "image_thumbnails_retrieved":available,
        "license_skipped_or_failed":len(records)-available,
        "pages":(len(records)+COLUMNS*ROWS-1)//(COLUMNS*ROWS),
        "sex_diagnoses":0
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
