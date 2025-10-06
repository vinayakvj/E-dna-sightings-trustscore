#!/usr/bin/env python3
"""
Download FishBase Parquet datasets with caching + resume.

Usage:
  python download_fishbase.py
  python download_fishbase.py --out data_cache
  python download_fishbase.py --overwrite --validate
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict
import requests
from tqdm import tqdm

# Optional validation dependencies:
try:
    import pyarrow.parquet as pq  # type: ignore
    HAVE_PYARROW = True
except Exception:
    HAVE_PYARROW = False


FISHBASE_URLS: Dict[str, str] = {
    "fishbase_species.parquet": (
        "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/"
        "data/fb/v24.07/parquet/species.parquet?download=true"
    ),
    "fishbase_families.parquet": (
        "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/"
        "data/fb/v24.07/parquet/families.parquet?download=true"
    ),
    "fishbase_synonyms.parquet": (
        "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/"
        "data/fb/v24.07/parquet/synonyms.parquet?download=true"
    ),
}


def _head_content_length(url: str, timeout: int = 20) -> int | None:
    try:
        r = requests.head(url, allow_redirects=True, timeout=timeout)
        if r.status_code == 200:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl is not None else None
    except requests.RequestException:
        pass
    return None


def _download_with_resume(url: str, dest: Path, overwrite: bool = False) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and not overwrite:
        # If the file already exists and overwrite is False, skip.
        return

    temp = dest.with_suffix(dest.suffix + ".part")
    resume_pos = temp.stat().st_size if temp.exists() else 0

    # Try to get total size for progress bar
    total_size = _head_content_length(url)
    headers = {}
    if resume_pos and total_size and resume_pos < total_size:
        headers["Range"] = f"bytes={resume_pos}-"

    with requests.get(url, stream=True, headers=headers, timeout=30) as r:
        # If server does not support Range, it will return 200 and we reset
        if r.status_code in (200, 206):
            mode = "ab" if r.status_code == 206 and resume_pos else "wb"
            if r.status_code == 200 and resume_pos:
                # Server ignored Range; start fresh
                mode = "wb"
                resume_pos = 0

            chunk_size = 1024 * 1024  # 1 MiB
            if total_size is None:
                # Try to infer from response header if HEAD failed
                cl = r.headers.get("Content-Length")
                total_size = int(cl) + resume_pos if cl else None

            with open(temp, mode) as f, tqdm(
                total=total_size if total_size is not None else None,
                initial=resume_pos,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=dest.name,
            ) as pbar:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
        else:
            raise RuntimeError(
                f"Download failed: HTTP {r.status_code} for {url}"
            )

    temp.replace(dest)


def validate_parquet(path: Path) -> None:
    if not HAVE_PYARROW:
        print(f"[validate] Skipped (pyarrow not installed): {path.name}")
        return
    try:
        # Quick metadata read + tiny scan to ensure file is readable
        meta = pq.ParquetFile(path)
        schema = meta.schema_arrow
        table = meta.read_row_groups([0]) if meta.num_row_groups > 0 else None
        nrows = meta.metadata.num_rows if meta.metadata else None
        print(f"[validate] OK: {path.name} | rows={nrows} | columns={[f.name for f in schema]}")
    except Exception as e:
        raise RuntimeError(f"[validate] Failed to read {path}: {e}") from e


def main() -> None:
    parser = argparse.ArgumentParser(description="Download FishBase Parquet datasets.")
    parser.add_argument(
        "--out",
        type=str,
        default="data_cache",
        help="Output directory for downloaded files (default: data_cache)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download files even if they already exist.",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate Parquet files after download (requires pyarrow).",
    )
    args = parser.parse_args()

    out_dir = Path(args.out).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading to: {out_dir}")

    for filename, url in FISHBASE_URLS.items():
        dest = out_dir / filename
        try:
            _download_with_resume(url, dest, overwrite=args.overwrite)
            size = dest.stat().st_size if dest.exists() else 0
            print(f"✓ {filename} ({size:,} bytes)")
            if args.validate:
                validate_parquet(dest)
        except Exception as e:
            # Clean up partial file if something went wrong
            part = dest.with_suffix(dest.suffix + ".part")
            if part.exists():
                try:
                    part.unlink()
                except Exception:
                    pass
            print(f"✗ Failed: {filename}\n  URL: {url}\n  Error: {e}")

    print("Done.")


if __name__ == "__main__":
    # Improve reliability on flaky networks
    requests.adapters.DEFAULT_RETRIES = 2
    os.environ.setdefault("NO_PROXY", "*")
    main()

