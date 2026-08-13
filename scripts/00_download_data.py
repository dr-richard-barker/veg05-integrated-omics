#!/usr/bin/env python3
"""
Download processed RNA-seq data and microbiome FASTQ from NASA OSDR API.

Usage:
    python 00_download_data.py                    # Download everything
    python 00_download_data.py --study OSD-767    # RNA-seq only
    python 00_download_data.py --study OSD-766    # Microbiome only
"""

import argparse
import json
import os
import sys
import urllib.request
import time

OSDR_API_BASE = "https://osdr.nasa.gov/osdr/data/osd"
DOWNLOAD_BASE = "https://osdr.nasa.gov"

# Files to download from each study
RNA_SEQ_FILES = [
    "GLDS-709_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
    "OSD-767_metadata_OSD-767-ISA.zip",
    "GLDS-709_rna_seq_bulkRNASeq_v2_runsheet.csv",
    "GLDS-709_rna_seq_qc_metrics_GLbulkRNAseq.csv",
    "GLDS-709_rna_seq_software_versions_GLbulkRNAseq.md",
]

MICROBIOME_METADATA_FILES = [
    "OSD-766_metadata_OSD-766-ISA.zip",
]


def get_file_listings(osd_id):
    """Fetch complete file listing from OSDR API."""
    url = f"{OSDR_API_BASE}/files/{osd_id}/?all_files=true&page_size=1000"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = json.load(r)
    return data["studies"][f"OSD-{osd_id}"]["study_files"]


def download_file(remote_url, dest_path, max_retries=3):
    """Download a single file with retry logic."""
    full_url = DOWNLOAD_BASE + remote_url
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(full_url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=300) as r:
                content = r.read()
            with open(dest_path, "wb") as f:
                f.write(content)
            size_mb = os.path.getsize(dest_path) / 1e6
            print(f"  Downloaded: {os.path.basename(dest_path)} ({size_mb:.2f} MB)")
            return True
        except Exception as e:
            print(f"  Attempt {attempt+1} failed for {os.path.basename(dest_path)}: {e}")
            if attempt < max_retries - 1:
                time.sleep(5)
    return False


def download_rnaseq(data_dir):
    """Download processed RNA-seq files from OSD-767."""
    print("\n=== Downloading RNA-seq data (OSD-767) ===")
    dest_dir = os.path.join(data_dir, "rnaseq")
    os.makedirs(dest_dir, exist_ok=True)

    files = get_file_listings(767)
    downloaded = 0
    for target_name in RNA_SEQ_FILES:
        match = [f for f in files if f["file_name"] == target_name]
        if match:
            dest = os.path.join(dest_dir, target_name)
            if download_file(match[0]["remote_url"], dest):
                downloaded += 1
        else:
            print(f"  WARNING: {target_name} not found in OSD-767 file listing")

    print(f"RNA-seq: {downloaded}/{len(RNA_SEQ_FILES)} files downloaded")


def download_microbiome(data_dir):
    """Download microbiome FASTQ + metadata from OSD-766."""
    print("\n=== Downloading microbiome data (OSD-766) ===")
    fastq_dir = os.path.join(data_dir, "microbiome", "fastq")
    os.makedirs(fastq_dir, exist_ok=True)

    files = get_file_listings(766)

    # Download metadata first
    meta_dir = os.path.join(data_dir, "metadata")
    os.makedirs(meta_dir, exist_ok=True)
    for target_name in MICROBIOME_METADATA_FILES:
        match = [f for f in files if f["file_name"] == target_name]
        if match:
            dest = os.path.join(meta_dir, target_name)
            download_file(match[0]["remote_url"], dest)

    # Download all FASTQ files
    fastq_files = [f for f in files if f["file_name"].endswith(".fastq.gz")]
    print(f"\nDownloading {len(fastq_files)} FASTQ files...")

    downloaded = 0
    failed = []
    for i, f in enumerate(fastq_files):
        dest = os.path.join(fastq_dir, f["file_name"])
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            downloaded += 1
            continue
        if download_file(f["remote_url"], dest):
            downloaded += 1
        else:
            failed.append(f["file_name"])

        if (i + 1) % 50 == 0:
            print(f"  Progress: {i+1}/{len(fastq_files)} files processed")

    print(f"\nMicrobiome: {downloaded}/{len(fastq_files)} FASTQ files downloaded")
    if failed:
        print(f"  Failed: {len(failed)} files")
        for fn in failed[:10]:
            print(f"    - {fn}")


def main():
    parser = argparse.ArgumentParser(description="Download VEG-05 data from NASA OSDR")
    parser.add_argument("--study", choices=["OSD-767", "OSD-766", "all"], default="all")
    parser.add_argument("--data-dir", default="data")
    args = parser.parse_args()

    data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), args.data_dir)
    os.makedirs(data_dir, exist_ok=True)

    if args.study in ("OSD-767", "all"):
        download_rnaseq(data_dir)
    if args.study in ("OSD-766", "all"):
        download_microbiome(data_dir)

    print("\n=== Download complete ===")


if __name__ == "__main__":
    main()
