#!/usr/bin/env python3
"""
Parse ISA-Tab metadata from OSD-767 (RNA-seq) and OSD-766 (microbiome),
build standardized sample metadata CSVs, and create a crosswalk matching
RNA-seq samples to microbiome samples by Flight/Ground x Plant x Tissue x Light.
"""

import csv
import os
import re
import zipfile
import json
import urllib.request
from collections import Counter

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "data")
META_DIR = os.path.join(DATA_DIR, "metadata")
FASTQ_DIR = os.path.join(DATA_DIR, "microbiome", "fastq")


def parse_rnaseq_metadata():
    """Parse RNA-seq metadata from ISA zip + runsheet."""
    print("=== Parsing RNA-seq metadata (OSD-767) ===")

    # Use the runsheet which has clean factor columns
    runsheet_path = os.path.join(DATA_DIR, "rnaseq", "GLDS-709_rna_seq_bulkRNASeq_v2_runsheet.csv")
    samples = []

    with open(runsheet_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_name = row["Sample Name"]
            # Parse: VEG-05-Flt-SN01-Adv-Root-Red_L008
            base = sample_name.replace("_L008", "")
            m = re.match(r'VEG-05-(Flt|Gnd)-SN(\d+)-(Adv-Root|Leaf)-(Red|Blue)', base)
            if not m:
                print(f"  WARNING: could not parse {sample_name}")
                continue

            flight = "Flight" if m.group(1) == "Flt" else "Ground"
            plant = f"SN{m.group(2)}"
            tissue = m.group(3)
            light = m.group(4)

            samples.append({
                "sample_id": sample_name,
                "study": "OSD-767",
                "organism": "Solanum lycopersicum",
                "cultivar": "Red Robin",
                "flight": flight,
                "plant": plant,
                "tissue": tissue,
                "light": light,
                "treatment": f"{flight}_{light}",
                "paired_end": row.get("paired_end", "True"),
                "has_ercc": row.get("has_ERCC", "False"),
            })

    print(f"  Parsed {len(samples)} RNA-seq samples")
    print(f"  Flight/Ground: {dict(Counter(s['flight'] for s in samples))}")
    print(f"  Tissue: {dict(Counter(s['tissue'] for s in samples))}")
    print(f"  Light: {dict(Counter(s['light'] for s in samples))}")
    print(f"  Treatment: {dict(Counter(s['treatment'] for s in samples))}")

    out_path = os.path.join(META_DIR, "sample_metadata_rnaseq.csv")
    fields = ["sample_id", "study", "organism", "cultivar", "flight", "plant",
              "tissue", "light", "treatment", "paired_end", "has_ercc"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(samples)
    print(f"  Written: {out_path}")
    return samples


def parse_microbiome_metadata():
    """Parse microbiome metadata from FASTQ filenames."""
    print("\n=== Parsing microbiome metadata (OSD-766) ===")

    fastq_files = [f for f in os.listdir(FASTQ_DIR) if f.endswith(".fastq.gz") and "_R1_" in f]
    print(f"  Found {len(fastq_files)} R1 FASTQ files (samples)")

    samples = []
    unparseable = []

    for fn in sorted(fastq_files):
        # Determine amplicon type
        amplicon = "16S" if "_16S-" in fn else ("ITS" if "_ITS-" in fn else "unknown")

        # Extract descriptive part: Amplicon_{type}-{num}-{desc}_S{xx}_L001_R1_raw.fastq.gz
        m = re.search(r'Amplicon_(16S|ITS)-\d+-(.+?)_S\d+_L001', fn)
        if not m:
            unparseable.append(fn)
            continue

        desc = m.group(2)
        amplicon = m.group(1)

        # Parse flight/ground
        if re.search(r'VEG-05-F-|VEG-05F-', desc):
            flight = "Flight"
        elif re.search(r'VEG-05-G-', desc):
            flight = "Ground"
        else:
            flight = "unknown"

        # Parse plant
        pm = re.search(r'SN-?(\d+)', desc)
        plant = f"SN{pm.group(1)}" if pm else "NA"

        # Parse tissue/compartment
        if re.search(r'AdvRoot', desc, re.I):
            compartment = "AdvRoot"
        elif re.search(r'(?<!\w)root(?!s)', desc, re.I):
            compartment = "root"
        elif re.search(r'leaf', desc, re.I):
            compartment = "leaf"
        elif re.search(r'wick', desc, re.I):
            compartment = "wick"
        elif re.search(r'soil', desc, re.I):
            # soil-1, soil-2, soil-3
            sm = re.search(r'soil-(\d+)', desc)
            compartment = "soil"
            soil_rep = sm.group(1) if sm else "1"
        elif re.search(r'water', desc, re.I):
            compartment = "water"
        elif re.search(r'swab', desc, re.I):
            compartment = "swab"
        elif re.search(r'fruit|grn-frt', desc, re.I):
            compartment = "fruit"
        else:
            compartment = "unknown"

        # Parse light treatment
        if re.search(r'blue', desc, re.I):
            light = "Blue"
        elif re.search(r'red', desc, re.I):
            light = "Red"
        else:
            light = "NA"

        # Soil replicate
        soil_rep = ""
        if compartment == "soil":
            sm = re.search(r'soil-(\d+)', desc)
            soil_rep = sm.group(1) if sm else "1"

        # Fruit harvest info
        harvest = ""
        fruit_num = ""
        if compartment == "fruit":
            hm = re.search(r'Hrvst-(\d+)', desc)
            harvest = hm.group(1) if hm else ""
            fm = re.search(r'Fruit-(\d+)', desc)
            fruit_num = fm.group(1) if fm else ""

        # Swab location
        swab_loc = ""
        if compartment == "swab":
            if re.search(r'Locker', desc, re.I):
                sm = re.search(r'Locker-(\d+)-(\d+)', desc)
                swab_loc = f"Locker-{sm.group(1)}-{sm.group(2)}" if sm else "Locker"
            elif re.search(r'pillow', desc, re.I):
                swab_loc = "pillow-exterior"
            elif re.search(r'Blue-Swab|Red-Swab', desc):
                swab_loc = "chamber"

        # Build sample ID (without R1/R2 and _raw)
        sample_id = re.sub(r'_R1_raw\.fastq\.gz$', '', fn)
        sample_id = re.sub(r'^GLDS-672_Amplicon_(16S|ITS)-\d+-', '', sample_id)

        samples.append({
            "sample_id": sample_id,
            "study": "OSD-766",
            "amplicon": amplicon,
            "flight": flight,
            "plant": plant,
            "compartment": compartment,
            "light": light,
            "treatment": f"{flight}_{light}" if light != "NA" else flight,
            "soil_replicate": soil_rep,
            "harvest": harvest,
            "fruit_number": fruit_num,
            "swab_location": swab_loc,
            "fastq_r1": fn,
            "fastq_r2": fn.replace("_R1_", "_R2_"),
            "description": desc,
        })

    print(f"  Parsed: {len(samples)} samples")
    print(f"  Unparseable: {len(unparseable)}")
    if unparseable:
        for u in unparseable[:5]:
            print(f"    {u}")

    print(f"  Amplicon: {dict(Counter(s['amplicon'] for s in samples))}")
    print(f"  Flight/Ground: {dict(Counter(s['flight'] for s in samples))}")
    print(f"  Compartment: {dict(Counter(s['compartment'] for s in samples))}")
    print(f"  Light: {dict(Counter(s['light'] for s in samples))}")

    out_path = os.path.join(META_DIR, "sample_metadata_microbiome.csv")
    fields = ["sample_id", "study", "amplicon", "flight", "plant", "compartment",
              "light", "treatment", "soil_replicate", "harvest", "fruit_number",
              "swab_location", "fastq_r1", "fastq_r2", "description"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(samples)
    print(f"  Written: {out_path}")
    return samples


def build_crosswalk(rnaseq_samples, micro_samples):
    """Build crosswalk matching RNA-seq to microbiome samples."""
    print("\n=== Building sample crosswalk ===")

    # Integration tissues: Leaf and AdvRoot (present in both datasets)
    # Microbiome compartments that match RNA-seq tissues
    tissue_map = {"Leaf": "leaf", "Adv-Root": "AdvRoot"}

    crosswalk = []

    for rna in rnaseq_samples:
        rna_tissue = rna["tissue"]
        micro_compartment = tissue_map.get(rna_tissue)
        if not micro_compartment:
            continue

        # Find matching microbiome samples (same flight, plant, light, compartment)
        # Match on both 16S and ITS
        for amp in ["16S", "ITS"]:
            matches = [m for m in micro_samples
                       if m["flight"] == rna["flight"]
                       and m["plant"] == rna["plant"]
                       and m["light"] == rna["light"]
                       and m["compartment"] == micro_compartment
                       and m["amplicon"] == amp]
            for match in matches:
                crosswalk.append({
                    "rnaseq_sample_id": rna["sample_id"],
                    "rnaseq_tissue": rna_tissue,
                    "rnaseq_flight": rna["flight"],
                    "rnaseq_plant": rna["plant"],
                    "rnaseq_light": rna["light"],
                    "rnaseq_treatment": rna["treatment"],
                    "microbiome_sample_id": match["sample_id"],
                    "microbiome_amplicon": amp,
                    "microbiome_compartment": micro_compartment,
                    "microbiome_flight": match["flight"],
                    "microbiome_plant": match["plant"],
                    "microbiome_light": match["light"],
                    "match_key": f"{rna['flight']}_{rna['plant']}_{micro_compartment}_{rna['light']}_{amp}",
                })

    print(f"  Crosswalk entries: {len(crosswalk)}")
    print(f"  Unique RNA-seq samples with matches: {len(set(c['rnaseq_sample_id'] for c in crosswalk))}")
    print(f"  Unique microbiome samples matched: {len(set(c['microbiome_sample_id'] for c in crosswalk))}")

    # Summary by tissue and amplicon
    by_tissue = Counter((c["rnaseq_tissue"], c["microbiome_amplicon"]) for c in crosswalk)
    print(f"  By tissue x amplicon: {dict(by_tissue)}")

    # Check for RNA-seq samples without microbiome matches
    rna_matched = set(c["rnaseq_sample_id"] for c in crosswalk)
    rna_unmatched = [s for s in rnaseq_samples if s["sample_id"] not in rna_matched and s["tissue"] in tissue_map]
    if rna_unmatched:
        print(f"\n  RNA-seq samples WITHOUT microbiome match ({len(rna_unmatched)}):")
        for s in rna_unmatched:
            print(f"    {s['sample_id']} ({s['flight']}, {s['plant']}, {s['tissue']}, {s['light']})")

    out_path = os.path.join(META_DIR, "sample_crosswalk.csv")
    fields = ["rnaseq_sample_id", "rnaseq_tissue", "rnaseq_flight", "rnaseq_plant",
              "rnaseq_light", "rnaseq_treatment", "microbiome_sample_id",
              "microbiome_amplicon", "microbiome_compartment", "microbiome_flight",
              "microbiome_plant", "microbiome_light", "match_key"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(crosswalk)
    print(f"  Written: {out_path}")

    return crosswalk


def main():
    os.makedirs(META_DIR, exist_ok=True)
    rnaseq_samples = parse_rnaseq_metadata()
    micro_samples = parse_microbiome_metadata()
    crosswalk = build_crosswalk(rnaseq_samples, micro_samples)
    print("\n=== Metadata parsing complete ===")


if __name__ == "__main__":
    main()
