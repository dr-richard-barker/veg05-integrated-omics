#!/usr/bin/env python3
"""
FAPROTAX collapse: Map taxonomy to functional categories.
"""
import sys
import csv
import re
import os
from collections import defaultdict

def parse_faprotax_db(db_path):
    """Parse FAPROTAX.txt into function -> list of taxon patterns."""
    func_to_taxa = defaultdict(list)
    current_func = None
    
    with open(db_path) as f:
        for line in f:
            line = line.rstrip('\n')
            # Skip comments and empty lines
            if line.startswith('#'):
                continue
            if not line.strip():
                # Blank line separates groups
                current_func = None
                continue
            
            # Split on tab - but member taxa may have trailing tabs/comments
            # First, strip trailing comments
            if '#' in line:
                line = line[:line.index('#')]
            
            parts = line.split('\t')
            # Clean parts
            parts = [p.strip() for p in parts]
            
            # If first part starts with *, it's a taxon pattern
            if parts[0].startswith('*') or parts[0].startswith("'"):
                if current_func:
                    taxon = parts[0].strip("'\"")
                    func_to_taxa[current_func].append(taxon)
            elif parts[0].startswith('add_group:') or parts[0].startswith('subtract_group:') or parts[0].startswith('intersect_group:'):
                # Set operations - skip for simplicity
                if current_func:
                    func_to_taxa[current_func].append(parts[0])
            else:
                # It's a new function definition
                # The function name may have metadata after a tab
                current_func = parts[0]
    
    return func_to_taxa

def match_taxon(pattern, taxonomy):
    """Check if a taxon pattern matches a taxonomy string."""
    # Handle set operations
    if pattern.startswith('add_group:') or pattern.startswith('subtract_group:') or pattern.startswith('intersect_group:'):
        return False  # Skip set operations for now
    
    # Convert glob pattern to regex (* = .*)
    regex_pattern = re.escape(pattern).replace(r'\*', '.*')
    regex = re.compile(regex_pattern, re.IGNORECASE)
    return bool(regex.search(taxonomy))

def collapse_taxonomy(input_file, db_path, output_file):
    """Collapse OTU table by FAPROTAX functions."""
    func_to_taxa = parse_faprotax_db(db_path)
    print(f"Parsed {len(func_to_taxa)} FAPROTAX functions")
    
    # Read input
    with open(input_file) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        
        sample_names = header[1:]
        
        rows = list(reader)
        print(f"Read {len(rows)} taxa")
        
        # For each taxon, find matching functions
        func_counts = defaultdict(lambda: [0.0] * len(sample_names))
        func_taxa = defaultdict(list)
        matched_count = 0
        
        for row in rows:
            taxonomy = row[0]
            counts = [float(x) if x else 0 for x in row[1:]]
            
            matched_funcs = set()
            for func, patterns in func_to_taxa.items():
                for pattern in patterns:
                    if match_taxon(pattern, taxonomy):
                        matched_funcs.add(func)
                        break
            
            if matched_funcs:
                matched_count += 1
                for func in matched_funcs:
                    for i, c in enumerate(counts):
                        func_counts[func][i] += c
                    func_taxa[func].append(taxonomy)
        
        print(f"Matched {matched_count}/{len(rows)} taxa to at least one function")
        print(f"Total functions with matches: {len(func_counts)}")
        
        # Write output (function table)
        with open(output_file, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['function'] + sample_names)
            for func in sorted(func_counts.keys()):
                row = [func] + [str(int(c)) if c == int(c) else str(round(c, 2)) for c in func_counts[func]]
                writer.writerow(row)
        
        print(f"Wrote {len(func_counts)} functions to {output_file}")
        
        # Write report
        report_file = output_file.replace('.tsv', '_report.tsv')
        with open(report_file, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['function', 'n_taxa', 'total_reads', 'representative_taxa'])
            for func in sorted(func_counts.keys(), key=lambda x: sum(func_counts[x]), reverse=True):
                total = sum(func_counts[func])
                taxa_list = func_taxa[func][:3]
                writer.writerow([func, len(func_taxa[func]), int(total), '; '.join(taxa_list)])
        
        print(f"Wrote report to {report_file}")
        
        # Print top functions
        print("\nTop 20 functions by total reads:")
        sorted_funcs = sorted(func_counts.keys(), key=lambda x: sum(func_counts[x]), reverse=True)
        for func in sorted_funcs[:20]:
            total = sum(func_counts[func])
            print(f"  {func}: {int(total)} reads ({len(func_taxa[func])} taxa)")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='FAPROTAX collapse')
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-g', '--database', required=True)
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()
    
    collapse_taxonomy(args.input, args.database, args.output)
