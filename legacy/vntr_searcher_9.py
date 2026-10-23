#!/usr/bin/env python3

import argparse
import subprocess
import os
import glob
from pathlib import Path
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Search for VNTR regions using mePCR')
    parser.add_argument('--sts', required=True, help='Path to the STS primers file')
    parser.add_argument('--fasta', required=True, help='Path pattern to FASTA file(s)')
    parser.add_argument('--outdir', required=True, help='Output directory for results')
    parser.add_argument('--prefix', required=True, help='File prefix')
    parser.add_argument('-M', default=500, type=int, help='Margin (default: 500)')
    parser.add_argument('-N', default=3, type=int, help='Number of mismatches allowed (default: 3)')
    parser.add_argument('--clwclean', action='store_true', help='Create additional CLW files without sequences containing N')
    parser.add_argument('--pattern', help='Pattern file to use for count summary generation')
    return parser.parse_args()

def run_mepcr(sts_file, fasta_file, max_size, mismatches):
    """Run mePCR command for a single FASTA file and return results."""
    try:
        cmd = [
            'me-PCR',
            sts_file,
            fasta_file,
            f'M={max_size}',
            f'N={mismatches}'
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        return parse_mepcr_output(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error processing {fasta_file}: {e}")
        print(f"Standard error: {e.stderr}")
        return None
    except Exception as e:
        print(f"Unexpected error processing {fasta_file}: {e}")
        return None

def parse_mepcr_output(output_text):
    """Parse mePCR output and return structured data."""
    results = []
    for line in output_text.strip().split('\n'):
        if line:
            parts = line.strip().split('\t')
            # Ensure we have at least the required columns
            if len(parts) >= 3:
                result = {
                    'contig': parts[0],
                    'position': parts[1],
                    'primer_id': parts[2],
                    'strand': parts[-1].strip('()')
                }
                
                # Parse the position
                start, end = map(int, parts[1].split('..'))
                result['start'] = start
                result['end'] = end
                
                results.append(result)
    return results

def extract_sequence(fasta_file, contig, start, end, strand):
    """Extract sequence from FASTA file for given coordinates."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == contig:
            seq = record.seq[start-1:end]
            if strand == '-':
                seq = seq.reverse_complement()
            return seq
    return None

def filter_N_sequences(records):
    """Filter out sequences containing N."""
    clean_records = []
    filtered_count = 0
    
    for record in records:
        if 'N' not in str(record.seq).upper():
            clean_records.append(record)
        else:
            filtered_count += 1
    
    return clean_records, filtered_count

def process_alignment(records, input_file, align_file, clustal_file, primer_id, clean=False):
    """Process alignment and conversion to CLUSTAL format."""
    prefix = "clean " if clean else ""
    
    if len(records) > 1:
        try:
            # First try MUSCLE v5 syntax
            cmd = ['muscle', '-align', input_file, '-output', align_file]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            # If v5 fails, try MUSCLE v3 syntax
            if result.returncode != 0:
                cmd = ['muscle', '-in', input_file, '-out', align_file, '-quiet']
                result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"Created {prefix}alignment for {primer_id}")
                
                # Convert to CLUSTAL format using seqret
                seqret_cmd = [
                    'seqret',
                    '-sequence', align_file,
                    '-outseq', clustal_file,
                    '-osformat', 'clustal'
                ]
                seqret_result = subprocess.run(seqret_cmd, capture_output=True, text=True)
                
                if seqret_result.returncode == 0:
                    print(f"Created {prefix}CLUSTAL format for {primer_id}")
                else:
                    print(f"Error creating {prefix}CLUSTAL format for {primer_id}: {seqret_result.stderr}")
            else:
                print(f"Error creating {prefix}alignment for {primer_id}: {result.stderr}")
        except Exception as e:
            print(f"Error in {prefix}alignment process for {primer_id}: {e}")
    else:
        # If only one sequence, copy the unaligned file as the aligned files
        import shutil
        shutil.copy(input_file, align_file)
        
        # Create CLUSTAL format for single sequence
        seqret_cmd = [
            'seqret',
            '-sequence', input_file,
            '-outseq', clustal_file,
            '-osformat', 'clustal'
        ]
        try:
            subprocess.run(seqret_cmd, capture_output=True, text=True, check=True)
            print(f"Single sequence for {prefix}{primer_id} - converted to CLUSTAL format")
        except subprocess.CalledProcessError as e:
            print(f"Error creating {prefix}CLUSTAL format for single sequence {primer_id}: {e.stderr}")

def write_sequences_to_files(results_by_primer, outdir, clwclean=False):
    """Write sequences to separate FASTA files for each primer pair."""
    primer_files = {}
    
    # Create a directory for sequence files
    seq_dir = os.path.join(outdir, "sequences")
    os.makedirs(seq_dir, exist_ok=True)
    
    # Write sequences for each primer pair
    for primer_id, primer_results in results_by_primer.items():
        fasta_file = os.path.join(seq_dir, f"{primer_id}.fasta")
        align_file = os.path.join(seq_dir, f"{primer_id}_aligned.fasta")
        clustal_file = os.path.join(seq_dir, f"{primer_id}_aligned.clw")
        
        if clwclean:
            clean_fasta = os.path.join(seq_dir, f"{primer_id}_clean.fasta")
            clean_align = os.path.join(seq_dir, f"{primer_id}_clean_aligned.fasta")
            clean_clustal = os.path.join(seq_dir, f"{primer_id}_clean_aligned.clw")
        
        records = []
        for result in primer_results:
            header = f"{result['fasta_name']}_{result['count']}"
            seq_record = SeqRecord(
                Seq(str(result['sequence'])),
                id=header,
                description=""
            )
            records.append(seq_record)
        
        # Write unaligned sequences
        SeqIO.write(records, fasta_file, "fasta")
        
        # Process regular alignment
        process_alignment(records, fasta_file, align_file, clustal_file, primer_id)
        
        # Process clean alignment if requested
        if clwclean and records:
            clean_records, filtered_count = filter_N_sequences(records)
            if clean_records:
                SeqIO.write(clean_records, clean_fasta, "fasta")
                print(f"Filtered {filtered_count} sequences containing N from {primer_id}")
                process_alignment(clean_records, clean_fasta, clean_align, clean_clustal, primer_id, clean=True)
            else:
                print(f"No clean sequences (without N) found for {primer_id}")

def create_summary_table(results_by_primer, fasta_files, output_prefix, outdir, clwclean=False):
    """Create summary tables of sequence lengths for both regular and clean sequences."""
    # Initialize dictionaries to store the lengths
    length_data = defaultdict(dict)
    clean_length_data = defaultdict(dict) if clwclean else None
    
    # Get all unique fasta file names
    fasta_names = [Path(f).stem for f in fasta_files]
    
    # Process each primer's results
    for primer_id, results in results_by_primer.items():
        # Initialize all fasta files with None for this primer
        for fasta_name in fasta_names:
            length_data[fasta_name][primer_id] = None
            if clwclean:
                clean_length_data[fasta_name][primer_id] = None
        
        # Fill in the lengths where we have data
        for result in results:
            fasta_name = result['fasta_name']
            sequence = str(result['sequence'])
            new_length = len(sequence)
            
            # Regular lengths
            current_value = length_data[fasta_name][primer_id]
            if current_value is None:
                length_data[fasta_name][primer_id] = str(new_length)
            else:
                length_data[fasta_name][primer_id] = f"{current_value}/{new_length}"
            
            # Clean lengths (sequences without N)
            if clwclean and 'N' not in sequence.upper():
                current_clean_value = clean_length_data[fasta_name][primer_id]
                if current_clean_value is None:
                    clean_length_data[fasta_name][primer_id] = str(new_length)
                else:
                    clean_length_data[fasta_name][primer_id] = f"{current_clean_value}/{new_length}"
    
    # Convert to DataFrames
    df = pd.DataFrame.from_dict(length_data, orient='index')
    df.sort_index(inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    
    # Save regular lengths
    output_file = os.path.join(outdir, f"{output_prefix}_lengths_summary.tsv")
    df.to_csv(output_file, sep='\t')
    print(f"\nRegular sequence lengths summary saved to: {output_file}")
    
    # Process and save clean lengths if requested
    if clwclean:
        clean_df = pd.DataFrame.from_dict(clean_length_data, orient='index')
        clean_df.sort_index(inplace=True)
        clean_df = clean_df.reindex(sorted(clean_df.columns), axis=1)
        
        clean_output_file = os.path.join(outdir, f"{output_prefix}_lengths_summary_clean.tsv")
        clean_df.to_csv(clean_output_file, sep='\t')
        print(f"Clean sequence lengths summary saved to: {clean_output_file}")
        
        return df, clean_df
    
    return df, None

def read_pattern_file(pattern_file):
    """
    Read patterns from a TSV file. Support multiple patterns per marker separated by semicolons.
    
    Args:
        pattern_file: Path to pattern file
        
    Returns:
        List of pattern results
    """
    results = []
    
    try:
        with open(pattern_file, 'r') as f:
            # Skip header
            header = f.readline().strip()
            
            for line in f:
                # Skip lines starting with #
                if not line.startswith("#"):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:  # Need at least marker and pattern
                        primer_id = parts[0]
                        pattern_str = parts[1]
                        
                        # Check if the pattern is valid
                        if pattern_str == "N/A":
                            results.append({
                                "primer_id": primer_id,
                                "has_repeats": False,
                                "reason": "No pattern in file"
                            })
                            continue
                        
                        # Split pattern string by semicolons to get multiple patterns
                        patterns = [p.strip() for p in pattern_str.split(';')]
                        
                        # Try to parse pattern length, min repeats, max repeats, and score
                        try:
                            # For multiple patterns, we'll use the first pattern's length as default
                            pattern_length = int(parts[2]) if len(parts) > 2 else len(patterns[0])
                            min_repeats = int(parts[3]) if len(parts) > 3 else 1
                            max_repeats = int(parts[4]) if len(parts) > 4 else 10
                            score = float(parts[5]) if len(parts) > 5 else 1.0
                            
                            results.append({
                                "primer_id": primer_id,
                                "has_repeats": True,
                                "patterns": patterns,  # Now storing a list of patterns
                                "pattern_length": pattern_length,
                                "min_count": min_repeats,
                                "max_count": max_repeats,
                                "score": score
                            })
                        except (ValueError, IndexError):
                            # If there's an error, still try to use the patterns
                            results.append({
                                "primer_id": primer_id,
                                "has_repeats": True,
                                "patterns": patterns,  # Now storing a list of patterns
                                "pattern_length": len(patterns[0]),
                                "min_count": 1,
                                "max_count": 10,
                                "score": 1.0
                            })
        
        print(f"Loaded {len(results)} pattern definitions from {pattern_file}")
    except Exception as e:
        print(f"Error reading pattern file {pattern_file}: {e}")
    
    return results

def count_patterns_in_fasta(pattern_results, fasta_files, output_prefix, outdir, is_clean=False):
    """
    Count patterns in extracted FASTA files for each marker and sample.
    Shows individual counts for each entry separated by ";".
    Handles multiple patterns per marker, summing their occurrences.
    Missing samples are represented as "NA".
    Pattern matching is case-insensitive.
    
    Args:
        pattern_results: List of pattern analysis results
        fasta_files: List of original FASTA file paths (for sample names)
        output_prefix: Prefix for output files
        outdir: Output directory
        is_clean: Whether this is for clean sequences
        
    Returns:
        DataFrame with pattern counts
    """
    # Get sample names from FASTA file names
    sample_names = [Path(f).stem for f in fasta_files]
    print(f"Found {len(sample_names)} sample names from FASTA files")
    
    # Initialize a dictionary to store pattern counts
    # Structure: {sample_name: {marker_id: counts}}
    pattern_counts = defaultdict(dict)
    
    # Initialize with None for ALL sample/marker combinations
    valid_markers = [r["primer_id"] for r in pattern_results if r["has_repeats"]]
    for sample_name in sample_names:
        for marker_id in valid_markers:
            pattern_counts[sample_name][marker_id] = None
    
    # Organize patterns by marker ID for quick lookup and convert to lowercase
    patterns_by_marker = {r["primer_id"]: [p.lower() for p in r["patterns"]] for r in pattern_results if r["has_repeats"]}
    
    # Process sequences directory
    seq_dir = os.path.join(outdir, "sequences")
    
    # For each marker with valid patterns
    for marker_id, patterns in patterns_by_marker.items():
        # Look for the marker-specific FASTA file
        marker_fasta = os.path.join(seq_dir, f"{marker_id}.fasta")
        if not os.path.exists(marker_fasta):
            print(f"  No FASTA file found for marker {marker_id}")
            continue
            
        print(f"Processing marker {marker_id}")
        
        # Read the sequences for this marker
        try:
            records = list(SeqIO.parse(marker_fasta, "fasta"))
            
            # Group records by sample
            sample_records = defaultdict(list)
            
            for record in records:
                # Extract sample name from record ID by removing the last part after the last underscore
                parts = record.id.rsplit('_', 1)
                sample_name = parts[0]  # Everything before the last underscore
                
                # For clean sequences, skip if it has 'N'
                if is_clean and 'N' in str(record.seq).upper():
                    continue
                    
                sample_records[sample_name].append(record)
            
            # Process each sample's records
            for sample_name, records in sample_records.items():
                # Count pattern occurrences in each sequence separately
                counts = []
                
                for record in records:
                    # Convert sequence to lowercase for case-insensitive matching
                    seq = str(record.seq).lower()
                    
                    # Count occurrences of all patterns and sum them
                    total_count = 0
                    
                    for pattern in patterns:
                        # Count non-overlapping occurrences for this pattern
                        pos = 0
                        while True:
                            pos = seq.find(pattern, pos)
                            if pos == -1:
                                break
                            total_count += 1
                            pos += len(pattern)  # Move past this occurrence
                    
                    counts.append(str(total_count))
                
                # Join the counts with ";" 
                if counts:
                    pattern_counts[sample_name][marker_id] = ";".join(counts)
                
        except Exception as e:
            print(f"  Error counting patterns in {marker_fasta}: {e}")
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(pattern_counts, orient='index')
    
    # Sort rows and columns
    df.sort_index(inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    
    # Save to file
    suffix = "_clean" if is_clean else ""
    output_file = os.path.join(outdir, f"{output_prefix}_count_summary{suffix}.tsv")
    
    # Write to file with explicit NA representation
    df.to_csv(output_file, sep='\t', na_rep="NA")
    
    print(f"Pattern counts summary saved to: {output_file}")
    print(f"  Wrote data for {len(df.index)} samples and {len(df.columns)} markers")
    
    return df

def count_length_differences_in_fasta(pattern_results, fasta_files, output_prefix, outdir, is_clean=False):
    """
    Count length differences relative to the longest sequence for each marker and sample.
    Calculates (longest_length - current_length) / pattern_length for each sequence.
    Shows individual calculations for each entry separated by ";".
    Missing samples are represented as "NA".
    
    Args:
        pattern_results: List of pattern analysis results
        fasta_files: List of original FASTA file paths (for sample names)
        output_prefix: Prefix for output files
        outdir: Output directory
        is_clean: Whether this is for clean sequences
        
    Returns:
        DataFrame with length difference calculations
    """
    # Get sample names from FASTA file names
    sample_names = [Path(f).stem for f in fasta_files]
    print(f"Found {len(sample_names)} sample names from FASTA files")
    
    # Initialize a dictionary to store length difference calculations
    # Structure: {sample_name: {marker_id: calculations}}
    length_diff_data = defaultdict(dict)
    
    # Initialize with None for ALL sample/marker combinations
    valid_markers = [r["primer_id"] for r in pattern_results if r["has_repeats"]]
    for sample_name in sample_names:
        for marker_id in valid_markers:
            length_diff_data[sample_name][marker_id] = None
    
    # Organize pattern lengths by marker ID for quick lookup
    pattern_lengths = {r["primer_id"]: r["pattern_length"] for r in pattern_results if r["has_repeats"]}
    
    # Process sequences directory
    seq_dir = os.path.join(outdir, "sequences")
    
    # For each marker with valid patterns
    for marker_id in valid_markers:
        if marker_id not in pattern_lengths:
            print(f"  No pattern length found for marker {marker_id}")
            continue
            
        pattern_length = pattern_lengths[marker_id]
        if pattern_length <= 0:
            print(f"  Invalid pattern length ({pattern_length}) for marker {marker_id}")
            continue
        
        # Look for the marker-specific FASTA file
        marker_fasta = os.path.join(seq_dir, f"{marker_id}.fasta")
        if not os.path.exists(marker_fasta):
            print(f"  No FASTA file found for marker {marker_id}")
            continue
            
        print(f"Processing length differences for marker {marker_id} (pattern length: {pattern_length})")
        
        # Read the sequences for this marker
        try:
            records = list(SeqIO.parse(marker_fasta, "fasta"))
            
            # Group records by sample and collect all sequence lengths
            sample_records = defaultdict(list)
            all_lengths = []
            
            for record in records:
                # Extract sample name from record ID by removing the last part after the last underscore
                parts = record.id.rsplit('_', 1)
                sample_name = parts[0]  # Everything before the last underscore
                
                # For clean sequences, skip if it has 'N'
                if is_clean and 'N' in str(record.seq).upper():
                    continue
                    
                seq_length = len(record.seq)
                sample_records[sample_name].append((record, seq_length))
                all_lengths.append(seq_length)
            
            # Find the shortest sequence length across all samples for this marker
            if not all_lengths:
                print(f"  No valid sequences found for marker {marker_id}")
                continue
                
            shortest_length = min(all_lengths)
            print(f"  Shortest sequence length for {marker_id}: {shortest_length}")
            
            # Process each sample's records
            for sample_name, sample_data in sample_records.items():
                # Calculate length differences for each sequence
                calculations = []
                
                for record, seq_length in sample_data:
                    # Calculate (current_length-shortest_length) / pattern_length
                    length_diff = seq_length-shortest_length
                    calculation = length_diff / pattern_length
                    calculations.append(f"{calculation:.3f}")
                
                # Join the calculations with ";" 
                if calculations:
                    length_diff_data[sample_name][marker_id] = ";".join(calculations)
                
        except Exception as e:
            print(f"  Error calculating length differences in {marker_fasta}: {e}")
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(length_diff_data, orient='index')
    
    # Sort rows and columns
    df.sort_index(inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    
    # Save to file
    suffix = "_clean" if is_clean else ""
    output_file = os.path.join(outdir, f"{output_prefix}_count_lengthdiff{suffix}.tsv")
    
    # Write to file with explicit NA representation
    df.to_csv(output_file, sep='\t', na_rep="NA")
    
    print(f"Length difference summary saved to: {output_file}")
    print(f"  Wrote data for {len(df.index)} samples and {len(df.columns)} markers")
    
    return df

def main():
    """Main function to process multiple FASTA files with mePCR."""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    
    # Check if sequences directory already exists
    sequences_dir = os.path.join(args.outdir, "sequences")
    should_run_mepcr = not os.path.exists(sequences_dir)  # Renamed from run_mepcr to should_run_mepcr
    
    # Get list of FASTA files
    fasta_files = glob.glob(args.fasta)
    if not fasta_files:
        print(f"No FASTA files found matching pattern: {args.fasta}")
        return
    
    # Dictionary to store results by primer pair
    results_by_primer = defaultdict(list)
    
    if should_run_mepcr:  # Also changed here
        # Process each FASTA file
        successful = 0
        failed = 0
        
        print("Sequence directory not found. Running full analysis with mePCR...")
        
        for fasta_file in fasta_files:
            print(f"Processing: {fasta_file}")
            fasta_name = Path(fasta_file).stem
            
            # Run mePCR
            results = run_mepcr(args.sts, fasta_file, args.M, args.N)
            
            if results:
                successful += 1
                
                # Count occurrences for each primer in this file
                primer_counts = defaultdict(int)
                
                # Process results
                for result in results:
                    primer_id = result['primer_id']
                    primer_counts[primer_id] += 1
                    
                    # Extract sequence
                    sequence = extract_sequence(
                        fasta_file,
                        result['contig'],
                        result['start'],
                        result['end'],
                        result['strand']
                    )
                    
                    if sequence:
                        result['sequence'] = sequence
                        result['fasta_name'] = fasta_name
                        result['count'] = primer_counts[primer_id]
                        results_by_primer[primer_id].append(result)
            else:
                failed += 1
        
        # Write sequences and create alignments
        write_sequences_to_files(results_by_primer, args.outdir, clwclean=args.clwclean)
        
        print("\nProcessing complete!")
        print(f"Successfully processed: {successful} files")
        print(f"Failed to process: {failed} files")
    else:
        print(f"Sequence directory found at {sequences_dir}")
        print("Skipping mePCR processing and using existing sequence files.")
        
        # Load information about existing primer files
        for primer_file in glob.glob(os.path.join(sequences_dir, "*.fasta")):
            if "_aligned" not in primer_file and "_clean" not in primer_file:
                primer_id = Path(primer_file).stem
                print(f"Found existing primer file: {primer_id}")
                
                # Count entries for each sample in this primer file
                try:
                    for record in SeqIO.parse(primer_file, "fasta"):
                        # Extract sample name from record ID by removing the last part after the last underscore
                        parts = record.id.rsplit('_', 1)
                        sample_name = parts[0]
                        
                        # Create a dummy result entry with the necessary fields
                        result = {
                            'primer_id': primer_id,
                            'fasta_name': sample_name,
                            'sequence': record.seq,
                            'count': 1  # This might not be accurate, but it's not used for length summary
                        }
                        
                        results_by_primer[primer_id].append(result)
                except Exception as e:
                    print(f"Error reading primer file {primer_file}: {e}")
    
    # Create and save summary tables
    print("\nGenerating length summary tables...")
    summary_df, clean_summary_df = create_summary_table(
        results_by_primer, 
        fasta_files,
        args.prefix,
        args.outdir,
        clwclean=args.clwclean
    )
    
    # If pattern file is provided, generate count summaries
    if args.pattern:
        print(f"\nUsing provided pattern file: {args.pattern}")
        pattern_results = read_pattern_file(args.pattern)
        
        # Generate standard count summary
        print("\nGenerating pattern count summary...")
        count_patterns_in_fasta(pattern_results, fasta_files, args.prefix, args.outdir, is_clean=False)
        
        # Generate length difference summary
        print("\nGenerating length difference summary...")
        count_length_differences_in_fasta(pattern_results, fasta_files, args.prefix, args.outdir, is_clean=False)
        
        # If clean analysis requested, process clean sequences too
        if args.clwclean:
            print("\nGenerating clean pattern count summary...")
            count_patterns_in_fasta(pattern_results, fasta_files, args.prefix, args.outdir, is_clean=True)
            
            print("\nGenerating clean length difference summary...")
            count_length_differences_in_fasta(pattern_results, fasta_files, args.prefix, args.outdir, is_clean=True)
    else:
        print("\nNo pattern file provided. Skipping pattern count summary generation.")
    
    # Print summary of outputs
    print(f"\nResults saved in: {args.outdir}")
    if should_run_mepcr:  # Changed here as well
        print(f"Sequence files and alignments saved in: {sequences_dir}")
    print("Summary files saved:")
    print(f"- {args.prefix}_lengths_summary.tsv")
    
    if args.clwclean:
        print(f"- {args.prefix}_lengths_summary_clean.tsv")
    
    if args.pattern:
        print(f"- {args.prefix}_count_summary.tsv")
        print(f"- {args.prefix}_count_lengthdiff.tsv")
        if args.clwclean:
            print(f"- {args.prefix}_count_summary_clean.tsv")
            print(f"- {args.prefix}_count_lengthdiff_clean.tsv")

if __name__ == "__main__":
    main()