#!/usr/bin/env python3
import pysam
import re
from collections import defaultdict
import argparse

def parse_cigar(cigar_string):
    """Parse CIGAR string to find soft-clipped regions"""
    # Find soft clips at start and end
    start_clip = re.match(r'^(\d+)S', cigar_string)
    end_clip = re.search(r'(\d+)S$', cigar_string)
    
    start_clip_len = int(start_clip.group(1)) if start_clip else 0
    end_clip_len = int(end_clip.group(1)) if end_clip else 0
    
    return start_clip_len, end_clip_len

def extract_softclipped_sequence(read, start_clip_len, end_clip_len):
    """Extract the soft-clipped portions of the sequence"""
    seq = read.query_sequence
    if not seq:
        return None, None
    
    start_seq = seq[:start_clip_len] if start_clip_len > 0 else None
    end_seq = seq[-end_clip_len:] if end_clip_len > 0 else None
    
    return start_seq, end_seq

def analyze_softclipped_reads(bam_file, min_clip_length=10, output_file="fusion_candidates.txt"):
    """Analyze soft-clipped reads for potential fusion events"""
    
    # Define canonical chromosomes (adjust based on your reference)
    canonical_chromosomes = {
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY',
        # Also include versions without 'chr' prefix
        '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
        '21', '22', 'X', 'Y'
    }
    
    # Dictionary to store potential fusion partners
    fusion_candidates = defaultdict(list)
    read_pairs = defaultdict(list)
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
                
            cigar_string = read.cigarstring
            if not cigar_string or 'S' not in cigar_string:
                continue
            
            start_clip, end_clip = parse_cigar(cigar_string)
            
            # Only consider reads with significant soft-clipping AND on canonical chromosomes
            if (start_clip >= min_clip_length or end_clip >= min_clip_length) and read.reference_name in canonical_chromosomes:
                start_seq, end_seq = extract_softclipped_sequence(read, start_clip, end_clip)
                
                read_info = {
                    'read_name': read.query_name,
                    'chromosome': read.reference_name,
                    'position': read.reference_start,
                    'cigar': cigar_string,
                    'start_clip_len': start_clip,
                    'end_clip_len': end_clip,
                    'start_clip_seq': start_seq,
                    'end_clip_seq': end_seq,
                    'mapping_quality': read.mapping_quality
                }
                
                # Store by read name to find pairs
                read_pairs[read.query_name].append(read_info)
    
    # Group similar fusion events and count supporting reads
    fusion_events = defaultdict(list)
    
    for read_name, reads in read_pairs.items():
        if len(reads) == 2:  # Paired reads
            read1, read2 = reads[0], reads[1]
            
            # Check that both reads are on canonical chromosomes
            if read1['chromosome'] in canonical_chromosomes and read2['chromosome'] in canonical_chromosomes:
                
                # Check if reads map to different chromosomes (inter-chromosomal)
                if read1['chromosome'] != read2['chromosome']:
                    # Create a fusion event key (sort chromosomes for consistency)
                    chr1, pos1, chr2, pos2 = read1['chromosome'], read1['position'], read2['chromosome'], read2['position']
                    if chr1 > chr2:  # Sort to ensure consistent grouping
                        chr1, pos1, chr2, pos2 = chr2, pos2, chr1, pos1
                    
                    # Group by approximate breakpoint (within 1000bp window)
                    breakpoint_key = f"{chr1}:{pos1//1000}_{chr2}:{pos2//1000}"
                    
                    fusion_events[breakpoint_key].append({
                        'read_name': read_name,
                        'chr1': chr1,
                        'pos1': pos1,
                        'chr2': chr2,
                        'pos2': pos2,
                        'type': 'inter-chromosomal'
                    })
                
                # Check for intra-chromosomal but >100MB apart
                elif abs(read1['position'] - read2['position']) > 100000:  # >100K apart
                    chr1 = read1['chromosome']
                    pos1, pos2 = sorted([read1['position'], read2['position']])  # Sort positions
                    
                    # Group by approximate breakpoint (within 1000bp window)
                    breakpoint_key = f"{chr1}:{pos1//1000}_{chr1}:{pos2//1000}"
                    
                    fusion_events[breakpoint_key].append({
                        'read_name': read_name,
                        'chr1': chr1,
                        'pos1': pos1,
                        'chr2': chr1,
                        'pos2': pos2,
                        'type': 'intra-chromosomal',
                        'distance': abs(pos2 - pos1)
                    })
    
    # Filter fusion events with at least 2 supporting reads
    potential_fusions = []
    for breakpoint_key, supporting_reads in fusion_events.items():
        if len(supporting_reads) >= 10:  # At least 10 supporting reads
            # Use the first read as representative for the fusion event
            representative = supporting_reads[0]
            representative['supporting_reads'] = len(supporting_reads)
            representative['read_names'] = [r['read_name'] for r in supporting_reads]
            potential_fusions.append(representative)
    
    # Write results
    with open(output_file, 'w') as f:
        f.write("Chr1\tPos1\tChr2\tPos2\tType\tDistance\tSupportingReads\tReadNames\n")
        for fusion in potential_fusions:
            distance = fusion.get('distance', 'N/A')
            read_names = ','.join(fusion['read_names'][:5])  # Show first 5 read names
            if len(fusion['read_names']) > 5:
                read_names += f",... (+{len(fusion['read_names'])-5} more)"
            
            f.write(f"{fusion['chr1']}\t{fusion['pos1']}\t"
                   f"{fusion['chr2']}\t{fusion['pos2']}\t{fusion['type']}\t{distance}\t"
                   f"{fusion['supporting_reads']}\t{read_names}\n")
    
    print(f"Found {len(potential_fusions)} potential fusion candidates with ≥10 supporting reads")
    print(f"Results written to {output_file}")
    
    return potential_fusions

def main():
    parser = argparse.ArgumentParser(description='Detect potential fusions from soft-clipped reads with minimum read support')
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('--min-clip-length', type=int, default=10, 
                       help='Minimum soft-clip length to consider (default: 10)')
    parser.add_argument('--min-reads', type=int, default=5,
                       help='Minimum number of supporting reads for a fusion event (default: 10)')
    parser.add_argument('--output', default='fusion_candidates.txt',
                       help='Output file name (default: fusion_candidates.txt)')
    
    args = parser.parse_args()
    
    analyze_softclipped_reads(args.bam_file, args.min_clip_length, args.output)

if __name__ == "__main__":
    main()
