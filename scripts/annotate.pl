#!/usr/bin/env perl
#
# Copyright (c) 2025 Abhinav Mishra
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Author: Abhinav Mishra <mishraabhinav36@gmail.com>
# Date: 2025-01-31

use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;
use List::Util qw(min max);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

# Global data structures for genomic features
my (%genes, %tss, %cpg, %repeats, %sorted_cpg, %sorted_repeats);

# Detect number of CPU cores for parallel processing
my $num_cores = `lscpu -p | grep -v '^#' | wc -l`;
chomp($num_cores);
my $pm = Parallel::ForkManager->new($num_cores);

##############################
# Data Loading Subroutines
##############################

### Load GFF3 Gene Annotations
# Accept features of type "gene" or "ncRNA_gene" and extract ID, biotype, and Name.
sub load_gene_annotation {
    my ($gff_file) = @_;
    my $fh;
    if ($gff_file =~ /\.gz$/) {
        $fh = IO::Uncompress::Gunzip->new($gff_file)
          or die "Cannot open gzipped GFF3 file ($gff_file): $GunzipError\n";
    } else {
        open($fh, "<", $gff_file) or die "Cannot open GFF3 file ($gff_file): $!\n";
    }
    
    while (<$fh>) {
        chomp;
        next if /^#/;
        my @fields = split "\t", $_;
        # Only consider features of type "gene" or "ncRNA_gene"
        next unless @fields >= 9 && ($fields[2] eq "gene" || $fields[2] eq "ncRNA_gene");
        
        my ($chr, $start, $end, $strand) = @fields[0,3,4,6];
        my ($gene_id)       = $fields[8] =~ /ID=([^;]+)/;
        my ($gene_type)     = $fields[8] =~ /biotype=([^;]+)/;
        my ($gene_synonyms) = $fields[8] =~ /Name=([^;]+)/;
        
        $genes{$chr}{$start} = { 
            end           => $end,
            strand        => $strand,
            gene_id       => $gene_id,
            gene_type     => ($gene_type // "Unknown"),
            gene_synonyms => ($gene_synonyms // "N/A")
        };
    }
    close $fh;
}

### Load TSS Data (skip header; expected columns:
### gene_id, transcript_id, start, end, chr, gene_type, gene_synonyms)
sub load_tss_data {
    my ($tss_file) = @_;
    my $fh;
    if ($tss_file =~ /\.gz$/) {
        $fh = IO::Uncompress::Gunzip->new($tss_file)
          or die "Cannot open gzipped TSS file ($tss_file): $GunzipError\n";
    } else {
        open($fh, "<", $tss_file) or die "Cannot open TSS file ($tss_file): $!\n";
    }
    my $header = <$fh>;  # skip header line
    while (<$fh>) {
        chomp;
        my @cols = split "\t", $_;
        next unless @cols >= 7;
        my ($gene_id, $transcript_id, $start, $end, $chr, $gene_type, $gene_synonyms) = @cols[0,1,2,3,4,5,6];
        $chr = ($chr =~ /^chr/) ? $chr : "chr$chr";
        $tss{$chr}{$start} = { 
            gene_id       => $gene_id,
            transcript_id => $transcript_id,
            start         => $start,
            end           => $end,
            gene_type     => $gene_type,
            gene_synonyms => ($gene_synonyms // "N/A")
        };
    }
    close $fh;
}

### Load CpG Island Data (no header; expected 11 columns):
### bin, chrom, chromStart, chromEnd, name, length, cpgNum, gcNum, perCpg, perGc, obsExp
sub load_cpg_data {
    my ($cpg_file) = @_;
    my $fh;
    if ($cpg_file =~ /\.gz$/) {
        $fh = IO::Uncompress::Gunzip->new($cpg_file)
          or die "Cannot open gzipped CpG file ($cpg_file): $GunzipError\n";
    } else {
        open($fh, "<", $cpg_file) or die "Cannot open CpG file ($cpg_file): $!\n";
    }
    # Process every line (no header)
    while (<$fh>) {
        chomp;
        next if $_ =~ /^\s*$/;  # skip blank lines
        my @fields = split "\t", $_;
        next unless @fields >= 11;
        my $chr   = $fields[1];         # column 2: chrom
        $chr = ($chr =~ /^chr/) ? $chr : "chr$chr";
        my $start = $fields[2];         # column 3: chromStart
        my $end   = $fields[3];         # column 4: chromEnd
        # Use perGc (column 10, index 9) as GC content and obsExp (column 11, index 10) as CpG ratio
        my $gc_content = $fields[9];
        my $obs_exp    = $fields[10];
        push @{$cpg{$chr}}, { start => $start, end => $end, gc_content => $gc_content, obs_exp => $obs_exp };
    }
    close $fh;
    
    # Sort each chromosome’s CpG islands by start coordinate
    foreach my $chr (keys %cpg) {
        @{$sorted_cpg{$chr}} = sort { $a->{start} <=> $b->{start} } @{$cpg{$chr}};
    }
}

### Load RepeatMasker Data (expected columns: chrom, start, end, repeat name, repeat class, strand)
sub load_repeatmasker_data {
    my ($repeat_file) = @_;
    my $fh;
    if ($repeat_file =~ /\.gz$/) {
        $fh = IO::Uncompress::Gunzip->new($repeat_file)
          or die "Cannot open gzipped RepeatMasker file ($repeat_file): $GunzipError\n";
    } else {
        open($fh, "<", $repeat_file) or die "Cannot open RepeatMasker file ($repeat_file): $!\n";
    }
    while (<$fh>) {
        chomp;
        my @fields = split "\t", $_;
        next unless @fields >= 6;
        my ($chr, $start, $end, $repeat_name, $repeat_class, $rstrand) = @fields[0,1,2,3,4,5];
        $chr = ($chr =~ /^chr/) ? $chr : "chr$chr";
        push @{$repeats{$chr}}, { start => $start, end => $end, name => $repeat_name, class => $repeat_class, strand => $rstrand };
    }
    close $fh;
    
    # Sort repeats for each chromosome by start coordinate
    foreach my $chr (keys %repeats) {
        @{$sorted_repeats{$chr}} = sort { $a->{start} <=> $b->{start} } @{$repeats{$chr}};
    }
}

############################################
# Annotation Retrieval Subroutines
############################################

### TSS Annotation: Return the closest TSS (gene ID, distance, gene type, gene synonyms)
sub check_proximity_to_tss {
    my ($chr, $pos) = @_;
    return ("None", "N/A", "N/A", "N/A") unless exists $tss{$chr};

    my $closest_gene = "None";
    my $min_distance = 1e9;
    my ($gene_type, $gene_synonyms) = ("N/A", "N/A");

    foreach my $tss_pos (keys %{$tss{$chr}}) {
        my $distance = abs($tss_pos - $pos);
        if ($distance < $min_distance) {
            $min_distance = $distance;
            $closest_gene = $tss{$chr}{$tss_pos}->{gene_id};
            $gene_type  = $tss{$chr}{$tss_pos}->{gene_type};
            $gene_synonyms = $tss{$chr}{$tss_pos}->{gene_synonyms};
        }
    }
    return ($closest_gene, $min_distance, $gene_type, $gene_synonyms);
}

### GFF Gene Annotation: Return the overlapping gene annotation (gene ID, type, synonyms)
# Only consider overlaps if the matched region (ms_start to ms_end) overlaps a gene
# and the gene’s strand matches the matched sequence strand.
sub get_gff_gene_annotation {
    my ($chr, $ms_start, $ms_end, $ms_strand) = @_;
    # Map matched sequence strand: F -> "+", R -> "-"
    my $mapped_strand = ($ms_strand eq "F") ? "+" : ($ms_strand eq "R" ? "-" : $ms_strand);
    my $best_gene;
    my $best_overlap = 0;
    if (exists $genes{$chr}) {
        foreach my $gene_start (keys %{$genes{$chr}}) {
            my $gene = $genes{$chr}{$gene_start};
            # Check if the matched sequence overlaps the gene region
            if ($ms_end >= $gene_start && $ms_start <= $gene->{end}) {
                # And check strand match
                if ($gene->{strand} eq $mapped_strand) {
                    # Calculate overlap length
                    my $overlap = min($ms_end, $gene->{end}) - max($ms_start, $gene_start) + 1;
                    if ($overlap > $best_overlap) {
                        $best_overlap = $overlap;
                        $best_gene = $gene;
                    }
                }
            }
        }
    }
    if ($best_gene) {
        return ($best_gene->{gene_id}, $best_gene->{gene_type}, $best_gene->{gene_synonyms});
    }
    else {
        return ("None", "N/A", "N/A");
    }
}

### CpG Island Annotation: Return whether the matched region overlaps a CpG island,
### and if so, its GC content and observed/expected ratio.
sub is_in_cpg_island {
    my ($chr, $ms_start, $ms_end) = @_;
    return ("No", 0, 0) unless exists $sorted_cpg{$chr};
    foreach my $region (@{$sorted_cpg{$chr}}) {
        # Check if the matched region overlaps the CpG island region
        if ($ms_end >= $region->{start} && $ms_start <= $region->{end}) {
            return ("Yes", $region->{gc_content}, $region->{obs_exp});
        }
        last if $region->{start} > $ms_end;
    }
    return ("No", 0, 0);
}

### Repeat Region Annotation: Return the repeat annotation if the matched region overlaps a repeat,
### and only if the repeat’s strand (from the file) matches the matched sequence’s strand.
sub is_in_repeat_region {
    my ($chr, $ms_start, $ms_end, $ms_strand) = @_;
    # Map matched sequence strand
    my $mapped_strand = ($ms_strand eq "F") ? "+" : ($ms_strand eq "R" ? "-" : $ms_strand);
    return ("No", "None", 0) unless exists $sorted_repeats{$chr};
    foreach my $region (@{$sorted_repeats{$chr}}) {
        # Check for overlap
        if ($ms_end >= $region->{start} && $ms_start <= $region->{end}) {
            # Check strand match (if available)
            if (exists $region->{strand} and $region->{strand} ne $mapped_strand) {
                next;
            }
            my $length = $region->{end} - $region->{start};
            return ($region->{name}, $region->{class}, $length);
        }
        last if $region->{start} > $ms_end;
    }
    return ("No", "None", 0);
}

############################################
# Processing Function for Matched Sequences
############################################

### Process matched sequences in parallel.
### Each line in the matched sequences file is assumed to have at least 9 tab‐separated fields:
### start, end, strand, matched sequence, marker_id, length, occurrences, upstream, downstream.
### The chromosome is assumed (or provided elsewhere); here we assume "chr1" if not available.
sub process_matched_sequences {
    my ($input_file) = @_;
    open my $fh, "<", $input_file or die "Cannot open input file ($input_file): $!\n";

    my $output_file = "matchedseqs_annotate.txt";
    open my $out_fh, ">", $output_file or die "Cannot open output file ($output_file): $!\n";
    # Output header with separate columns for each annotation source
    print $out_fh join("\t",
        "Start", "End", "Strand", "MatchedSeq", "Occurrences", "Chromosome",
        "TSS_Gene", "TSS_Distance", "TSS_Gene_Type", "TSS_Gene_Synonyms",
        "GFF_Gene", "GFF_Gene_Type", "GFF_Gene_Synonyms",
        "InCpG", "CpG_GC_Content", "CpG_ObsExp",
        "Repeat_Name", "Repeat_Class", "Repeat_Length"
    ), "\n";
    close $out_fh;

    my $temp_dir = "/tmp/perl_parallel";
    mkdir $temp_dir unless -d $temp_dir;

    while (<$fh>) {
        chomp;
        my @fields = split "\t", $_;
        next unless scalar @fields >= 9;
        # Extract matched sequence details:
        # Expected: start, end, strand, matched seq, marker_id, length, occurrences, upstream, downstream
        my ($start, $end, $strand, $seq, $marker_id, $length, $occurrences, $upstream, $downstream) = @fields;
        next unless ($start =~ /^\d+$/ && $end =~ /^\d+$/);
        # Assume chromosome is provided or use default ("chr1"); adjust if needed.
        my $chr = "chr1";

        $pm->start and next;

        # Retrieve annotations using the full region and strand:
        my ($tss_gene, $tss_distance, $tss_type, $tss_synonyms) = check_proximity_to_tss($chr, $start);
        my ($gff_gene, $gff_type, $gff_synonyms) = get_gff_gene_annotation($chr, $start, $end, $strand);
        my ($in_cpg, $gc_content, $cpg_ratio) = is_in_cpg_island($chr, $start, $end);
        my ($repeat_name, $repeat_class, $repeat_length) = is_in_repeat_region($chr, $start, $end, $strand);

        my $temp_file = "$temp_dir/process_$$.txt";
        open my $temp_fh, ">>", $temp_file or die "Cannot open temp file ($temp_file): $!\n";
        print $temp_fh join("\t",
            $start, $end, $strand, $seq, $occurrences, $chr,
            $tss_gene, $tss_distance, $tss_type, $tss_synonyms,
            $gff_gene, $gff_type, $gff_synonyms,
            $in_cpg, $gc_content, $cpg_ratio,
            $repeat_name, $repeat_class, $repeat_length
        ), "\n";
        close $temp_fh;

        $pm->finish;
    }
    close $fh;
    $pm->wait_all_children;

    # Merge temporary files into the final output file
    opendir my $dir, $temp_dir or die "Cannot open temp dir ($temp_dir): $!\n";
    open $out_fh, ">>", $output_file or die "Cannot open output file ($output_file): $!\n";
    while (my $file = readdir($dir)) {
        next unless $file =~ /^process_\d+\.txt$/;
        open my $in_fh, "<", "$temp_dir/$file" or next;
        while (<$in_fh>) { print $out_fh $_; }
        close $in_fh;
        unlink "$temp_dir/$file";
    }
    closedir $dir;
    close $out_fh;
}

############################################
# Main Execution
############################################

my ($input_file, $gff_file, $tss_file, $cpg_file, $repeat_file) = @ARGV;
die "Usage: $0 <matchedseqs.txt> <gff3[.gz]> <tss[.gz]> <cpg[.gz]> <repeatmasker[.gz]>\n"
    unless @ARGV == 5;

load_gene_annotation($gff_file);
load_tss_data($tss_file);
load_cpg_data($cpg_file);
load_repeatmasker_data($repeat_file);
process_matched_sequences($input_file);

print "Done.\n";
exit 0;