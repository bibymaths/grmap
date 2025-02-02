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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
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

# Hashes & arrays for genomic features
my (%genes, %tss, %cpg, %repeats, %sorted_cpg, %sorted_repeats);

# Detect number of CPU cores
my $num_cores = `lscpu -p | grep -v '^#' | wc -l`;
chomp($num_cores);
my $pm = Parallel::ForkManager->new($num_cores);

# Load GFF3 Gene Annotations
sub load_gene_annotation {
    my ($gff_file) = @_;
    open my $fh, "-|", "zcat $gff_file" or die "Cannot open GFF3 file: $!\n";

    while (<$fh>) {
        next if /^#/;
        my @fields = split "\t";
        next unless @fields >= 9 && $fields[2] eq "gene";

        my ($chr, $start, $end, $strand) = @fields[0, 3, 4, 6];
        my ($gene_id) = $fields[8] =~ /ID=([^;]+)/;
        my ($gene_type) = $fields[8] =~ /biotype=([^;]+)/;  # Extract gene type
        my ($gene_synonyms) = $fields[8] =~ /Name=([^;]+)/; # Extract gene name/synonyms

        $genes{$chr}{$start} = { 
            end => $end, strand => $strand, gene_id => $gene_id,
            gene_type => ($gene_type // "Unknown"), gene_synonyms => ($gene_synonyms // "N/A")
        };
    }
    close $fh;
}

# Load TSS Data
sub load_tss_data {
    my ($tss_file) = @_;
    open my $fh, "-|", "zcat $tss_file" or die "Cannot open TSS file: $!\n";
    <$fh>;  # Skip header

    while (<$fh>) {
        my @cols = split "\t";
        next unless @cols >= 5;

        my ($gene_id, $transcript_id, $start, $end, $chr) = @cols[0, 1, 2, 3, 4];
        $chr = "chr$chr";
        $tss{$chr}{$start} = { gene_id => $gene_id, transcript_id => $transcript_id, start => $start, end => $end };
    }
    close $fh;
}

# Load CpG Island Data
sub load_cpg_data {
    my ($cpg_file) = @_;
    open my $fh, "-|", "zcat $cpg_file" or die "Cannot open CpG file: $!\n";

    while (<$fh>) {
        my ($chr, $start, $end, $gc_content, $obs_exp_cpg) = split "\t"; # Add GC content & CpG ratio if available
        push @{$cpg{$chr}}, { start => $start, end => $end, gc_content => $gc_content, obs_exp => $obs_exp_cpg };
    }
    close $fh;

    # Sort CpG islands for fast binary search
    foreach my $chr (keys %cpg) {
        @{$sorted_cpg{$chr}} = sort { $a->{start} =~ /(\d+)/ <=> $b->{start} =~ /(\d+)/ } @{$cpg{$chr}};
    }
}

# Load RepeatMasker Data
sub load_repeatmasker_data {
    my ($repeat_file) = @_;
    open my $fh, "-|", "zcat $repeat_file" or die "Cannot open RepeatMasker file: $!\n";

    while (<$fh>) {
        my ($chr, $start, $end, $repeat_name, $repeat_class) = split "\t";
        push @{$repeats{$chr}}, { start => $start, end => $end, name => $repeat_name, class => $repeat_class };
    }
    close $fh;

    foreach my $chr (keys %repeats) {
        @{$sorted_repeats{$chr}} = sort { $a->{start} <=> $b->{start} } @{$repeats{$chr}};
    }
}

# Find closest TSS
sub check_proximity_to_tss {
    my ($chr, $pos) = @_;
    return ("None", "N/A") unless exists $tss{$chr};

    my $closest_gene = "None";
    my $min_distance = 1e9;

    foreach my $tss_pos (keys %{$tss{$chr}}) {
        my $distance = abs($tss_pos - $pos);
        if ($distance < $min_distance) {
            $min_distance = $distance;
            $closest_gene = $tss{$chr}{$tss_pos}->{gene_id};
        }
    }
    return ($closest_gene, $min_distance);
}

# Find CpG properties
sub is_in_cpg_island {
    my ($chr, $pos) = @_;
    return ("No", 0, 0) unless exists $sorted_cpg{$chr};
    foreach my $region (@{$sorted_cpg{$chr}}) {
        return ("Yes", $region->{gc_content}, $region->{obs_exp}) if $pos >= $region->{start} && $pos <= $region->{end};
        last if $region->{start} > $pos;
    }
    return ("No", 0, 0);
}

# Find Repeat Region properties
sub is_in_repeat_region {
    my ($chr, $pos) = @_;
    return ("No", "None", 0) unless exists $sorted_repeats{$chr};
    foreach my $region (@{$sorted_repeats{$chr}}) {
        return ($region->{name}, $region->{class}, $region->{end} - $region->{start}) if $pos >= $region->{start} && $pos <= $region->{end};
        last if $region->{start} > $pos;
    }
    return ("No", "None", 0);
}

# Process matched sequences in parallel
sub process_matched_sequences {
    my ($input_file) = @_;
    open my $fh, "<", $input_file or die "Cannot open input file: $!\n";

    my $output_file = "matchedseqs_annotate.txt";
    open my $out_fh, ">", $output_file or die "Cannot open output file: $!\n";
    print $out_fh join("\t", "Start", "End", "Strand", "MatchedSeq", "Occurrences",
                      "Chromosome", "ClosestGene", "DistToTSS", "GeneType", "GeneSynonyms",
                      "InCpG", "GCContent", "CpG_Ratio", "RepeatName", "RepeatClass", "RepeatLength"), "\n";
    close $out_fh;

    my $temp_dir = "/tmp/perl_parallel";
    mkdir $temp_dir unless -d $temp_dir;

    while (<$fh>) {
        chomp;
        my @fields = split /\t/;
        next unless scalar @fields >= 9;

        my ($start, $end, $strand, $seq, $marker_id, $length, $occurrences, $upstream, $downstream) = @fields;
        next unless ($start =~ /^\d+$/ && $end =~ /^\d+$/);
        my $chr = "chr1";

        $pm->start and next;  

        my ($closest_gene, $tss_distance) = check_proximity_to_tss($chr, $start);
        my ($in_cpg, $gc_content, $cpg_ratio) = is_in_cpg_island($chr, $start);
        my ($repeat_name, $repeat_class, $repeat_length) = is_in_repeat_region($chr, $start);

        my $temp_file = "$temp_dir/process_$$.txt";
        open my $temp_fh, ">>", $temp_file or die "Cannot open temp file: $!\n";
        print $temp_fh join("\t", 
            $start, $end, $strand, $seq, $occurrences, 
            $chr, 
            $closest_gene // "None",              
            $tss_distance // "N/A",               
            $genes{$chr}{$start}->{gene_type} // "N/A",      
            $genes{$chr}{$start}->{gene_synonyms} // "N/A",  
            $in_cpg // "No",                       
            $gc_content // "0",                    
            $cpg_ratio // "0",                     
            $repeat_name // "No",                  
            $repeat_class // "None",               
            $repeat_length // "0"                   
        ), "\n";
        close $temp_fh;

        $pm->finish;
    }
    close $fh;
    $pm->wait_all_children;  

    opendir my $dir, $temp_dir or die "Cannot open temp dir: $!\n";
    open $out_fh, ">>", $output_file or die "Cannot open output file: $!\n";
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
# Load data
my ($input_file, $gff_file, $tss_file, $cpg_file, $repeat_file) = @ARGV;
die "Usage: $0 <matchedseqs.txt> <gff3.gz> <tss.gz> <cpg.gz> <repeatmasker.gz>\n" unless @ARGV == 5;
load_gene_annotation($gff_file);
load_tss_data($tss_file);
load_cpg_data($cpg_file);
load_repeatmasker_data($repeat_file);
process_matched_sequences($input_file);   
print "Done.\n";
exit 0; 