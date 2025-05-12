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
use bytes;
use Time::HiRes qw(gettimeofday tv_interval);
 
#-----------------------------------------------------
# Record initial time and memory usage
#-----------------------------------------------------
my $start_time      = [gettimeofday()];
my $start_timestamp = scalar(localtime());
my $start_mem       = get_memory_usage();

#----------------------------------------
# Usage: match.pl <reads.fasta.gz> <ref.fasta.gz> [<max_queries>]
#----------------------------------------
sub usage {
  die <<"USAGE";
Usage: $0 <reads_file> <reference_file> [<max_queries>]
   <reads_file>     FASTA or FASTA.GZ of your markers
   <reference_file> FASTA or FASTA.GZ of your reference
   <max_queries>    (optional) number of reads to load; default = all
USAGE
}

# require at least two arguments
usage() if @ARGV < 2;

# assign files directly
my ($reads_file, $ref_file, $num_queries) = @ARGV;
$num_queries //= 0;   # zero or undef -> load _all_ sequences
#my $result_file = 'matchedseqs.txt';
my $OUT = *STDOUT;    # send everything to STDOUT

#-----------------------------------------------------
# Read Reference Genome
#-----------------------------------------------------
my $DNA = read_fasta($ref_file);
die "Reference genome is empty!\n" unless length($DNA);

#-----------------------------------------------------
# Load marker sequences (queries)
# Now each query is stored as a hash with key 'id' and 'seq'
#-----------------------------------------------------
my @queries = $num_queries
  ? read_n_sequences($reads_file, $num_queries)
  : read_n_sequences($reads_file);   # load all if no limit
my $actual_queries = scalar @queries;

#-----------------------------------------------------
# Detect number of cores
#-----------------------------------------------------
my $num_cores = `lscpu -p | grep -v '^#' | wc -l`;
#my $num_cores = 2;
chomp($num_cores);

#-----------------------------------------------------
# Split queries among child processes (forking)
#-----------------------------------------------------
my $total      = scalar @queries;
my $base_chunk = int($total / $num_cores);
my $remainder  = $total % $num_cores;

my @child_pids;
my $start_index = 0;
for (my $i = 0; $i < $num_cores; $i++) {
    my $chunk_size = $base_chunk + ($i < $remainder ? 1 : 0);
    last if $chunk_size <= 0;
    my @chunk = @queries[$start_index .. $start_index + $chunk_size - 1];
    $start_index += $chunk_size;

    my $pid = fork();
    defined $pid or die "Cannot fork: $!";
    if ($pid == 0) {
        process_chunk(\@chunk, $DNA);
        exit(0);
    }
    else {
        push @child_pids, $pid;
    }
}

# Wait for all children to finish.
foreach my $pid (@child_pids) {
    waitpid($pid, 0);
}

#-----------------------------------------------------
# Record end time and memory usage
#-----------------------------------------------------
my $end_time      = [gettimeofday()];
my $end_timestamp = scalar(localtime());
my $elapsed       = tv_interval($start_time, $end_time);
my $end_mem       = get_memory_usage();

#-----------------------------------------------------
# Combine temporary files into final output
#-----------------------------------------------------

#open(my $OUT, '>', $result_file) or die "Cannot open $result_file: $!";

print $OUT "# ReferenceFile: $ref_file\n";
print $OUT "# QueryFile: $reads_file\n";
print $OUT "# Total Queries Processed: $actual_queries\n";
print $OUT "# Query Option (number of queries): $num_queries\n";
print $OUT "# Number of Cores Used: $num_cores\n";
print $OUT "# Run Start: $start_timestamp\n";
print $OUT "# Run End:   $end_timestamp\n";
print $OUT "# Elapsed Time (s): $elapsed\n";
print $OUT "# Memory Usage at Start: " . format_memory_usage($start_mem) . "\n";
print $OUT "# Memory Usage at End:   " . format_memory_usage($end_mem) . "\n";
print $OUT "Start\tEnd\tStrand\tMatchedSequence\tMarkerID\tLength\tOccurrenceCount\tUpstreamContext\tDownstreamContext\n";

my $total_matched = 0;
foreach my $child_pid (@child_pids) {
    my $tmp_file = "tmp_child_$child_pid.txt";
    if ( -e $tmp_file ) {
        open(my $fh, '<', $tmp_file) or die "Cannot open $tmp_file: $!";
        while (<$fh>) {
            print $OUT $_;
            $total_matched++;
        }
        close $fh;
        unlink $tmp_file or warn "Couldn't remove $tmp_file: $!";
    }
}
#-----------------------------------------------------
#close $OUT or die "Cannot close $result_file: $!";
#print "\nDone. $total_matched match(es) reported.\n";
#print "See '$result_file' for the matched markers.\n";
#-----------------------------------------------------
exit 0;

#=====================================================
# SUBROUTINES
#=====================================================

#-----------------------------------------------------
# Subroutine: get_memory_usage
# Returns the current memory usage in KB by calling ps.
# (Works on Linux systems.)
#-----------------------------------------------------
sub get_memory_usage {
    my $pid = $$;
    my $mem = `ps -o rss= -p $pid`;
    chomp($mem);
    return $mem;  # in KB
}

#-----------------------------------------------------
# Subroutine: format_memory_usage
# Converts a memory value in KB into a human-readable string,
# using KB, MB, or GB as appropriate.
#-----------------------------------------------------
sub format_memory_usage {
    my $kb = shift;
    if ($kb < 1024) {
        return "$kb KB";
    }
    my $mb = $kb / 1024;
    if ($mb < 1024) {
        return sprintf("%.2f MB", $mb);
    }
    my $gb = $mb / 1024;
    return sprintf("%.2f GB", $gb);
}

#-----------------------------------------------------
# Subroutine: read_fasta
# Reads a (possibly gzipped) FASTA file into a single scalar.
#-----------------------------------------------------
sub read_fasta {
    my ($file) = @_;
    my $seq = '';
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, '-|', "gunzip -c $file") or die "Can't open $file: $!";
    }
    else {
        open($fh, '<', $file) or die "Can't open $file: $!";
    }
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^>/;
        $seq .= $line;
    }
    close $fh;
    return $seq;
}

#-----------------------------------------------------
# Subroutine: read_n_sequences
# Reads up to $max sequences from a FASTA file.
# Returns an array of hash references with keys 'id' and 'seq'.
#-----------------------------------------------------
sub read_n_sequences {
    my ($file, $max) = @_;
    my @list;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, '-|', "gunzip -c $file") or die "Can't open $file: $!";
    }
    else {
        open($fh, '<', $file) or die "Can't open $file: $!";
    }
    my $seq = '';
    my $id  = '';
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>(.*)/) {
            if ($seq ne '') {
                push @list, { id => $id, seq => $seq };
                # only stop early if the caller really passed a max>0
                last if defined($max) && $max > 0 && @list == $max;
                $seq = '';
            }
            $id = $1;  # capture header (without '>')
        }
        else {
            $seq .= $line;
        }
    }
    # always push the final sequence
    if ($seq ne '') {
        push @list, { id => $id, seq => $seq };
    }
    close $fh;
    return @list;
}

#-----------------------------------------------------
# Subroutine: process_chunk
# Each child process handles its chunk of queries.
# For each marker, both forward and reverse-complement searches
# are performed using a combined loop over the search patterns.
# For each match, additional metadata is saved:
#   - Start and end coordinates in the reference.
#   - Strand (F for forward, R for reverse).
#   - The matched sequence.
#   - Marker ID (from the FASTA header).
#   - Marker length.
#   - Occurrence count (total matches for that marker).
#   - Upstream (20 bp) and downstream (20 bp) context.
#-----------------------------------------------------
sub process_chunk {
    my ($chunk_ref, $DNA_ref) = @_;
    my @output_lines;
    my $context    = 20;              # bases of flanking context
    my $ref_length = length($DNA_ref);  # cache reference length

    foreach my $marker_ref (@$chunk_ref) {
        my $marker    = $marker_ref->{seq};
        my $marker_id = $marker_ref->{id} // "Unknown";
        my $len       = length($marker);
        # Compute reverse complement once.
        my $rev = revcomp($marker);
        # Build an array of search patterns and associated strand.
        my @patterns = ( { pattern => $marker, strand => "F" } );
        push @patterns, { pattern => $rev, strand => "R" } if ($marker ne $rev);
        my @matches;
        foreach my $p (@patterns) {
            my $pattern = $p->{pattern};
            my $strand  = $p->{strand};
            my $pos = 0;
            while ( (my $found = index($DNA_ref, $pattern, $pos)) != -1 ) {
                my $end = $found + $len - 1;
                my $upstream   = $found >= $context ? substr($DNA_ref, $found - $context, $context) : substr($DNA_ref, 0, $found);
                my $downstream = ($ref_length - $end - 1) >= $context ? substr($DNA_ref, $end + 1, $context) : substr($DNA_ref, $end + 1);
                push @matches, {
                    start      => $found,
                    end        => $end,
                    strand     => $strand,
                    matched    => $pattern,
                    upstream   => $upstream,
                    downstream => $downstream
                };
                $pos = $found + 1;
            }
        }
        my $occurrence_count = scalar(@matches);
        foreach my $match (@matches) {
            my $line = join("\t",
                $match->{start},
                $match->{end},
                $match->{strand},
                $match->{matched},
                $marker_id,
                $len,
                $occurrence_count,
                $match->{upstream},
                $match->{downstream}
            );
            push @output_lines, $line;
        }
    }
    my $tmp_file = "tmp_child_" . $$ . ".txt";
    open(my $fh, '>', $tmp_file) or die "Cannot open temporary file $tmp_file: $!";
    print $fh "$_\n" for @output_lines;
    close $fh;
}

#-----------------------------------------------------
# Subroutine: revcomp
# Returns the reverse complement of a DNA sequence.
#-----------------------------------------------------
sub revcomp {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}