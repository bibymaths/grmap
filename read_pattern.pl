#!/usr/bin/env perl
use strict;
use warnings; 
use diagnostics;
use bytes;  # Force byte semantics for DNA sequences
use Time::HiRes qw(gettimeofday tv_interval);
 
#-----------------------------------------------------
# Usage message and command-line argument processing
#-----------------------------------------------------
sub usage {
    print <<"END_USAGE";
Usage: $0 <reads_choice> <query_choice>

  <reads_choice> options:
      1   Use illumina_reads_40.fasta.gz
      2   Use illumina_reads_60.fasta.gz
      3   Use illumina_reads_80.fasta.gz
      4   Use illumina_reads_100.fasta.gz

  <query_choice> options:
      1   Load 1,000 queries
      2   Load 10,000 queries
      3   Load 100,000 queries
      4   Load 1,000,000 queries

Example:
  $0 2 3
    => Use illumina_reads_60.fasta.gz and load 100,000 queries.

END_USAGE
    exit 1;
}

# If no arguments or help flag provided, show usage.
if ( @ARGV != 2 or $ARGV[0] eq '--help' or $ARGV[0] eq '-h' ) {
    usage();
}

my ($reads_choice, $query_choice) = @ARGV; 

#-----------------------------------------------------
# Record initial time and memory usage
#-----------------------------------------------------
my $start_time      = [gettimeofday()];
my $start_timestamp = scalar(localtime());
my $start_mem       = get_memory_usage();

#-----------------------------------------------------
# File and Parameter Setup
#-----------------------------------------------------
my $ref_file   = 'hg38_partial.fasta.gz';
my %reads_map  = (
    1 => 'illumina_reads_40.fasta.gz',
    2 => 'illumina_reads_60.fasta.gz',
    3 => 'illumina_reads_80.fasta.gz',
    4 => 'illumina_reads_100.fasta.gz',
);
my @query_opts = ( 1000, 10_000, 100_000, 1_000_000 );
my $result_file = 'MatchedMarkers.txt'; 
 
# Validate command-line arguments.
exists $reads_map{$reads_choice} or die "Invalid read file choice: $reads_choice\n";
$query_choice =~ /^[1-4]$/ or die "Invalid query choice: $query_choice\n";
my $num_queries = $query_opts[$query_choice - 1];

#-----------------------------------------------------
# Read Reference Genome
#-----------------------------------------------------
my $DNA = read_fasta($ref_file);
die "Reference genome is empty!\n" unless length($DNA);

#-----------------------------------------------------
# Load marker sequences (queries)
# Now each query is stored as a hash with key 'id' and 'seq'
#-----------------------------------------------------
my $reads_file = $reads_map{$reads_choice};
my @queries = read_n_sequences($reads_file, $num_queries);
my $actual_queries = scalar @queries;

#-----------------------------------------------------
# Detect number of cores (using external command)
#-----------------------------------------------------
my $num_cores = `lscpu -p | grep -v '^#' | wc -l`;
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
# Combine temporary files into final output file with header
#-----------------------------------------------------
open(my $OUT, '>', $result_file) or die "Cannot open $result_file: $!";
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
print $OUT "# Columns: Start\tEnd\tStrand\tMatchedSequence\tMarkerID\tLength\tOccurrenceCount\tUpstreamContext\tDownstreamContext\n";

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
close $OUT or die "Cannot close $result_file: $!";
print "\nDone. $total_matched match(es) reported.\n";
print "See '$result_file' for the matched markers.\n";
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
                last if @list == $max;
                $seq = '';
            }
            $id = $1;  # capture header (without '>')
        }
        else {
            $seq .= $line;
        }
    }
    if ($seq ne '' and @list < $max) {
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