#!/usr/bin/env perl
use strict;
use warnings;

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

#-----------------------------------------------------
# Read Reference Genome
#-----------------------------------------------------
print "Loading reference genome...\n";
my $DNA = read_fasta($ref_file);
die "Reference genome is empty!\n" unless length($DNA);
print "Reference genome loaded (", length($DNA), " bp).\n";

#-----------------------------------------------------
# Ask user for which read file and number of queries
#-----------------------------------------------------
print "\nWhich read file do you want to use?\n",
      "1. illumina_reads_40.fasta.gz\n",
      "2. illumina_reads_60.fasta.gz\n",
      "3. illumina_reads_80.fasta.gz\n",
      "4. illumina_reads_100.fasta.gz\n",
      "\nPlease enter option number: ";
chomp(my $reads_choice = <STDIN>);
exists $reads_map{$reads_choice} or die "Invalid choice. Exiting.\n";

print "\nSelect number of queries:\n",
      "1. 1,000\n",
      "2. 10,000\n",
      "3. 100,000\n",
      "4. 1,000,000\n",
      "Please enter option number: ";
chomp(my $query_choice = <STDIN>);
$query_choice =~ /^[1-4]$/ or die "Invalid query choice. Exiting.\n";
my $num_queries = $query_opts[$query_choice - 1];

#-----------------------------------------------------
# Load marker sequences (queries)
#-----------------------------------------------------
my $reads_file = $reads_map{$reads_choice};
print "\nLoading up to $num_queries markers from $reads_file ...\n";
my @queries = read_n_sequences($reads_file, $num_queries);
my $actual_queries = scalar @queries;
print "Loaded $actual_queries marker sequences.\n";

#-----------------------------------------------------
# Detect number of cores (using external commands)
#-----------------------------------------------------
my $num_cores;
if ($^O eq 'linux') {
    $num_cores = `nproc`;
    chomp($num_cores);
} elsif ($^O eq 'darwin') {
    $num_cores = `sysctl -n hw.ncpu`;
    chomp($num_cores);
} else {
    $num_cores = 1;
}
$num_cores = 1 unless $num_cores =~ /^\d+$/;
print "\nDetected $num_cores core(s). Splitting work among them...\n";

#-----------------------------------------------------
# Split queries among child processes (forking)
#-----------------------------------------------------
my $total = scalar @queries;
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
# Combine temporary files into final output file
#-----------------------------------------------------
print "\nCombining results into '$result_file'...\n";
open(my $OUT, '>', $result_file) or die "Cannot open $result_file: $!";
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
print "\nDone. $total_matched marker(s) matched at least once.\n";
print "See '$result_file' for the matched markers.\n";
exit 0;

#-----------------------------------------------------
# Subroutine: process_chunk
# Each child process handles its chunk of queries.
#-----------------------------------------------------
sub process_chunk {
    my ($chunk_ref, $DNA_ref) = @_;
    my @matched;
    foreach my $marker (@$chunk_ref) {
        my $count = index_count($DNA_ref, $marker);
        if ($count > 0) {
            push @matched, $marker;
        }
    }
    my $tmp_file = "tmp_child_" . $$ . ".txt";
    open(my $fh, '>', $tmp_file) or die "Cannot open temporary file $tmp_file: $!";
    print $fh "$_\n" for @matched;
    close $fh;
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
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            if ($seq ne '') {
                push @list, $seq;
                last if @list == $max;
                $seq = '';
            }
            next;
        }
        else {
            $seq .= $line;
        }
    }
    push @list, $seq if $seq ne '' and @list < $max;
    close $fh;
    return @list;
}

#-----------------------------------------------------
# Subroutine: index_count
# Counts occurrences of $small in $big using index().
#-----------------------------------------------------
sub index_count {
    my ($big, $small) = @_;
    return 0 if length($small) > length($big);
    my $count = 0;
    my $pos   = 0;
    while (1) {
        $pos = index($big, $small, $pos);
        last if $pos == -1;
        $count++;
        $pos++;  # move past current match
    }
    return $count;
}