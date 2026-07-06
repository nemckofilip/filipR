#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# -----------------------------------------------------------
# Script: extract_umi_read1.pl
# Purpose:
#   Extract first N bases from Read 1 (UMI), append to header
#   of both R1 & R2, and trim R1 accordingly.
#
# Example:
#   perl extract_umi_read1.pl --r1 r1.fq.gz --r2 r2.fq.gz \
#       --out1 r1_umi.fq.gz --out2 r2_umi.fq.gz \
#       --umi-len 8 --threads 4
# -----------------------------------------------------------

my ($r1_in, $r2_in, $r1_out, $r2_out, $umi_len, $threads, $help);
$umi_len = 8;
$threads = $ENV{SLURM_CPUS_PER_TASK} // 4;

GetOptions(
    "r1=s"       => \$r1_in,
    "r2=s"       => \$r2_in,
    "out1=s"     => \$r1_out,
    "out2=s"     => \$r2_out,
    "umi-len=i"  => \$umi_len,
    "threads=i"  => \$threads,
    "help|h"     => \$help,
) or die "Error parsing options. Use --help.\n";

if ($help || !$r1_in || !$r1_out) {
    print <<"USAGE";
Usage:
  # Paired-end
  perl extract_umi_read1.pl --r1 <R1.fq.gz> --r2 <R2.fq.gz> \\
       --out1 <out1.fq.gz> --out2 <out2.fq.gz> \\
       [--umi-len 8] [--threads 4]

  # Single-end (R1 only)
  perl extract_umi_read1.pl --r1 <R1.fq.gz> --out1 <out1.fq.gz> \\
       [--umi-len 8] [--threads 4]

Description:
  Extracts first N bases from R1 (UMI), appends them to the R1 (and, in
  paired-end mode, R2) read headers, trims R1 by N bases, and outputs
  gzipped FASTQ.

Options:
  --r1, --r2       Input FASTQ files (gzipped). --r2 is optional (single-end if omitted).
  --out1, --out2   Output FASTQ files (gzipped). --out2 required only in paired-end mode.
  --umi-len        Number of bases to extract from R1 (default 8)
  --threads        Threads for pigz/gzip compression (default SLURM value or 4)
  --help           Show this help message
USAGE
    exit 1;
}

# Paired-end if R2 provided; require both --r2 and --out2 together
my $paired = (defined $r2_in || defined $r2_out) ? 1 : 0;
if ($paired && (!defined $r2_in || !defined $r2_out)) {
    die "For paired-end mode, provide both --r2 and --out2.\n";
}

# -----------------------------------------------------------
# Compression binary detection
# -----------------------------------------------------------
my $have_pigz = system("command -v pigz >/dev/null 2>&1") == 0 ? 1 : 0;
my $zcat_cmd  = $have_pigz ? "pigz -p $threads -dc"   : "zcat";
my $gzip_cmd  = $have_pigz ? "pigz -p $threads -c"    : "gzip -c";
warn "Note: pigz not found; using standard gzip (slower).\n"
    unless $have_pigz;

# -----------------------------------------------------------
# Open input files
# -----------------------------------------------------------
open(my $IN1, "$zcat_cmd '$r1_in' |") or die "Cannot open R1: $r1_in\n";
open(my $OUT1, "| $gzip_cmd > '$r1_out'") or die "Cannot write R1 output: $r1_out\n";

my ($IN2, $OUT2);
if ($paired) {
    open($IN2, "$zcat_cmd '$r2_in' |") or die "Cannot open R2: $r2_in\n";
    open($OUT2, "| $gzip_cmd > '$r2_out'") or die "Cannot write R2 output: $r2_out\n";
}

# -----------------------------------------------------------
# Processing loop
# -----------------------------------------------------------
my $count = 0;

while (my $h1 = <$IN1>) {
    my $s1 = <$IN1>;
    my $p1 = <$IN1>;
    my $q1 = <$IN1>;

    my ($h2, $s2, $p2, $q2);
    if ($paired) {
        $h2 = <$IN2>;
        $s2 = <$IN2>;
        $p2 = <$IN2>;
        $q2 = <$IN2>;

        die "Error: R2 ended before R1 at read $count\n"
            unless defined $h2;
    }

    chomp($s1);
    chomp($q1);

    # Validate read length
    if (length($s1) < $umi_len) {
        die "Error: R1 is shorter than UMI length ($umi_len bp) at read $count\n";
    }

    # Extract UMI (first N bases of R1)
    my $umi = substr($s1, 0, $umi_len);

    # Trim R1
    my $s1_trim = substr($s1, $umi_len);
    my $q1_trim = substr($q1, $umi_len);

    # Fix header (R1) and write trimmed R1
    chomp($h1);
    my ($id1, $desc1) = split(/\s+/, $h1, 2);
    my $new_h1 = $id1 . "_" . $umi . (defined $desc1 ? " $desc1" : "");
    print $OUT1 "$new_h1\n$s1_trim\n$p1$q1_trim\n";

    # Write R2 unchanged, with UMI appended to header
    if ($paired) {
        chomp($h2);
        my ($id2, $desc2) = split(/\s+/, $h2, 2);
        my $new_h2 = $id2 . "_" . $umi . (defined $desc2 ? " $desc2" : "");
        print $OUT2 "$new_h2\n$s2$p2$q2";
    }

    $count++;
}

close $IN1;
close $OUT1;
if ($paired) {
    close $IN2;
    close $OUT2;
}

print "Done. Processed $count " . ($paired ? "read pairs" : "reads") . ".\n";