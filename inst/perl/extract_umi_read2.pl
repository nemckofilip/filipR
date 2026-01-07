#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# -----------------------------------------------------------
# Script: extract_umi.pl
# Purpose:
#   Extract first N bases from Read 2 (UMI), append to header
#   of both R1 & R2, and trim R2 accordingly.
#
# Example:
#   perl extract_umi.pl --r1 r1.fq.gz --r2 r2.fq.gz \
#       --out1 r1_umi.fq.gz --out2 r2_umi.fq.gz \
#       --umi-len 8 --threads 4
# -----------------------------------------------------------

my ($r1_in, $r2_in, $r1_out, $r2_out, $umi_len, $threads, $help);

$umi_len = 8;  # default UMI length
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

if ($help || !$r1_in || !$r2_in || !$r1_out || !$r2_out) {
    print <<"USAGE";
Usage:
  perl extract_umi.pl --r1 <R1.fq.gz> --r2 <R2.fq.gz> \\
       --out1 <out1.fq.gz> --out2 <out2.fq.gz> \\
       [--umi-len 8] [--threads 4]

Description:
  Extracts first N bases from R2 (UMI), appends them to both R1 & R2
  read headers, trims R2 by N bases, and outputs gzipped FASTQ.

Options:
  --r1, --r2       Input FASTQ files (gzipped)
  --out1, --out2   Output FASTQ files (gzipped)
  --umi-len        Number of bases to extract from R2 (default 8)
  --threads        Threads for pigz/gzip compression (default SLURM value or 4)
  --help           Show this help message
USAGE
    exit 1;
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
open(my $IN2, "$zcat_cmd '$r2_in' |") or die "Cannot open R2: $r2_in\n";

# Open output files
open(my $OUT1, "| $gzip_cmd > '$r1_out'") or die "Cannot write R1 output: $r1_out\n";
open(my $OUT2, "| $gzip_cmd > '$r2_out'") or die "Cannot write R2 output: $r2_out\n";

# -----------------------------------------------------------
# Processing loop
# -----------------------------------------------------------
my $count = 0;

while (my $h1 = <$IN1>) {
    my $s1 = <$IN1>;
    my $p1 = <$IN1>;
    my $q1 = <$IN1>;

    my $h2 = <$IN2>;
    my $s2 = <$IN2>;
    my $p2 = <$IN2>;
    my $q2 = <$IN2>;

    die "Error: R2 ended before R1 at read $count\n"
        unless defined $h2;

    chomp($s2);
    chomp($q2);

    # Validate read length
    if (length($s2) < $umi_len) {
        die "Error: R2 is shorter than UMI length ($umi_len bp) at read $count\n";
    }

    # Extract UMI (first N bases)
    my $umi = substr($s2, 0, $umi_len);

    # Trim R2
    my $s2_trim = substr($s2, $umi_len);
    my $q2_trim = substr($q2, $umi_len);

    # Fix headers
    chomp($h1);
    chomp($h2);

    my ($id1, $desc1) = split(/\s+/, $h1, 2);
    my ($id2, $desc2) = split(/\s+/, $h2, 2);

    my $new_h1 = $id1 . "_" . $umi . (defined $desc1 ? " $desc1" : "");
    my $new_h2 = $id2 . "_" . $umi . (defined $desc2 ? " $desc2" : "");

    # Write outputs
    print $OUT1 "$new_h1\n$s1$p1$q1";
    print $OUT2 "$new_h2\n$s2_trim\n$p2$q2_trim\n";

    $count++;
}

close $IN1;
close $IN2;
close $OUT1;
close $OUT2;

print "Done. Processed $count read pairs.\n";

