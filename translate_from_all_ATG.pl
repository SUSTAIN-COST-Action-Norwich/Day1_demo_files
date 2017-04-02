#!/usr/bin/perl -w
use strict;
use warnings;

#--------------------------------------------------------------------------
# translate_from_all_ATG.pl
#--------------------------------------------------------------------------
#
#
# This script reads in a nucleotide sequences in FASTA format from 
# "input_fasta_file" and translates from every ATG it finds on a forward and reverse
# directions for each sequence with at least 210 nucleotides or 70 amino acids.
# If this value need to be changed, modify $minimum_length to a new value.
# Then it reverse-complement the nucleotide sequences and repeats the 
# translation process.
#
# As a consequence, multiple peptides with overlapping sequences are created.
# 
# The script compensates for the 3rd codon "N" where possible.
#
#
# Usage: translate_from_all_ATG.pl <input_fasta_file> <output_file>
#
# @joewinnz
# 2015
#
#--------------------------------------------------------------------------


# Change the following value if different aa length cutoff is needed 
# for translation
my $minimum_length = 50;

#
# hash for translating DNA
#
my %DNAtoAA = ('GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TGT' => 'C',
			'TGC' => 'C', 'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
			'TTT' => 'F', 'TTC' => 'F', 'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G',
			'GGG' => 'G', 'CAT' => 'H', 'CAC' => 'H', 'ATT' => 'I', 'ATC' => 'I',
			'ATA' => 'I', 'AAA' => 'K', 'AAG' => 'K', 'TTG' => 'L', 'TTA' => 'L',
			'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 'ATG' => 'M',
			'AAT' => 'N', 'AAC' => 'N', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
			'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'CGT' => 'R', 'CGC' => 'R',
			'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 'TCT' => 'S',
			'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
			'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T', 'GTT' => 'V',
			'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'TGG' => 'W', 'TAT' => 'Y',
			'TAC' => 'Y', 'TAA' => '*', 'TAG' => '*', 'TGA' => '*',
			'ACN' => 'T', 'CCN' => 'P', 'CGN' => 'R', 'CTN' => 'L',
			'GCN' => 'A', 'GGN' => 'G', 'GTN' => 'V', 'TCN' => 'S');


#--------------------------------------------------------------------------
#
# Main
#
#--------------------------------------------------------------------------

my $usage = 'Usage: translate_from_all_ATG.pl <input_fasta_file> <output_file>';
unless ($ARGV[0] && $ARGV[1]) {die "Aborted\. Requires input and output files!!!\n$usage\n"};
my $input_file = $ARGV[0];
my $output_file = $ARGV[1];
my %nucl_sequences = ();

Read_fasta_file (\%nucl_sequences, $input_file);
open (OUTPUT, ">$output_file") || die "Can't create $output_file: $!\n";
Output_translations (\%nucl_sequences, $minimum_length, "forward", *OUTPUT);
Reverse_complement (\%nucl_sequences);
Output_translations (\%nucl_sequences, $minimum_length, "reverse", *OUTPUT);
close OUTPUT;

exit;


#--------------------------------------------------------------------------
# Subroutines
#--------------------------------------------------------------------------

sub Read_fasta_file {

#
#	Reads in a fasta input file containing amino acid or nucleotide sequences
#	and store them in a hash passed as a reference
#

	my ($sequences_ref, $file_name) = @_;
	my $total_seq_count = 0;
	my $header = "";
	my $seq = "";
	my $unique_ID = "";
	my $first = "YES";
	open (FASTA, "$file_name") || die "Can't open $file_name: $!\n";
	while (<FASTA>) {
		chomp;
		if (/^>/) {
			if ($first eq "YES") {
				$header = $_;
				$first = "NO";
			} elsif ($first eq "NO") {
				($unique_ID) = $header =~ /^>(\S+)\s?/;
				$seq = $header."\n".$seq."\n";
			if (exists $$sequences_ref{$unique_ID}) {
				print "$unique_ID is redundant and thus skipped\n";
			} else {
				$$sequences_ref{$unique_ID} = $seq;
				$total_seq_count++;
			}
			$seq = "";
			$unique_ID = "";
			$header = $_;
			} 
		} else {
			$seq .= $_;
		}
	}	
	#capture the last sequence entry...
	($unique_ID) = $header =~ /^>(\S+)\s?/;
	$$sequences_ref{$unique_ID} = $header."\n".$seq."\n";
	$total_seq_count++;
	close FASTA;
	print "There are $total_seq_count sequences in \"$file_name\" file\n";
}


#-----------------------------------------------------------------------------------


sub Output_translations {

#
#	Reads in nucleotide sequences from a hash passed as a reference and
#	print out the amino acid translations from every ATG it finds until
#	a stop codon or end of the nucleotide sequence as long as the translated
#	peptide is equal or greater than $min_length specified by the caller.
#

	my ($sequences_ref, $min_length, $orientation, $fh_output) = @_;
	if ($orientation eq "forward") { $orientation = "F" }
	elsif ($orientation eq "reverse") {$orientation = "R"}
	else {die "Can\'t understand the orientation"}
	my ($header, $seq, $this_seq, $ID, $description);
	foreach my $key(sort keys(%$sequences_ref)) {
	    $this_seq = $$sequences_ref{$key};
		if ($this_seq =~ /(>.*?)\n(.*)\n/) {
			$header = $1;
			$seq = $2;
		}
		#if ($orientation eq "R") {
		#	print "\nTranslating in reverse\.\.\. $header\n";
		#} else {
		#	print "\nTranslating\.\.\. $header\n";
		#}
		$seq =~ tr/[actgn]/[ACTGN]/;
		my $input_seq = $seq;
		my $input_seq_len = length($seq);
		my $peptide_num = 0;
		while ($seq =~ /(ATG\w+?$)/) {
			my $this_part = $1;
			my ($start, $end) = (0,0);
			$start = $input_seq_len-length($this_part) + 1;
			my $seq_length = length ($this_part);
			if ($seq_length < $min_length * 3) {		
				last;
			}
			my $peptide = Translate($this_part);
			my $pep_len = length($peptide);
			if ($peptide =~ /\*/) {
				$pep_len--;
			}
			if ($pep_len >= $min_length) {
				$peptide_num++;
				if ($header =~ /^(>\S+)( .*)$/) {
					$ID = $1;
					$description = $2;
				} elsif ($header =~ /^(>\S+$)/) {
					$ID = $1;
					$description = "";
				} else {
					die "The sequence is not in fasta format\n";
				}
				$end = $start + $pep_len * 3 - 1;
				my ($end_R, $start_R) = (0,0);
				if ($orientation eq "R"){
					$end_R = $input_seq_len - $start + 1;
					$start_R = $input_seq_len - $end + 1;
					# print ("$end_R \- $start_R\n");
					my $temp_seq = reverse $input_seq;
					$temp_seq =~ tr/ATGC/TACG/;
					# print substr($temp_seq, $start_R-1, $pep_len * 3), "\n";
					print $fh_output ($ID, "_pep_", $orientation, $peptide_num, " \[", $end_R, " - ", $start_R, "\] ", $description, "\n");
					print $fh_output ("$peptide\n");
				} else {
					# print ("$start \- $end\n");
					# print substr($input_seq, $start-1, $pep_len * 3), "\n";
					print $fh_output ($ID, "_pep_", $orientation, $peptide_num, " \[", $start, " - ", $end, "\] ", $description, "\n");
					print $fh_output ("$peptide\n");
				}
				$seq = substr($this_part, 3);
			} else {
				$seq = substr($this_part, 3);
			}
		}
	}
	return;		
}

#--------------------------------------------------------------------------

sub Translate {

#
#	Translates the nucleotide sequence passed in the argument by looking up 
#	the globally defined %DNAtoAA hash table for the genetic code and
#	returns the peptide sequence
#

	my $this_seq = shift;
	$this_seq =~ tr/[a-z]/[A-Z]/;
	my $pept;
	for (my $y = 0; $y <= (length($this_seq) - 3); $y += 3) {
		if (!defined $DNAtoAA{substr($this_seq, $y, 3)}) {
			#$pept .= "X";
			last;
		} else {
			$pept .= $DNAtoAA{substr($this_seq, $y, 3)};
			if ($DNAtoAA{substr($this_seq, $y, 3)} eq '*') {last}
		}
	}
	return $pept;
}

#--------------------------------------------------------------------------

sub Reverse_complement {

#
#	Reads in the nucleotide sequences from a hash passed as a reference and
#	replace the sequences with their reverse-complement
#

	my $sequences_ref = shift;
	foreach my $key(keys(%$sequences_ref)) {
		my $this_seq = $$sequences_ref{$key};
		$this_seq =~ /(>.*?)\n(.*)\n/;
		my $header = $1;
		my $seq = $2;
		$seq = reverse $seq;
		$seq =~ tr/ATGC/TACG/;
		$this_seq = $header."\n".$seq."\n";
		$$sequences_ref{$key} = $this_seq;
	}
}

#--------------------------------------------------------------------------
