#!/usr/bin/perl -w
use strict;
use warnings;

#--------------------------------------------------------------------------
# print_all_orfs_from_ATG.pl
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
# Usage: print_all_orfs_from_ATG.pl <input_nucleotide_fasta_file> <output_file>
#
# Joe Win, The Sainsbury Laboratory, Norwich NR4 7UH, UK
# joe.win@tsl.ac.uk
# 2011©
#
#--------------------------------------------------------------------------


# Change the following value if different nt length cutoff is needed 
# for orf length
my $minimum_length = 210;

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

my $usage = 'Usage: print_all_orfs_from_ATG.pl <input_nucleotide_fasta_file> <output_file>';
unless ($ARGV[0] && $ARGV[1]) {die "Aborted\. Requires input and output files!!!\n$usage\n"};
my $input_file = $ARGV[0];
my $output_file = $ARGV[1];
my %nucl_sequences = ();

Read_fasta_file (\%nucl_sequences, $input_file);
open (OUTPUT, ">$output_file") || die "Can't create $output_file: $!\n";
Output_orfs (\%nucl_sequences, $minimum_length, "forward", *OUTPUT);
Reverse_complement (\%nucl_sequences);
Output_orfs (\%nucl_sequences, $minimum_length, "reverse", *OUTPUT);
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


sub Output_orfs {

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
	my $orf_num = 0;
	foreach my $key(sort keys(%$sequences_ref)) {
	    $this_seq = $$sequences_ref{$key};
		if ($this_seq =~ /(>.*?)\n(.*)\n/) {
			$header = $1;
			$seq = $2;
		}
		#print "Translating\.\.\. $header\n";
		$this_seq =~ tr/[actgn]/[ACTGN]/;
		while ($this_seq =~ /[ATCG]*?(ATG\w+?$)/) {
			my $this_part = $1;
			my $seq_length = length ($this_part);
			if ($seq_length < $min_length) {		
				last;
			}
			my $orf = Print_codons($this_part);
			if (length($orf) >= $min_length) {
				$orf_num++;
				if ($header =~ /^(>\S+)( .*)$/) {
					$ID = $1;
					$description = $2;
				} elsif ($header =~ /^(>\S+$)/) {
					$ID = $1;
					$description = "";
				} else {
					die "The sequence is not in fasta format\n";
				}
				print $fh_output ($ID,"_",$orientation,$orf_num,$description,"\n");
				print $fh_output ("$orf\n");
				$this_seq = substr($this_part, 3);
			} else {
				$this_seq = substr($this_part, 3);
			}
		}
	}
	return;
		
}

#--------------------------------------------------------------------------

sub Print_codons {

#
#	Translates the nucleotide sequence passed in the argument by looking up 
#	the globally defined %DNAtoAA hash table for the genetic code and
#	returns the peptide sequence
#

	my $this_seq = shift;
	$this_seq =~ tr/[a-z]/[A-Z]/;
	my $nucl;
	for (my $y = 0; $y < (length($this_seq) - 3); $y += 3) {
		if (!defined $DNAtoAA{substr($this_seq, $y, 3)}) {
			#$pept .= "X";
			last;
		} else {
			$nucl .= substr($this_seq, $y, 3);
			if ($DNAtoAA{substr($this_seq, $y, 3)} eq '*') {last}
		}
	}
	return $nucl;
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
