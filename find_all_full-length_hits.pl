#!/usr/bin/perl -w
#
#
use warnings;
use strict;
#
#
#
#--------------------------------------------------------------------------
# find_full-length_hits.pl
#--------------------------------------------------------------------------
#
# Parses the NCBI-BLAST+ output file into a tabular form with
# following annotations in order:
# Query, Query Length, Matching Subject, Subject Length, Bit Score,
# E-value, Percent Match, Number of Matching aa/Nucleotides, Length of
# the HSP. (You can add more if you want by modifying the sub.)
#
# Extracts the top hits.
#
# Usage: ./find_full-length_hits.pl blast_result_file > output_tab_file
#
# Modified for Ronny to extract details of "full-length" hits.
#
# @joewinnz
# Dec 2015
#
#
#--------------------------------------------------------------------------
#


# "Global" scalars
my $query_orf_start_count = 0;
my $subj_orf_start_count = 0;
my $matched_orf_start_count = 0;


my $storedResult = "";
my $noHitsCount = 0;
my $hitsCount = 0;
my $totalQuery = 0;

# Change the value of the following scalar to get the desired number
# of HSPs from hits (e.g. $numberHSPs = 5; to get 5 top HSPs).
my $HSP_count = 1000;
my $first_query = 1;

print "Query\tQuery Length\tQuery Start\tStart aa\tQuery End\tSubject\tSubject Start\tSubject start aa\tSubject End\tPercent query covered\tFrame\n";
while (<>) {
# This snippet removes presence of ">" in the middle of lines. This is
# sometimes true for query/subject descriptions. Not needed for most files.
#	if (!(/^>/)){
#		s/>/ /g;
#	}
    
# Reads in the blast result file line by line and processes each query
    if (/^Query=/) {
    	$totalQuery++;
    	if ($first_query) {
    		$first_query = 0;
    		$storedResult = $_;
    		next;
    	} else {        
        	$storedResult .= ">\n";
        	printHSPs ($storedResult, $HSP_count);
        	if ($storedResult =~ /^.*No hits found.*/m) {
            	$noHitsCount++;
        	} else {
            	$hitsCount++;
        	}
        	$storedResult = "";
        	$storedResult = $_;
        }
    } else {
    	$storedResult .= $_;
    }
}
$storedResult .= ">\n";
printHSPs ($storedResult, $HSP_count);
if ($storedResult =~ /^.*No hits found.*/m) {
	$noHitsCount++;
} else {
	$hitsCount++;
}


print "\n\nTotal number of queries in the BLAST output file = $totalQuery \n";
print "Query starts at aa position 1 = $query_orf_start_count\n";
print "Subject aa \"M\" at hsp position 1 = $subj_orf_start_count\n";
print "Both Query and Subject start with \"M\" = $matched_orf_start_count\n";

exit;


#--------------------------------------------------------------------------
#Subroutine
#--------------------------------------------------------------------------
#
#
sub printHSPs {
	# This is where the real hard gruff is done.
	# Extracts and print annotations for HSPs contained in each blast result.
    # Needs two parameters: (1) a string containing blast result from single
    # query, and (2) an integer - number of HSPs to process
    # Only collects info from top HSP from the subjects with multiple HSPs.
	
    my ($oneResult, $numberToPrint) = @_;
    my ($query) = $oneResult =~ /^Query= (.*)$/m;
    my ($queryLength) = $oneResult =~ /\n\nLength=(\d+)/m;
    my @HSPs = ();
    my @lines = split /\n/, $oneResult;
    my $temp = "";
    my $seqFound = 0;
	my $collecting = 0;
		
    foreach my $line(@lines) {
    	$line .= "\n";
    	if ($line =~ /^>.*/) {   	
    		my $header = $line;
    		if ($collecting) {
    			push @HSPs, $temp;
    			$collecting = 0;
    			$temp = "";
    		}
    		if (!$collecting) {
    			$collecting = 1;
    			$temp .= $line;
    		}
    	} else {
    		$temp .= $line;
    	}
    }
    		
    my $i = @HSPs;
    my $j = 0;
    $numberToPrint = $i if ($i < $numberToPrint);
    if ($i == 0) {
        #print "$query\t$queryLength\tNo hits found\n";
        #print "No hits found\n";
    } else {
        for ($j = 0; $j < $numberToPrint; $j++) {
            my ($HSP) = $HSPs[$j];
            my ($subjHeader) = $HSP =~ /^> (.*)Length/ms;
            $subjHeader =~ s/\n/ /g;
            $subjHeader =~ s/\s+/ /g;
            my ($annotation) = $subjHeader =~ /^\S+\s(.*)/;
            my ($subjLength) = $HSP =~ /^Length=(\d+)\n\n Score/m;
            my ($bits) = $HSP =~ /^\s+Score =\s+([\d\.\+\-e]+) bits/m;
            $bits = "1$bits" if $bits =~ /^e/;
            my ($eVal) = $HSP =~ /Expect\S*? = ([\d\.\+\-e]+)/;
            $eVal = "1$eVal" if $eVal =~ /^e/;
            my ($frame) = $HSP =~ /^ Frame = (\S+)\n/m;
            my ($match, $total, $percent)
                = $HSP =~ /Identities = (\d+)\/(\d+) \((\d+)%\)/;
            my @hsp_lines = split (/\n/, $HSP);
            my $subj_start = 0;
            my $query_start = 0;
            my $subj_end = 0;
            my $query_end = 0;
            my $first_score_line = 1;
            my $orf_start = "NO";
            my $subj_orf_start = "NO";
            my $query_start_aa = "";
            my $subj_start_aa = "";
            foreach my $hsp_line(@hsp_lines) {
            	if ($hsp_line =~ /^\sScore/) {
            		if ($first_score_line == 1) {
            			$first_score_line = 0;
            		} else {
            			last;
            		}
            	} 
            	next if ($hsp_line =~ /^\s+/);
            	if ($hsp_line =~ /^Query\s+(\d+)\s+/) {
            		if ($query_start == 0) {
            			$hsp_line =~ /^Query\s+(\d+)\s+/;
            			$query_start = $1;
            			if ($query_start == 1) {
            				$orf_start = "YES";
            				$hsp_line =~ /^Query\s+\d+\s+(\S+)\s+\d+$/;
            				my $query_seq_line = $1;
            				$query_start_aa = substr ($query_seq_line, 0, 1);
            				$query_orf_start_count++;
            			}
            			$hsp_line =~ /^Query\s+\d+.*\s+(\d+)$/;
            			$query_end = $1;
            		} else {
            			$hsp_line =~ /^Query\s+\d+.*\s+(\d+)$/;
            			$query_end = $1;
            		}
            	}
            	if ($hsp_line =~ /^Sbjct\s+(\d+)\s+/) {
            		if ($subj_start == 0) {
            			$hsp_line =~ /^Sbjct\s+(\d+)\s+/;
            			$subj_start = $1;
            			$hsp_line =~ /^Sbjct\s+\d+\s+(\S+)\s+\d+$/;
            			my $subj_seq_line = $1;
            			$subj_start_aa = substr ($subj_seq_line, 0, 1);
            			if ($subj_start_aa eq "M") {
            				$subj_orf_start_count++;
            				$subj_orf_start = "YES";
            			}
            			$hsp_line =~ /^Sbjct\s+\d+.*\s+(\d+)$/;
            			$subj_end = $1;
            		} else {
            			$hsp_line =~ /^Sbjct\s+\d+.*\s+(\d+)$/;
            			$subj_end = $1;
            		}
            	}
            }
            my $percent_covered = (($query_end - $query_start +1)/$queryLength)*100;
            if (($orf_start eq "YES") && ($subj_orf_start eq "YES")) {
            	print "$query\t$queryLength\t$query_start\t$query_start_aa\t$query_end\t$subjHeader\t$subj_start\t$subj_start_aa\t$subj_end\t";
            	printf("%.2f", $percent_covered);
            	print "\t$frame";
            	print "\n";
            	$matched_orf_start_count++;
            }           
        }
    }    
}
#--------------------------------------------------------------------------
