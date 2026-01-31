#!usr/local/bin/perl
use warnings;
use strict;

open IN,"customAnnot.gbff";

my $specie=0;
my $acc="";
my $line="";
my $UTR5_STOP="";
my $CDS_STOP="";
my %acc_utr5_stop=();
my %acc_cds_stop=();

while (<IN>){
    chomp $_;
    if ($_ =~ /SOURCE\s+Mus musculus .+/){
	$specie=1;
	}
    if ($_=~/ACCESSION/){
	my @line = split /\s+/, $_;
	$acc = $line[1];
    }
    if ($specie==1 && $_ =~ /CDS\s+\d+\..\d+/){
	my @line=split /\s+/, $_;
	my $line2=$line[2];
	my @line2 = split /\.\./, $line2;
	$UTR5_STOP =$line2[0]-1;
	$CDS_STOP = $line2[1]
 }
    if ($specie==1 && $_=~/\/\//){
	$acc_utr5_stop{$acc}=$UTR5_STOP;
	$acc_cds_stop{$acc}=$CDS_STOP;
	$acc="";
	$specie=0;
	$UTR5_STOP="";
	$CDS_STOP="";
    }
}
close IN;

##

open IN , "customFasta.fa";

open ALL_STATS, ">customAnnot.txt";
print ALL_STATS "id\tUTR5_len\tcds_stop\ttotal_length\n";


my $utr5_stop=0;
my $cds_stop=0;
my $id="";
my $tmp_utr="";
my $seq="";
my $tracker=0;
my $tmp_utr5="";

my $tmp_cds="";
my $tmp_utr3="";


while (<IN>){
    chomp $_;
    if($_ =~/\>.*/ && $tracker==0){
	my @line = split /\s/, $_;
	$id = $line[0];
        $id =~s/\>//g;
	$id =~s /(\w\_\w+)\.\w+/$1/g;
    }
    unless ($_ =~/\>.*/){
	$seq=join "", $seq,$_;
	$tracker=1;
    }
    if($_ =~/\>.*/ && $tracker==1){
	$utr5_stop=$acc_utr5_stop{$id};
	$cds_stop=$acc_cds_stop{$id};
	unless($utr5_stop eq ""){
	    $tmp_utr5 = substr ($seq, 0 , $utr5_stop);
	    
	    ##tricky with start and stops in relation to positions in the vector
	    $tmp_cds = substr($seq, $utr5_stop, $cds_stop-$utr5_stop);
	    ##Get full length;
	    my @tmpSeq= split "", $seq;
	    my $tmpLength = @tmpSeq;
	    my $utr3Length= $tmpLength - $cds_stop;
	    $tmp_utr3 = substr($seq, $cds_stop, $utr3Length);

 
	    print ALL_STATS "$id\t$utr5_stop\t$cds_stop\t$tmpLength\n";
	}
	
	my @line = split /\s/, $_;
	$id = $line[0];
        $id =~s/\>//g;
	$id =~s /(\w\_\w+)\.\w+/$1/g;

	$tracker =1;
	$utr5_stop = 0;
	$cds_stop = 0;
	$seq="";
	$tmp_utr5="";
	$tmp_cds="";
	$tmp_utr3="";
    }
}
