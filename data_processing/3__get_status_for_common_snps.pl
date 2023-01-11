#!/usr/bin/perl

use strict;
use Getopt::Long;
my $config;

main();
my $hash;

sub main{
    my $file = $ARGV[0];
    open (FILE, $file) or die "cannot open $file : $!";
    while(<FILE>){
	chomp;
	my @d = split "\t";
	my $nCol = scalar(@d);
	my ($c,$p,$i,$ref,$alt,$qual,$filt,$inf,$form)=($d[0],$d[1],$d[2],$d[3],$d[4],$d[5],$d[6],$d[7],$d[8]);
	next if($_=~/^##/);

	if ($_=~/^#CHROM/){
	    foreach my $n (9..$nCol-1){
		my $lpID = $d[$n];
		$hash->{orderedLPid}{$n}=$lpID;
	    }
	}

	else{
	    my $var = join("\t",$c,$p,$i,$ref,$alt);
	    $hash->{variant}{$var}++;

	    foreach my $n (9..$nCol-1){
		my @details = split ":", $d[$n];
		my ($gt,$numReadsVar,$numReadsRef,$gq) = ($details[0],$details[6],$details[4],$details[3]);
		$hash->{genotype}{$var}{$n}=$gt;
		$hash->{numReadsVar}{$var}{$n}=$numReadsVar;
		$hash->{GQ}{$var}{$n}=$gq;
		$hash->{numReadsRef}{$var}{$n}=$numReadsRef;
		$hash->{allIDsGT}{$n}{$var}=$gt;
		$hash->{variants}{$var}++;
	    }
	}
    }
    close FILE;

    open (OUT, ">genotype-by-ID.txt");
    print OUT "chr\tpos\tid\tref\talt\tcolumnID\tLPid\tgenotype\tnumReadsVar\tnumReadsRef\tGQ\n";
    foreach my $var (keys %{$hash->{variant}}){
	foreach my $col (keys %{$hash->{orderedLPid}}){
	    print OUT "$var\t$col\t$hash->{orderedLPid}{$col}\t$hash->{genotype}{$var}{$col}\t$hash->{numReadsVar}{$var}{$col}\n";
	}
    }

    my @vars = ("chr11\t89284793\t.\tG\tA","chr11\t89178528\t.\tC\tA","chr11\t89177653\t.\tC\tT");
    open (OUT2,">status-by-ID-testing.txt");
    print OUT2 "LPid\tArg402Gln\tSer192Tyr\trs4547091\n";
    foreach my $col (keys %{$hash->{allIDsGT}}){
	my $LPid = $hash->{orderedLPid}{$col};
	print OUT2 "$LPid";
	foreach my $variant (@vars){
	    my $status;
	    if (defined $hash->{allIDsGT}{$col}{$variant}){
		($status)=($hash->{allIDsGT}{$col}{$variant});
	    }
	    print OUT2 "\t$status";
	}
	print OUT2 "\n";

    }

}
#########################


##########################
sub configure{
    my ($config)=@_;
    my $args = shift;
    my $config= {};
    GetOptions($config, "proband=s",
	       "mother=s",
	       "father=s",
	       "help|h!", )
	|| warn "error : $!\n";
    my @f = ("proband","mother","father");
    foreach my $F (@f){
	if (!defined $config->{$F}){
	    print "\n ERROR:\n a file needs to be provided with the --$F flag\n\n";
	    usage();
	}
    }

    return($config);

}
####################
sub usage{
    die "\n\n Not using script correctly \n
Usage:
--proband
--mother
--father


";

}
