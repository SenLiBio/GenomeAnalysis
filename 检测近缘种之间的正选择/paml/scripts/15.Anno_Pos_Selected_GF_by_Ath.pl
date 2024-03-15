## This script is used to annotate the positive selected gene family using the gene description information of A. thaliana.
 
die "Usage: perl $0 positive_selected.txt Orthogroups.txt ath.gene.discription ath.tf\n" if @ARGV!=4;

use autodie;
use strict;

open POSGF,"<","$ARGV[0]";
open GF,"<","$ARGV[1]";
open ANNO,"<","$ARGV[2]";
open TF,"<","$ARGV[3]";
open POSANNO,">","positive_selected.gf.Ath_anno.txt";
open POSTF,">","positive_selected.gf.Ath_tf.txt";
open TFCOUNT,">","positive_selected.gf.Ath_tf.count.txt";

my @Positive_Selected_OGs;
while (<POSGF>) {
	s/\s+$//;
	my $OG=(split/\s+/)[0];
	$OG=~s/pal2aln_/OG/;
	$OG=~s/\-gb.*$//;
	push @Positive_Selected_OGs,$OG;
}
close POSGF;

my %Anno;
while (<ANNO>) {
	next if /^Locus/;
	s/\s+$//;
	my @eles=split/\s+/;
	shift @eles;
	my $Gene=shift @eles;
	my $Anno=join(" ",@eles);
	$Anno{$Gene}=$Anno;
}
close ANNO;

my %TF;
while (<TF>) {
	next if /^TF/;
	s/\s+$//;
	my ($Gene,$TF)=(split/\s+/)[0,2];
	$TF{$Gene}=$TF;
}
close TF;

my %TFNum;
while (<GF>) {
	s/\s+$//;
	my $Line=$_;
	my @eles=split/\s+/;
	my $OG=shift @eles;
	$OG=~s/://;
	if (grep {$OG eq $_} @Positive_Selected_OGs) {
		$Line=~m/Ath\|(.*?)\s+/;
		my $Ath_Gene=$1;
		print POSANNO "$OG\t$Ath_Gene\t$Anno{$Ath_Gene}\n";
		print POSTF "$OG\t$Ath_Gene\t$TF{$Ath_Gene}\n";
		$TFNum{$TF{$Ath_Gene}}++ if (exists $TF{$Ath_Gene});
	}
}
close GF;
close POSTF;
close POSANNO;


foreach my $TF (sort {$TFNum{$b}<=>$TFNum{$a}} keys %TFNum) {
	print TFCOUNT "$TF\t$TFNum{$TF}\n";
}