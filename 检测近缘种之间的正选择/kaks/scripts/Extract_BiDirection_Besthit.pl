# This script is used to extract inters-pecies bidirection besthit gene pairs .
# The input file is the inter-species blast output in outfmt 6.


die "Usage: perl $0 [BlastOut] >Outputfile\n" if @ARGV<1;
use strict;
use autodie;

open BLAST,"<","$ARGV[0]";
my %BestHit;
my %BestHitScore;
my @Genes;

Line: while (<BLAST>) {
	s/\s+$//;
	my ($RefGene,$QueryGene,$Score)=(split/\s+/)[0,1,-1];
	next Line if $RefGene eq $QueryGene;
	my $SpeOfRefGene=(split/\./,$RefGene)[0];
	my $SpeOfQueryGene=(split/\./,$QueryGene)[0];
	next Line if $SpeOfRefGene eq $SpeOfQueryGene;
	if (!exists $BestHit{$RefGene} || $Score>=$BestHitScore{$RefGene}) {
		$BestHit{$RefGene}=$QueryGene;
		$BestHitScore{$RefGene}=$Score;
	} else {
		next Line;
	}
	push @Genes,$RefGene unless (@Genes~~$RefGene);
}
close BLAST;

my %BiDirectionBestHit;
foreach my $Gene (@Genes) {
	if ($BestHit{$BestHit{$Gene}}=$Gene) {
		my ($Key,$Value)=(sort ($Gene,$BestHit{$Gene}))[0,1];
		$BiDirectionBestHit{$Key}=$Value;
	} 
}

foreach my $Gene (keys %BiDirectionBestHit) {
	print "$Gene\t$BiDirectionBestHit{$Gene}\n";
}