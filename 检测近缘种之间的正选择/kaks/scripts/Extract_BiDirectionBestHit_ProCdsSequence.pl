# This script is used to extract PROTEIN and CDS sequences of besthit gene pairs.
# The extracted protein and cds sequences will be outputed to two independent folders.

die "Usage: perl $0 [Protein 1] [Protein 2] [Cds 1] [Cds 2] [BestHit file]\n" if @ARGV<5;

open PEP,"<","$ARGV[0]";
open CDS,"<","$ARGV[1]";
open PAIR,"<","$ARGV[2]";

my ($GeneName,%PepSeq,%CdsSeq);
while (<PEP>) {
	s/\s+$//;
	if (/^>/) {
		s/^>//;
		$GeneName=$_;
	} else {
		$PepSeq{$GeneName}.=$_;
	}
}
close PEP;


while (<CDS>) {
	s/\s+$//;
	if (/^>/) {
		s/^>//;
		$GeneName=$_;
	} else {
		$CdsSeq{$GeneName}.=$_;
	}
}
close CDS;

`mkdir ProSeqs`;
`mkdir CdsSeqs`;
my $FileCount;
while (<PAIR>) {
	$FileCount++;
	my $ProOut="ProSeqs\/$FileCount\.pro";
	my $CdsOut="CdsSeqs\/$FileCount\.cds";
	s/\s+$//;
	my ($Gene1,$Gene2)=split/\s+/;
	open PRO,">","$ProOut";
	open CDS,">","$CdsOut";
	print PRO ">$Gene1\n$PepSeq{$Gene1}\n>$Gene2\n$PepSeq{$Gene2}\n";
	print CDS ">$Gene1\n$CdsSeq{$Gene1}\n>$Gene2\n$CdsSeq{$Gene2}\n";
	close PRO;
	close CDS;
}