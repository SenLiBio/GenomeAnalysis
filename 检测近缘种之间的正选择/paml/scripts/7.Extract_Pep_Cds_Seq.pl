use autodie;
use strict;

open CDS,"<","$ARGV[0]";
open PEP,"<","$ARGV[1]";
open GENELIST,"<","$ARGV[2]";

my ($Gene_Name,%Cds,%Pep);
while (<CDS>) {
	s/\s+$//;
	if (/^\>/) {
		s/^\>//;
		$Gene_Name=$_;
	} else {
		$Cds{$Gene_Name}.=$_;
	}
}
close CDS;


while (<PEP>) {
	s/\s+$//;
	if (/^\>/) {
		s/^\>//;
		$Gene_Name=$_;
	} else {
		$Pep{$Gene_Name}.=$_;
	}
}
close PEP;


while (<GENELIST>) {
	s/\s+$//;
	my @eles=split/\s+/;
	my $OG=shift @eles;
	open PEPSEQ,">","$OG.pep";
	open CDSSEQ,">","$OG.cds";
	foreach my $Gene (@eles) {
		print PEPSEQ "\>$Gene\n$Pep{$Gene}\n";
		print CDSSEQ "\>$Gene\n$Cds{$Gene}\n";
		warn "Cannot find the CDS sequence of Gene $Gene in $OG\n" if !exists $Cds{$Gene};
		warn "Cannot find the protein sequence of Gene $Gene in $OG\n" if !exists $Pep{$Gene};		
	}
}
close GENELIST;

