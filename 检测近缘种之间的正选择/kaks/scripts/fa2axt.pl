use strict;
use warnings;
use Bio::SeqIO;
use autodie;

my @fas=glob("pal*");
foreach my $fa (@fas) {
	my $axt="$fa.axt";
	open OUT,">","$axt";
	print OUT "$fa\n";
	my $obj=Bio::SeqIO->new (-format=>'fasta',-file=>"$fa");
	while (my $seq=$obj->next_seq) {
		my $sequence=$seq->seq;
		print OUT "$sequence\n";
	}
	close OUT;
}