use autodie;
use strict;

my $PosGFFile=shift @ARGV;
open POSGF,"<","$PosGFFile";
my $OrthogroupFile=shift @ARGV;
open GF,"<","$OrthogroupFile";
open ATHGENELIST,">","Ath.genelist.txt";
my @Spes=@ARGV;
open SPEGENELIST,">","Spe.genelist.txt";

my @Positive_Selected_OGs;
while (<POSGF>) {
		s/\s+$//;
		my $OG=(split/\s+/)[0];
		$OG=~s/pal2aln_/OG/;
		$OG=~s/\-gb.*$//;
		push @Positive_Selected_OGs,$OG;
}
close POSGF;

while (<GF>) {
	s/\s+$//;
	my $Line=$_;
	my @eles=split/\s+/;
	my $OG=shift @eles;
	$OG=~s/://;
	if (grep {$OG eq $_} @Positive_Selected_OGs) {
		foreach my $ele (@eles) {
			if ($ele=~m/Ath\|(.*)/) {
				my $Ath_Gene=$1;
				print ATHGENELIST "$Ath_Gene\n";
			} elsif ($ele=~m/(?<Gene>(?<Spe>\w{3})\|.*)/) {
				my $Spe_Gene=$+{Gene};
				my $Spe=$+{Spe};
				if (grep {$Spe eq $_} @Spes) {
					print SPEGENELIST "$Spe_Gene\n";
				}
			}
		}
	}
}
close GF;