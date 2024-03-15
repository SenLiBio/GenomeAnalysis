open M0,"<","$ARGV[0]";
open MA,"<","$ARGV[1]";

my @fas;

my %M0_lnL;
while (<M0>) {
	s/\s+$//;
	$_=~m/(?<fa>pal.*)\.m0:lnL.*(?<lnL>\-\d+).*/;
	my $fa=$+{fa};
	push @fas,$fa;
	my $M0_lnL=$+{lnL};
	$M0_lnL{$fa}=$M0_lnL;
}

my %Ma_lnL;
while (<MA>) {
	s/\s+$//;
	$_=~m/(?<fa>pal.*)\.ma:lnL.*(?<lnL>\-\d+).*/;
	my $fa=$+{fa};
	my $Ma_lnL=$+{lnL};
	$Ma_lnL{$fa}=$Ma_lnL;
}

my %delta_lnLs;
foreach my $fa (@fas) {
	my $delta_lnL=abs($Ma_lnL{$fa}-$M0_lnL{$fa});
	$delta_lnLs{$fa}=2*$delta_lnL;
}

my @fas=sort {$delta_lnLs{$b}<=>$delta_lnLs{$a}} keys %delta_lnLs;
open POS,">positive_selected.txt";
open NON,">non_positive_selected.txt";
foreach my $fa (@fas) {
	if ($delta_lnLs{$fa}>=2.7055) {
		print POS "$fa\t$delta_lnLs{$fa}\n";
	} else {
		print NON "$fa\t$delta_lnLs{$fa}\n";
	}
}