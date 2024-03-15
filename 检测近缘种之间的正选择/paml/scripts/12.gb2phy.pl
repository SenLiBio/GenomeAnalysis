my @gbs=glob("*gb");

foreach my $gb (@gbs) {
	open GB,"<","$gb";
	while (<GB>) {
		if (/^>/) {
			$_=~m/\>(\w{3})\|/;
			$spe=$1;
		} else {
			s/\s+//g;
			$seq{$spe}.=$_;
		}
	}
	$phy="$gb.phy";
	open PHY,">","$phy";
	my @spes=sort keys %seq;
	my $spe_num=@spes;
	my $length=length $seq{$spes[1]};
	print PHY "$spe_num $length\n";
	foreach my $spe (@spes) {
		print PHY "$spe  $seq{$spe}\n";
	}
	undef %seq;
}