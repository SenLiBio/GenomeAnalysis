while (<>) {
	if (/^>/) {
		$id=$_;
	} else {
		s/\s+$//;
		my $seq=$_;
		my $gap=grep /X/,(split //,$_);
print STDERR "$gap\n";
		my $length=length $seq;
		print "$id$seq\n" if ($gap/$length<=0.2);
	}
}
