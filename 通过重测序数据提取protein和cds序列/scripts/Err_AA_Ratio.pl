open IN1,"<","$ARGV[0]";
open IN2,"<","$ARGV[1]";
open LIST,">","genelist.ratio";
open HIST,">","genelist.ratio.list";

while (<IN1>) {
	if (/^>/) {
		s/\s+$//;
		s/^>//;
		$id=$_;
	} else {
		s/\s+$//;
		s/\*$//;
		my $seq1=$_;
		$seq1{$id}=$seq1;
	}
}
close IN1;

while (<IN2>) {
	if (/^>/) {
		s/\s+$//;
		s/^>//;
		$id=$_;
	} else {
		s/\s+$//;
		s/\*$//;
		my $seq2=$_;
		my $seq1=$seq1{$id};
		my @seq1=split//,$seq1;
		my @seq2=split//,$seq2;
		my $length=length $seq2;
		for (1..$length) {
			my $base1=shift @seq1;
			my $base2=shift @seq2;
			$SameBase++ if $base1 eq $base2;
		}
		my $ratio=$SameBase/$length;
		print LIST "$id\t$ratio\n";
		$num{$ratio}++;
		$hist=int($ratio*100+0.5)/100;
		$freq{$hist}++;
		$sum++;
		undef $SameBase;
	}
}

map {my $ratio=$freq{$_}/$sum;$accu+=$ratio;print HIST "$_\t$ratio\t$accu\n"} (sort {$b<=>$a} keys %freq);

