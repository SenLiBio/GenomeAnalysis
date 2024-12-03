open PHY,"<$ARGV[0]";
open OUT,">","mcmctree.phy";

@lines=(<PHY>);
$header=shift @lines;
my ($species,$base)=split/\s+/,$header;
$codonbase=$base/3;

for (1..$species) {
        my $seq=shift @lines;
        my ($spe,$codon)=split/\s+/,$seq;
        my @bases=split//,$codon;
        my $count=1;
        foreach my $base (@bases) {
                $codon1{$spe}.="$base" if $count%3==1;
                $codon2{$spe}.="$base" if $count%3==2;
                $codon3{$spe}.="$base" if $count%3==0;
                $count++;
        }
}

foreach my $num (1..3) {
        print OUT "$species $codonbase\n";
        my $hashname="codon" . "$num";
        foreach my $spe (sort keys %codon1) {
                print OUT "$spe  ${$hashname}{$spe}\n";
        }
}
