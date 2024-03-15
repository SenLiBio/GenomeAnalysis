die "Usage: perl $0 Orthogroups.GeneCount.tsv\n" if @ARGV!=1;

L: while (<>) {
        if (/^Or/) {next L}
        s/\s+$//;
        my @eles=split/\s+/;
        my $OG=shift @eles;
        my $ath=$eles[0];
        next if $ath !=1;
        for ($m=0;$m<$#eles;$m++) {
                next L if ($eles[$m]>3 || $eles[$m]==0);
        }
        push @OGs,$OG;
}
`mkdir lowcopy`;
my @fas=glob("*fa");
foreach my $fa (@fas) {
        $fa=~m/(OG\d+)\.fa/;
        my $OG=$1;
        `cp $fa lowcopy` if $OG~~@OGs;
}