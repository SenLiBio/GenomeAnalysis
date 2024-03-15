die "Usage:perl $0 Orthogroups.GeneCount.tsv Orthogroups.txt allpep.fa\n" if @ARGV != 3;



open GC,"<","$ARGV[0]";
open OG,"<","$ARGV[1]";
open FA,"<","$ARGV[2]";


L: while (<GC>) {
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

while (<FA>) {
        s/\s+$//;
        if (/^\>/) {
                s/^\>//;
                $Gene=$_;
        } else {
                $Seq{$Gene}.=$_;
        }
}

while (<OG>) {
        s/\s+$//;
        my @eles=split/\s+/;
        my $OG=shift @eles;
        $OG=~s/://;
        if (grep {$OG eq $_} @OGs) {
                open OUT,">","$OG.fa";
                foreach my $ele (@eles) {
                        print OUT "\>$ele\n$Seq{$ele}\n";
                }
                close OUT;
        }
}