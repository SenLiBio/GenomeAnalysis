my @axts=glob("*axt");
foreach my $axt (@axts) {
        ($kaks=$axt)=~s/axt/kaks/;
        `~/software/KaKs_Calculator2.0/bin/Linux/KaKs_Calculator -i $axt -o $kaks -c 1 -m YN`;
}