my @fas=glob("*fa");
open SHS,">","mkblastdb.sh";
foreach my $fa (@fas) {
        $fa=~m/(OG\d+)\.fa/;
        my $prefix=$1;
        print SHS "makeblastdb -in $fa -out $prefix -dbtype prot\n";
}

`split -l 1000 -d -a 1 mkblastdb.sh mkblastdb_`;

my @shs=glob("mkblastdb_*");
foreach my $sh (@shs) {
        `nohup bash $sh >>mkblastdb.log 2>>mkblastdb.err &`;
}