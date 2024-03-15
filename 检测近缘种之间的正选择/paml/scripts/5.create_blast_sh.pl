my @fas=glob("*fa");
open BLS,">","blast.sh";
foreach my $fa (@fas) {
		$fa=~m/(OG\d+)\.fa/;
		my $prefix=$1;
		print BLS "blastp -query $fa -db $prefix -outfmt 6 -num_threads 5 -evalue 1e-5 -out $prefix.out\n"
}

`split -l 1000 -d -a 1 blast.sh blast_`;

my @shs=glob("blast_*");
foreach my $sh (@shs) {
        `nohup bash $sh &`;
}