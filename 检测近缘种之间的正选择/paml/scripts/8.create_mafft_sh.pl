use warnings;
use strict;
my $homolog_dir=shift;
my @files=<$homolog_dir/*>;
my $output_sh=shift;
open (O,">$output_sh");

my $outdir=shift;
system("mkdir $outdir");
foreach my $file (@files){
    $file=~/OG(\d+)/;
    my $sign=$1;
    my $outfile="$outdir/mafft$1";
    print O "mafft --maxiterate 1000 --localpair $file >$outfile\n";
}
close O;

`split -l 500 -d -a 1 mafft.sh mafft_`;

my @shs=glob("mafft_*");
foreach my $sh (@shs) {
        `nohup bash $sh >>mafft.log 2>>mafft.err &`;
}