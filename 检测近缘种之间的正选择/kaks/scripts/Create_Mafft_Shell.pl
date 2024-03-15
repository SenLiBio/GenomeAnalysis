use warnings;
use strict;
my $homolog_dir=shift;
my @files=<$homolog_dir/*>;
my $output_sh=shift;
open (O,">$output_sh");

my $outdir=shift;
system("mkdir $outdir");
foreach my $file (@files){
    $file=~/(\d+)\.pro/;
    my $sign=$1;
    my $outfile="$outdir/mafft$sign";
    print O "mafft --maxiterate 1000 --localpair $file >$outfile\n";
}
close O;