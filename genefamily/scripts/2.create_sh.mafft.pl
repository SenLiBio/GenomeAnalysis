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