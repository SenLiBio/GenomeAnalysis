use warnings;
use strict;
my $pepdir=shift;
my @mafftfiles=<$pepdir/*>;
my $cdsdir=shift;
my $outdir=shift;
system ("mkdir $outdir");
my $output_sh=shift;
open(O,">$output_sh");

my $pal2alnpath="~/software/pal2nal.v14/pal2nal.pl";
foreach my $pepfile (@mafftfiles){
    $pepfile=~/mafft(\d+)/;
    my $sign=$1;
    print O "perl $pal2alnpath $pepfile $cdsdir/paltoaln$sign -output fasta >$outdir/pal2aln_$sign\n";
}
close O;