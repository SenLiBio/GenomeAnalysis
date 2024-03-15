use strict;
my $pepdir=shift;
my @files=<$pepdir/pal*>;
my $logdir=shift;
system ("mkdir $logdir");
my $output_sh=shift;
open(O,">$output_sh");

my $gblockspath="~/software/Gblocks_0.91b/Gblocks";
foreach my $pepfile (@files){
    print O "$gblockspath $pepfile -t=c\n";
}
close O;