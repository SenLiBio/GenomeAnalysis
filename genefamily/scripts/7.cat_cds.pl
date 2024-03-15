use warnings;
use strict;
use Bio::SeqIO;
my $trimmeddir=shift;
my @pepfiles=<$trimmeddir/*>;
my $pepdir=shift;
system("mkdir $pepdir");

my %pep;
my $count=0;
foreach my $file (sort @pepfiles){
    my $fa=Bio::SeqIO->new(-file=>$file,-format=>'fasta');
    while (my $seq_obj=$fa->next_seq){
        my $pepid=$seq_obj->display_name;
        my $pepseq=$seq_obj->seq;
    $pepid=~/(\w+)\|(.+)/;
    my $species=$1;
    $pep{$species}.=$pepseq;
    }
    $count++;
}
print "total files:$count\n";

my @species=sort keys %pep;
foreach my $species (@species){
    my $output="$pepdir/$species.fa";
    open (O,">$output");
    print O ">$species\n$pep{$species}\n";
    close O;
}