use warnings;
use strict;

use Bio::SeqIO;
my $mafftdir=shift;
my @pepfiles=<$mafftdir/*>;
my $cdsdir=shift;
my $outdir=shift;
system ("mkdir $outdir");

foreach my $pepfile (@pepfiles){
    $pepfile=~/mafft(\d+)/;
    my $sign=$1;
    my $cdsfile="$cdsdir/OG$sign.cds";

    my %pep;
    my $count=0;
    my $fa_P=Bio::SeqIO->new(-file=>$pepfile,-format=>'fasta');
    while (my $seq_obj=$fa_P->next_seq){
    my $pep_id=$seq_obj->display_name;
    $pep_id=~/^(\w+)\|(.+)/;
    my $species=$1;
    my $pepid=$2;
    $count++;
    $pep{$count}=$species;
    }

    my %nucle;
    my $fa_c=Bio::SeqIO->new(-file=>$cdsfile,-format=>'fasta');
    while (my $seq_obj=$fa_c->next_seq){
    my $nuc_id=$seq_obj->display_name;
    my $seq=$seq_obj->seq;
    $nuc_id=~/^(\w+)\|(.+)/;
    my $species=$1;
    my $nucid=$2;
    $nucle{$species}{$nuc_id}=$seq;
    }

    my $output="$outdir/paltoaln$sign";
    open(O,">$output");
    my @order=sort {$a<=>$b} keys %pep;
    foreach my $order (@order){
    my $species=$pep{$order};
    my @id=keys %{$nucle{$species}};
    foreach my $id (@id){
        print O ">$id\n$nucle{$species}{$id}\n";
    }
    }
    close O;
}