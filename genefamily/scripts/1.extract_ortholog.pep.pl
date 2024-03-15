use strict;
use Bio::SeqIO;


my $orthomclfile=shift;
open (I,"<$orthomclfile");
my $allpepfiles=shift;       
my $outdir=shift;
system("mkdir $outdir");

my %pep;
my $fa=Bio::SeqIO->new(-file=>$allpepfiles,-format=>'fasta');
while (my $seq_obj=$fa->next_seq){
    my $pep_name=$seq_obj->display_name;
    my $seq=$seq_obj->seq;
    $pep{$pep_name}=$seq;
}
print "pep.fa has done\n";


my %family;
while (<I>){
    chomp;
    my $line=$_;
    my @inf=split/\s+/,$line;
    my $genefamily=shift(@inf);
    $genefamily=~s/://;
    foreach my $gene (@inf){
        $gene=~/^(\w+)\|(.+)/;
        my $species=$1;
        my $id=$gene;
        $family{species}{$genefamily}{$species}++;
        $family{id}{$genefamily}{$species}=$gene;

    }
}
close I;

my @familyes=keys %{$family{species}};
foreach my $genefamily (@familyes){
    my @species=sort keys %{$family{species}{$genefamily}};
    my $count=0;
    foreach my $species (@species){
        $count+=$family{species}{$genefamily}{$species};
    }
    my $mark=scalar(@species);
    if ($count==$mark && $mark==10){

        my $output="$outdir/$genefamily";
        open (O,">$output");
        foreach my $species (@species){
            my $id=$family{id}{$genefamily}{$species};
            my $ortholog=$pep{$id};
            print O ">$id\n$ortholog\n";
        }
        close O;
    }
}