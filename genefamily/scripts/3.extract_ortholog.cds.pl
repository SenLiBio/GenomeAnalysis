use warnings;
use strict;
use Bio::SeqIO;

my $orthomclfile=shift;
open (I,"<$orthomclfile");
my $allnucledir=shift;
my @nuclefiles=<$allnucledir/*>;
my $outdir=shift;
system("mkdir $outdir");
my $genefamilylength_check=shift;
open (F,">$genefamilylength_check");

my %nucle;
foreach my $nuclefile (@nuclefiles){
    my $fa=Bio::SeqIO->new(-file=>$nuclefile,-format=>'fasta');
    while (my $seq_obj=$fa->next_seq){
    my $nucle_name=$seq_obj->display_name;
    my $seq=$seq_obj->seq;
    $nucle{$nucle_name}=$seq;
    }
    print "$nuclefile\n";
}

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
    my $id=$2;
    $family{species}{$genefamily}{$species}++;
    $family{id}{$genefamily}{$species}=$gene;

    }
}
close I;

my %familylength;
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
        my $ortholog=$nucle{$id};
        my $length=length($ortholog);
        if ($length<150){
        $familylength{$genefamily}++;
        }
        print O ">$id\n$ortholog\n";  
    }
    close O;
    }
}

my @familylength=keys %familylength;
foreach my $familylength (@familylength){
    print F "$familylength\n";
}
close F;