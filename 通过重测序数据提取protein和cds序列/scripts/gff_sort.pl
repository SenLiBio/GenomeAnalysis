use 5.010;
use strict;
use autodie;

my (@Chroms,%Genes,%Position,%Num,%mRNA,%CDS);

open GFF,"<","$ARGV[0]";
open SORTEDGFF,">","$ARGV[1]";
while (<GFF>) {
        my $line=$_;
        s/\s+$//;
        my ($Chrom,$Type,$Sta,$End,$Stat)=(split/\s+/)[0,2,3,4,-1];
        $Stat=~m/=(.*);/;
        my $GeneName=$1;
        $mRNA{$GeneName}=$line if $Type eq "mRNA";
        next if $Type ne 'CDS';
        my @Position=($Sta,$End);
        @Position=sort {$a<=>$b} @Position;
        push @Chroms,$Chrom unless $Chrom~~@Chroms;
        push @{$Genes{$Chrom}},$GeneName unless $GeneName~~$Genes{$Chrom};
        $Num{$GeneName}++;
        $Position{$GeneName}->[$Num{$GeneName}]=[@Position];
        $CDS{$GeneName}->[$Num{$GeneName}]=$line;
}

@Chroms=map {$_->[0]}
                sort {$a->[1]<=>$b->[1]}
                map {my $Chrom=$_;$Chrom=~m/(\d+)$/;my $temp=$1;[$Chrom,$1]} @Chroms;

foreach my $Chrom (@Chroms) {
        my %Numbers;
        my @GeneNames=@{$Genes{$Chrom}};
        foreach my $GeneName (@GeneNames) {
				my @temp=sort {$Position{$GeneName}->[$a][0]<=>$Position{$GeneName}->[$b][0]} (1..$Num{$GeneName});   
				$Numbers{$GeneName}=[@temp];
		}
        @GeneNames=sort {$Position{$a}->[$Numbers{$a}->[0]][0]<=>$Position{$b}->[$Numbers{$b}->[0]][0]} @GeneNames;
        foreach my $GeneName (@GeneNames) {
print "$GeneName\t$Position{$GeneName}->[@{$Numbers{$GeneName}}[0]][0]\n";
	       print SORTEDGFF "$mRNA{$GeneName}";
                foreach my $Number (@{$Numbers{$GeneName}}) {
                        print SORTEDGFF "$CDS{$GeneName}->[$Number]";
                }
        }
}
























