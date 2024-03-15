use 5.010;
use strict;
use autodie;



open GFF,"<","$ARGV[0]";
my (%GeneNames,%GeneCounts,%CDSCounts,%Position,%Direction);
GFFLINE: while (<GFF>) {
    my $line=$_;
    s/\s+$//;
    my ($Chrom,$Type,$Sta,$End,$Direction,$Stat)=(split/\s+/)[0,2,3,4,6,8];
    $Stat=~m/=(.*);/;
    my $GeneName=$1;
    if ($Type eq 'mRNA') {
        $GeneCounts{$Chrom}=-1 if !exists $GeneCounts{$Chrom};
        $GeneCounts{$Chrom}++;
        push @{$GeneNames{$Chrom}},$GeneName;
        $Direction{$GeneName}=$Direction;
    } else {
        my @Position=($Sta,$End);
        @Position=sort {$a<=>$b} @Position;
        $CDSCounts{$GeneName}=-1 if (!exists $CDSCounts{$GeneName});
        $CDSCounts{$GeneName}++;
        $Position{$GeneName}->[$CDSCounts{$GeneName}]=[@Position];
    }
}
close GFF;


open PILE,"<","$ARGV[1]";
my $LineCounts=0;
my (%CDSSequence,%FinishedCDSCounts,%FinishedGeneCounts);
PILELINE: while (<PILE>) {
    my $Line=$_;
    $LineCounts++;

    s/\s+$//;
    my ($Chrom,$Position,$Ref,$Depth,$Base,$BaseQuality,$MappingQuality)=(split/\s+/);


    $FinishedGeneCounts{$Chrom}=0 if (!exists $FinishedGeneCounts{$Chrom});
    next PILELINE if ($FinishedGeneCounts{$Chrom}>$GeneCounts{$Chrom});

    my $GeneName=@{$GeneNames{$Chrom}}[$FinishedGeneCounts{$Chrom}];
    next PILELINE if (!defined $GeneName);
	$FinishedCDSCounts{$GeneName}=0 if (!exists $FinishedCDSCounts{$GeneName});
    if ($FinishedCDSCounts{$GeneName}==$CDSCounts{$GeneName}+1) {
        $FinishedGeneCounts{$Chrom}++;
        my %CDS;
        if ($Direction{$GeneName} eq '+') {
            for (my $i = 0; $i <= $CDSCounts{$GeneName}; $i++) {
                $CDS{$GeneName}.=$CDSSequence{$GeneName}->[$i];
            }
        } else {
            for (my $i = $CDSCounts{$GeneName}; $i >= 0; $i--) {
#print "$i\t$GeneName\t$CDSSequence{$GeneName}->[$i]\n";
                my $Reverse_Complementary_Sequence = reverse_complementarity ($CDSSequence{$GeneName}->[$i]);
                $CDS{$GeneName}.=$Reverse_Complementary_Sequence;
            }
        }

        my $length=length $CDS{$GeneName};
        warn "The length of CDS sequence of gene \'$GeneName\' is $length, not an integer multiple of three at $LineCounts\n" if ($length % 3 != 0);
        open OUT,">","$GeneName\.fa";
        print OUT ">$GeneName\n$CDS{$GeneName}\n";

    } else {
        my $GeneName=@{$GeneNames{$Chrom}}[$FinishedGeneCounts{$Chrom}];
        my @pos=@{$Position{$GeneName}->[$FinishedCDSCounts{$GeneName}]};
        my ($sta,$end)=@pos[0,1];
        if ($Position<$sta) {
            next PILELINE;
        } else {
            if ($Depth<=5) {
#print "$LineCounts\t$GeneName\t'N'\t$FinishedCDSCounts{$GeneName}\t$CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}]\n";
                $CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}].='N';
                $FinishedCDSCounts{$GeneName}++ if $Position==$end;
            } else {

                $Base=~s/\$|\^[\000-\177]{1}//g;

                my @IndelCounts=($Base=~/[\+\-]([0-9]+)[ACGTNacgtn]+/g);
                my (@Bases,@BaseQualities,@MappingQualities);
                for (@IndelCounts) {
                    $Base=~s/[\+\-][0-9]+[ACGTNacgtn]{$_}//;
                }
                @Bases=split//,$Base;
                @BaseQualities=split//,$BaseQuality;
                @MappingQualities=split//,$MappingQuality;
                if (@Bases!=@BaseQualities || @Bases!=@MappingQualities) {
                    print "$Line\n";
                    print "\@Bases=","scalar(@Bases)",",\@BaseQualities=","scalar(@BaseQualities)","\n";
                    die("die at line $LineCounts \@read != \@baseq \n");
                }
                    my $EffectiveDepth=0;
                    my %AlleleCount='';
                  BASE: for (1..@Bases) {
                        my $Base=shift @Bases;
                        if ($Base eq '.' || $Base eq ',') {
                            $Base=$Ref;
                        } else {
                            $Base=uc($Base);
                        }
                        my $TheBaseQuality=shift @BaseQualities;
                        my $TheMappingQuality=shift @MappingQualities;
                        next BASE if (ord($TheBaseQuality)-33<=20);
                        next BASE if (ord ($TheMappingQuality)-33<=20);
                        $EffectiveDepth++;
                        $AlleleCount{$Base}++;
                  }

                if ($EffectiveDepth>=200 || $EffectiveDepth <= 5) {
                    $CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}].='N';
#print "$LineCounts\t$GeneName\t'N'\t$FinishedCDSCounts{$GeneName}\t$CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}]\n";
                } else {
                  my @ATCG=sort {$AlleleCount{$b}<=>$AlleleCount{$a}} keys %AlleleCount;
                    if (@ATCG>=3) {
                        if ($AlleleCount{$ATCG[2]}/$EffectiveDepth>=0.05) {
                            $CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}].='N';
#print "$LineCounts\t$GeneName\t'N'\t$FinishedCDSCounts{$GeneName}\t$CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}]\n";
                       } elsif ($AlleleCount{$ATCG[1]}/$EffectiveDepth>=0.1 && $AlleleCount{$ATCG[1]}>=2) {
                            my @Alleles=($ATCG[0],$ATCG[1]);
                            my $DegeneratedBase=&bases_degeneration (@Alleles);
#print "$LineCounts\t$GeneName\t$DegeneratedBase\t$FinishedCDSCounts{$GeneName}\t$CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}]\n";
                            $CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}].=$DegeneratedBase;
                        } else {
#print "$LineCounts\t$GeneName\t$ATCG[0]\t$FinishedCDSCounts{$GeneName}\t$CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}]\n";
                            $CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}].=$ATCG[0];
                        }
                    } elsif (@ATCG==2) {
                        if ($AlleleCount{$ATCG[1]}/$EffectiveDepth>=0.1 && $AlleleCount{$ATCG[1]}>=2) {
                            my @Alleles=($ATCG[0],$ATCG[1]);
                            my $DegeneratedBase=&bases_degeneration (@Alleles);
                            $CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}].=$DegeneratedBase;
#print "$LineCounts\t$GeneName\t$DegeneratedBase\t$FinishedCDSCounts{$GeneName}\t$CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}]\n";
                        } else {
                            $CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}].=$ATCG[0];
#print "$LineCounts\t$GeneName\t$ATCG[0]\t$FinishedCDSCounts{$GeneName}\t$CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}]\n";
                        }
                    } else {
#print "$LineCounts\t$GeneName\t$ATCG[0]\t$FinishedCDSCounts{$GeneName}\t$CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}]\n";
                        $CDSSequence{$GeneName}->[$FinishedCDSCounts{$GeneName}].=$ATCG[0];
                    }
                }
                $FinishedCDSCounts{$GeneName}++ if $Position==$end;
            }
        }
    } 
}







sub bases_degeneration {
    my @Bases=sort map {uc($_)} @_;
    if (@Bases==0) {
        warn "There's no avail bases";
    }
    my $Base=join("",@Bases);
    my %DegeneratedBase=(
            'AG' => 'R',
            'CT' => 'Y',
            'AC' => 'M',
            'GT' => 'K',
            'AT' => 'W',
            'CG' => 'S',
            'ACT' => 'H',
            'CGT' => 'B',
            'ACG' => 'V',
            'AGT' => 'D',
            'ACGT' => 'N',
    );
    return $DegeneratedBase{$Base};
}

sub reverse_complementarity {
    my $Seq= $_[0];
    $Seq=~s/A/a/g;
    $Seq=~s/C/b/g;
    $Seq=~s/T/A/g;
    $Seq=~s/G/C/g;
    $Seq=~s/a/T/g;
    $Seq=~s/b/G/g;
    $Seq=join("",(reverse (split//,$Seq)));
    return $Seq;
}
