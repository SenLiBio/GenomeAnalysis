use autodie;
use strict;

my @blastouts=glob("*out");


foreach my $blastout (@blastouts) {

        open BLASTOUT,"<","$blastout";
        open GENELIST,">>","Best_Hit.genelist";

        my (%Best_Hit_Score,%Best_Hit_Line,%Best_Hit_Query_Gene);

        L: while (<BLASTOUT>) {
                s/\s+$//;
                my $Line=$_;
                my ($Ref,$Query,$Score)=(split/\s+/)[0,1,-1];
                $Query=~m/^(\w{3})\|/;
                my $Query_Spe=$1;
                next L if (($Ref=~m/Ath/)==0);
                if ($Score>=$Best_Hit_Score{$Query_Spe} || !exists $Best_Hit_Score{$Query_Spe}) {
                        $Best_Hit_Score{$Query_Spe}=$Score;
                        $Best_Hit_Line{$Query_Spe}=$Line;
                        $Best_Hit_Query_Gene{$Query_Spe}=$Query;
                }
        }

        $blastout=~m/(\w+)\.out/;
        my $OG=$1;
        print GENELIST "$OG\t";
        my $Best_Hit_Blastout="$OG.besthit";
        open BHOUT,">","$Best_Hit_Blastout";

        foreach my $Query_Spe (keys %Best_Hit_Score) {
                print BHOUT "$Best_Hit_Line{$Query_Spe}\n";
                print GENELIST "$Best_Hit_Query_Gene{$Query_Spe}\t";
        }

        print GENELIST "\n";

        undef %Best_Hit_Score;
        undef %Best_Hit_Line;
        undef %Best_Hit_Query_Gene;
}