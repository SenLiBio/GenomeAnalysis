#Usage: perl *pl Length_threshlod IN.fa passed.fa unpassed.fa Length_of_scaffold.txt
use Bio::SeqIO;

$length_threshold=$ARGV[0];
open CACHE,">>CACHE";
open PASSED,">>$ARGV[2]";
open FILTERED,">>$ARGV[3]";
open LEN,">>$ARGV[4]";


my $in=Bio::SeqIO->new(-format=>'fasta',
                    -file=>"$ARGV[1]");

while (my $seq=$in->next_seq) {
                $id=$seq->id;
                $length=$seq->length;
                $base=$seq->seq;
#print "$id\t$length\t$base\n";
                if ($length<=$length_threshold) {
                        print FILTERED ">$id\n$base\n";
                } else {
                        print CACHE ">$id\n$base\n";
                }
}
close FILTERED;
close CACHE;

my $in=Bio::SeqIO->new(-format=>"fasta",
                                        -file=>"CACHE");

while (my $seq=$in->next_seq) {
                $id=$seq->id;
                $length=$seq->length;
                $base=$seq->seq;

                $id{$length}=$id;
                $base{$id}=$base;
}
@length=sort {$b<=>$a} keys %id;
foreach $length (@length) {
        print PASSED ">$id{$length}\n$base{$id{$length}}\n";
        print LEN "$id{$length}\t$length\n";
}

close LEN;
close PASSED;
close CACHE;
unlink CACHE;