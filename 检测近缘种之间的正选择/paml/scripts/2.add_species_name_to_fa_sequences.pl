use autodie;

my $name=$ARGV[0];
$name=~s/\.fa//;
open IN,"<$ARGV[0]";
open OUT,">bak";

while (<IN>) {
        if (/^>/) {
                s/\s+$//;
                $a=(split/\s+/,$_)[0];
                $a=~s/^>//;
                $a=">$name|" . "$a";
                print OUT "$a\n";
        } else {
                print OUT;
        }
}
close IN;
close OUT;
`mv bak $ARGV[0]`;