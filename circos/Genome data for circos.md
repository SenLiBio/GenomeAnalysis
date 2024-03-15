# Genome data for circos

# 1. karyotype.txt

## 1.1 统计染色体长度

```bash
perl sca_length_and_cumulative_percentage.pl PREFIX.fa >PREFIX.length
```

- sca_length_and_cumulative_percentage.pl
    
    ```perl
    use strict;
    use warnings;
    use Bio::SeqIO;
    
    die "Usage: perl sca_length_and_cumulative_percentage.pl PREFIX.fa >PREFIX.sca_length\n" if @ARGV==0;
    my $fa=Bio::SeqIO->new(-file=>$ARGV[0],-format=>'fasta');
    
    my %length;
    my $total_length;
    while (my $obj=$fa->next_seq()) {
    	my $id=$obj->id;
    	my $length=$obj->length;
    	$length{$id}=$length;
    	$total_length+=$length;
    }
    
    print "# rank\tchr\tlength\tcumulative_percentage\n";
    my @ids=sort {$length{$b}<=>$length{$a}} keys %length;
    my $cumulative_length;
    my $sort_index;
    foreach my $id (@ids) {
    	$sort_index++;
    	$cumulative_length+=$length{$id};
    	my $cumulative_percentage=100*$cumulative_length/$total_length;
    	print "$sort_index\t$id\t$length{$id}\t";
    	printf "%.2f%%\n",$cumulative_percentage;
    }
    ```
    

## 1.2 做karyotype.txt

可以自由选取两种配色方式：

1. 用ucsc的配色，比较单个基因组的时候推荐用这种方式；
2. 比较两个基因组的时候，可以每个基因组的染色体用一种颜色表示；

```bash
perl sca_length2circors_karyotype.pl PREFIX.length [species prefix] [SCA_NUM] >PREFIX.karyotype
```

- sca_length2circors_karyotype.pl
    
    ```perl
    use strict;
    use warnings;
    use autodie;
    die "Usage: perl sca_length2circors_karyotype.pl PREFIX.LENGTH PREFIX SCA_NUM >PREFIX.karyotype\n" if @ARGV!=3;
    open LEN,"<","$ARGV[0]";
    my $prefix=$ARGV[1];
    my $sca_num=$ARGV[2];
    
    while (<LEN>) {
    	next if /^#/;
    	chomp;
    	my ($rank,$sca,$length)=(split/\s+/)[0,1,2];
    	die if $rank>$sca_num;
    	$length-=1;
    	my $color=$rank%18;
    	$color=($color==0) ? ("18") : ("$color");
    	print "chr - $prefix$rank $prefix$rank 0 $length chr$color\n";
    }
    close LEN;
    ```
    

---

# 2. heatmap

基因密度、repeat可以用heatmap来画

## 2.1 基因密度

基因密度的文件可以用tbtools直接获取，这里对这个文件进行处理，制作heatmap的输入文件。

```bash
perl genedensity2circosheatmap.pl PREFIX.length PREFIX.genedensity [SCA_NUM] [species prefix] >PREFIX.genedensity.4circos
```

- genedensity2circosheatmap.pl
    
    ```perl
    use strict;
    use warnings;
    
    die "Usage: perl $0 PREFIX.length PREFIX.genedensity [target_sca_num] [SPECIES PREFIX] PREFIX.genedensity.4circos >prefix.gd.range\n" if @ARGV!=5;
    
    open LENGTH,"<$ARGV[0]";
    open GD,"<$ARGV[1]";
    my $target_sca_num=$ARGV[2];
    my $prefix=$ARGV[3];
    open OUT,">$ARGV[4]";
    my @target_scas;
    my %scaname;
    while (<LENGTH>) {
    	next if /^#/;
    	my ($index,$chr)=(split/\s+/,$_)[0,1];
    	push @target_scas,$chr if ($index<=$target_sca_num);
    	$scaname{$chr}="$prefix$index";
    }
    close LENGTH;
    my @gds;
    my @gd_on_target_sca;
    my @lines;
    while (<GD>) {
    	chomp;
    	my ($sca,$sta,$end,$gd)=(split/\s+/)[0..3];
    	$end-=1;
    	$gd=int $gd;
    	push @gds,$gd;
    	if ($sca~~@target_scas) {
    		$sca=$scaname{$sca};
    		my @name=($sca,$sta,$end,$gd);
    		my $name=\@name;
    		push @lines,$name;
    	}
    }
    
    @lines=sort {$a->[0] cmp $b->[0] or $a->[1]<=>$b->[1]} @lines;
    map {print OUT "@$_\n"} @lines;
    @gds=sort {$a<=>$b} @gds;
    my $max=pop @gds;
    my $min=shift @gds;
    print "$max\t$min\n";
    ```
    

## 2.2 重复序列密度

重复序列的密度可以直接从repeatmasker的输出文件统计里面的“N”的个数（当然更推荐的是将重复序列用“X”来mask，然后统计“X”的个数，因为基因组文件里面本身就可能含有"N")。

```bash
perl maskedgenome2circorsheatmap.pl PREFIX.length PREFIX.fa.masked [SCA_NUM] [block_size] [species prefix] >PREFIX.repeatdensity.4circos
```

- maskedgenome2circorsheatmap.pl
    
    ```perl
    use strict;
    use warnings;
    use Bio::SeqIO;
    die "Usage: perl $0 PREFIX.length PREFIX.fa.masked [target_sca_num] [block_size] [species prefix] PREFIX.repeatcontent >PREFIX.repeat.range\n" if @ARGV!=6;
    
    open LENGTH,"<$ARGV[0]";
    my $fastq=$ARGV[1];
    my $target_sca_num=$ARGV[2];
    my $block_size=$ARGV[3];
    my $prefix=$ARGV[4];
    open OUT,">$ARGV[5]";
    my @target_scas;
    my %scaname;
    my %length;
    while (<LENGTH>) {
    	next if /^#/;
    	my ($index,$chr,$length)=(split/\s+/,$_)[0,1,2];
    	push @target_scas,$chr if ($index<=$target_sca_num);
    	$scaname{$chr}="$prefix$index";
    	$length{$chr}=$length;
    }
    close LENGTH;
    my @values;
    my $maskedgenome=Bio::SeqIO->new(-file=>$fastq,-format=>'fasta');
    while (my $obj=$maskedgenome->next_seq()) {
    	my $id=$obj->id;
    	my $seq=$obj->seq;
    	if ($id~~@target_scas) {
    		my @seq=split//,$seq;
    		my $length=$length{$id};
    		my $block_num=int ($length/$block_size);
    		foreach my $num (1..$block_num) {
    			my %num;
    
    			for (1..$block_size) {
    				my $base=shift @seq;
    				$num{$base}++;
    			}
    			my $percent=$num{'N'}/$block_size;
    			$percent=int ($percent*=100);
    			push @values,$percent;
    			my $sta=($num-1)*$block_size;
    			my $end=$num*$block_size-1;
    			print OUT "$scaname{$id}\t$sta\t$end\t$percent\n";
    			undef %num;
    		}
    		my %num;
    		foreach my $base (@seq) {
    			$num{$base}++;
    		}
    		my $sta=$block_num*$block_size;
    		my $end=$length{$id}-1;
    		my $percent=$num{'N'}/$block_size;
    		$percent=int ($percent*=100);
    		push @values,$percent;
    		print OUT "$scaname{$id}\t$sta\t$end\t$percent\n";
    		undef %num;
    	}
    }
    @values=sort {$a<=>$b} @values;
    my $max=pop @values;
    my $min=shift @values;
    print "$max\t$min\n";
    ```
    

## 2.3 gc密度

```perl
perl [genome2circosgccontent.pl](http://genome2circosgccontent.pl/) PREFIX.length PREFIX.genome.fasta [SCA_NUM] [block_size] [species prefix] PREFIX.gccontent >PREFIX.gc.range
```

- genome2circosgccontent.pl
    
    ```perl
    use strict;
    use warnings;
    use Bio::SeqIO;
    die "Usage: perl $0 PREFIX.length PREFIX.genome.fasta [target_sca_num] [block_size] [species prefix] PREFIX.gccontent >PREFIX.gc.range\n" if @ARGV!=6;
    
    open LENGTH,"<$ARGV[0]";
    my $fastq=$ARGV[1];
    my $target_sca_num=$ARGV[2];
    my $block_size=$ARGV[3];
    my $prefix=$ARGV[4];
    open OUT,">$ARGV[5]";
    my @target_scas;
    my %scaname;
    my %length;
    while (<LENGTH>) {
    	next if /^#/;
    	my ($index,$chr,$length)=(split/\s+/,$_)[0,1,2];
    	push @target_scas,$chr if ($index<=$target_sca_num);
    	$scaname{$chr}="$prefix$index";
    	$length{$chr}=$length;
    }
    close LENGTH;
    my @values;
    my $maskedgenome=Bio::SeqIO->new(-file=>$fastq,-format=>'fasta');
    while (my $obj=$maskedgenome->next_seq()) {
    	my $id=$obj->id;
    	my $seq=$obj->seq;
    	if ($id~~@target_scas) {
    		my @seq=split//,$seq;
    		my $length=$length{$id};
    		my $block_num=int ($length/$block_size);
    		foreach my $num (1..$block_num) {
    			my %num;
    
    			for (1..$block_size) {
    				my $base=shift @seq;
    				$num{$base}++;
    			}
    			my $percent=($num{'G'}+$num{'C'})/$block_size;
    			$percent=int ($percent*=100);
    			push @values,$percent;
    			my $sta=($num-1)*$block_size;
    			my $end=$num*$block_size-1;
    			print OUT "$scaname{$id}\t$sta\t$end\t$percent\n";
    			undef %num;
    		}
    		my %num;
    		foreach my $base (@seq) {
    			$num{$base}++;
    		}
    		my $sta=$block_num*$block_size;
    		my $end=$length{$id}-1;
    		my $percent=($num{'G'}+$num{'C'})/@seq;
    		$percent=int ($percent*=100);
    		push @values,$percent;
    		print OUT "$scaname{$id}\t$sta\t$end\t$percent\n";
    		undef %num;
    	}
    }
    @values=sort {$a<=>$b} @values;
    my $max=pop @values;
    my $min=shift @values;
    print "$max\t$min\n";
    ```
    

---

# 3. links

## 3.1 从共线性文件提取links

注意，这里的gff文件用的是mcscan运行时用的同一个gff文件。

```bash
perl collinearity2circoslink.pl PREFIX.gff PREFIX.collinearity >PREFIX.inks.txt
```

- collinearity2circoslink.pl
    
    ```perl
    # Usage: perl collinearity2circoslink.pl PREFIX.gff PREFIX.collinearity >PREFIX.inks.txt
    
    use strict;
    use autodie;
    die "Usage: perl collinearity2circoslink.pl PREFIX.gff PREFIX.collinearity >PREFIX.inks.txt\n" if @ARGV==0;
    open GFF,"<","$ARGV[0]";
    open CLO,"<","$ARGV[1]";
    
    my (%chr_of_gene,%sta_of_gene,%end_of_gene);
    while (<GFF>) {
    	chomp;
    	my ($chr,$gene,$sta,$end)=split/\s+/;
    	$chr_of_gene{$gene}=$chr;
    	$sta_of_gene{$gene}=$sta;
    	$end_of_gene{$gene}=$end;
    }
    close GFF;
    
    while (<CLO>) {
    	next if /^#/;
    	chomp;
    	$_=~s/^\s*\d+\-\s*\d+:\s+//;
    	my ($gene1,$gene2)=(split/\s+/,$_)[0,1];
    	print "$chr_of_gene{$gene1} $sta_of_gene{$gene1} $end_of_gene{$gene1} $chr_of_gene{$gene2} $sta_of_gene{$gene2} $end_of_gene{$gene2}\n"; 
    }
    close CLO;
    ```