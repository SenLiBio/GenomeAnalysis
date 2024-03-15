# MCScanX+Kaks_calculator2

scripts: MCScanX+Ka%2027bf7/for_MCScanX.zip

# 1.分析流程

## 1.1 建库

```bash
makeblastdb -in bg.fa -out bg -dbtype prot
```

---

## 1.2 比对

```bash
blastp -query bg.fa -db bg -outfmt 6 -num_threads 5 -evalue 1e-5 -out bg.out
```

---

## 1.3 筛选

**filter_blast_results_by_evalue_identity.pl**

```perl
use 5.010;
use strict;
use autodie;
use warnings;

die "USAGE: perl filter_blast_results_by_evalue_identity.pl [identity] [e-value] [INPUTFILE] [OUTPUTFILE] \n" if (@ARGV < 1);

my $identity_threshold=$ARGV[0];
my $evalue_threshold=$ARGV[1];
open IN,"<$ARGV[2]" || die "Cannot open $ARGV[2]:$!";
open OUT,">$ARGV[3]" || die "Cannot open $ARGV[3]:$!";

while (<IN>) {
	my $line=$_;
	my ($identity,$evalue)=(split/\s+/,$line) [2,10];
	next if ($identity <= $identity_threshold || $evalue >= $evalue_threshold);
	print OUT "$line";
}
close IN;
close OUT;
```

```bash
perl filter_blast_results_by_evalue_identity.pl 40 1e-5 bg.out bg.blast
```

## 1.4 制作mcscanx用的gff文件

<aside>
💡 注意：这一步要注意基因名、染色体名称格式，并且在脚本中根据情况进行调整
**gff_to_MCScanXgff.pl**

```perl
use 5.010;
use strict;
use autodie;
use warnings;

die "USAGE: perl gff_to_MCScanXgff.pl [two-letter species contraction] [INPUTFILE] [OUTPUTFILE] \n" if (@ARGV < 1);

my $species_name_contraction=$ARGV[0];

open IN,"<$ARGV[1]" || die "Cannot open $ARGV[1]:$!";
open OUT,">$ARGV[2]" || die "Cannot open $ARGV[2]:$!";

while (<IN>) {
	if (/mRNA/) {
		chomp;
		my $line=$_;
		my ($sca,$sta,$end,$info)=(split/\s+/,$line)[0,3,4,8];
		$sca=~/(?<sca_num>\d+$)/;
		my $sca_num=$+{sca_num};
		$info=~/ID=(?<gene_name>\S+?);/;
		print OUT "$species_name_contraction$sca_num\t$+{gene_name}\t$sta\t$end\n";
	} else {
		next;
	}
}

close IN;
close OUT;
```

```bash
perl gff_to_MCScanXgff.pl bg Bgy.gff bgy.gff
```

---

## 1.5 运行

将blast文件和gff文件放在同一文件夹下，如“exampledir”，运行mcscanx

```bash
MCscanx exampledir/bgy
```

---

## 1.6 进一步统计

如果对两个物种进行种间比对，然后要检测种间的共线性，可以进行进一步统计

<aside>
💡 注意extract_interspecies_blocks.pl脚本运行可能有warning，不影响输出即可
**extract_interspecies_blocks.pl**

```perl
use 5.010;
use autodie;
use strict;

die "USAGE: perl extract_interspecies_blocks.pl [INPUTFILE] [OUTPUTFILE] \n" if (@ARGV < 1);

open IN,"<$ARGV[0]" || die "Cannot open $ARGV[0]:$!";
open OUT,">$ARGV[1]" || die "Cannot open $ARGV[1]:$!";

my @inputlines=<IN>;
for (1..8) {
	my $line=shift @inputlines;
	print OUT "$line";
}

shift @inputlines;
my $line=shift @inputlines;
$line=~/(\d+$)/;
my $num_of_all_genes=$1;
shift @inputlines;

my ($judge,%interspecies_block_name);
my $num_of_interspecies_genes=0;
my $num_of_interspecies_blocks=-1;
my $tsv;
my $num_of_gene_pairs_in_the_block;
my ($species1,$species2);
my @collinear_genes;

foreach my $line (@inputlines) {
	if ($line =~ m/^##/) {
		$line=~m/(\w{2})\d+\&(\w{2})\d+/;
		$species1=$1;
		$species2=$2;
		if ($species1 ne $species2) {
			$judge=1;
			$num_of_interspecies_blocks++;
			$line=~s/Alignment \d+/Alignment $num_of_interspecies_blocks/;
			$interspecies_block_name{$num_of_interspecies_blocks}=$line;
			$num_of_gene_pairs_in_the_block=0;
			$tsv->[$num_of_interspecies_blocks][0]=$line;
		} else {
			$judge=0;
		}
	} else {
		if ($judge==0) {
			next;
		} else {
			$line=~/^ *\d+\- *\d+:\s+(?<gene1>.*?)\s+(?<gene2>.*?)\s+/;
			push @collinear_genes,$+{gene1},$+{gene2};
			$num_of_gene_pairs_in_the_block++;
			$line=~s/^ *\d+/ $num_of_interspecies_blocks/;
			$tsv->[$num_of_interspecies_blocks][$num_of_gene_pairs_in_the_block]=$line;
		}
	}
}

close IN;

my %num;
map {$num{$_}++} @collinear_genes;
@collinear_genes=keys %num;
$num_of_interspecies_genes=@collinear_genes;

my $percentage_of_interspecies_genes=100*$num_of_interspecies_genes/$num_of_all_genes;
print OUT "# Number of collinear genes: $num_of_interspecies_genes, Percentage: ";
printf OUT "%5.2f\n", $percentage_of_interspecies_genes;
print OUT "# Number of all genes: $num_of_all_genes\n##########################################\n";

foreach my $num_of_block (0..$num_of_interspecies_blocks) {
	foreach my $num_of_gene_pairs_in_the_block (0..@{$tsv->[$num_of_block]}) {
		print OUT "$tsv->[$num_of_block][$num_of_gene_pairs_in_the_block]";
	}
}

close OUT;
```

**interspecies_genes_for_each_species.pl**

```perl
die "USAGE:perl interspecies_genes_for_each_species.pl [INPUTFILE] \n" if (@ARGV < 1); 
use strict;
my @species1_genes;
my @species2_genes;

while (<>) {
	next if /^#/;
	chomp;
	$_=~s/^\s*\d+\-\s*\d+:\s*//;
	my ($gene1,$gene2)=(split/\s+/,$_)[0,1];
	push @species1_genes,$gene1;
	push @species2_genes,$gene2;
}

my (%num_of_species1_gene,%num_of_species2_gene);
map {$num_of_species1_gene{$_}++} @species1_genes;
map {$num_of_species2_gene{$_}++} @species2_genes;

my $num_of_species1_gene=keys %num_of_species1_gene;
my $num_of_species2_gene=keys %num_of_species2_gene;
print "num_of_species1_gene:$num_of_species1_gene\nnum_of_species2_gene:$num_of_species2_gene\n";
```

```perl
perl extract_interspecies_blocks.pl bgy.bse.collinearity bgy.bse.interspecies.collinearity
perl interspecies_genes_for_each_species.pl bgy.bse.interspecies.collinearity
```

---

# 2. KaKs计算

## 2.1 提取共线性基因对应的cds和protein序列


💡 cds、pep分别是物种的cds和protein文件，根据需要选择是用种间共线性基因还是种内共线性基因来做

**extract_collinear_pep_cds_sequence.pl**

```perl
use strict;
use warnings;
use Bio::SeqIO;
use autodie;

die "Usage:perl $0 [CDS file] [PEP file] [COLLINEAR file] [CDS outdir] [PEP outdir]" if @ARGV<4;

my $cds=shift;
my $pep=shift;
my $col=shift;
my $cdsdir=shift;
my $pepdir=shift;
system ("mkdir $cdsdir");
system ("mkdir $pepdir");

my %pep;
my $pep=Bio::SeqIO->new (-file=>"$pep",-format=>"fasta");
while (my $seq=$pep->next_seq) {
	my $id=$seq->id;
	my $sequence=$seq->seq;
	$pep{$id}=$sequence;
}
print "$pep has done\n";

my %cds;
my $cds=Bio::SeqIO->new (-file=>"$cds",-format=>"fasta");
while (my $seq=$cds->next_seq) {
	my $id=$seq->id;
	my $sequence=$seq->seq;
	$cds{$id}=$sequence;
}
print "$cds has done\n";

open COL,"<","$col";
while (<COL>) {
	next if /^#/;
	$_=~/^ *(?<block>\d+)\- *(?<code>\d+):\s+(?<gene1>.*?)\s+(?<gene2>.*?)\s+/;
	my $outputpep="$pepdir/$+{block}_$+{code}.pep";
	my $outputcds="$cdsdir/$+{block}_$+{code}.cds";
	open PEPOUT,">","$outputpep";
	open CDSOUT,">","$outputcds";
	print PEPOUT ">$+{gene1}\n$pep{$+{gene1}}\n>$+{gene2}\n$pep{$+{gene2}}\n";
	print CDSOUT ">$+{gene1}\n$cds{$+{gene1}}\n>$+{gene2}\n$cds{$+{gene2}}\n";
	close PEPOUT;
	close CDSOUT;
}
close COL;
```

```bash
perl extract_collinear_pep_cds_sequence.pl cds pep collinearity cdsdir pepdir 
```

---

## 2.2 mafft比对蛋白质序列

**create_mafft_sh.pl**

```perl
use warnings;
use strict;

die "Usage: perl $0 [PEP dir] [SH output] [MAFFT outdir]" if @ARGV<3;

my $homolog_dir=shift;
my @files=<$homolog_dir/*>;
my $output_sh=shift;
open (O,">$output_sh");

my $outdir=shift;
system("mkdir $outdir");
foreach my $file (@files){
    $file=~/^.*\/(.*)\.pep/;
        my $sign=$1;
    my $outfile="$outdir/mafft$sign";
    print O "mafft --maxiterate 1000 --localpair $file >$outfile\n";
}
close O;
```

```bash
perl create_mafft_sh.pl pepdir mafft.sh mafftoutdir
bash mafft.sh
```

---

## 2.3 pal2nal 将比对好的蛋白序列转换成cds
💡 记得更改pal2nal软件的路径

**create_pal2nal.pl**

```perl
use warnings;
use strict;
my $pepdir=shift;
my @mafftfiles=<$pepdir/*>;
my $cdsdir=shift;
my $outdir=shift;
system ("mkdir $outdir");
my $output_sh=shift;
open(O,">$output_sh");

my $pal2alnpath="/data/users/lisen/software/pal2nal.v14";
foreach my $pepfile (@mafftfiles){
    $pepfile=~/mafft(\d+_\d+)/;
    my $sign=$1;
        my $cds="$sign.cds";
    print O "perl $pal2alnpath/pal2nal.pl $pepfile $cdsdir/$cds -nogap -output fasta >$outdir/pal2aln$sign\n";
}
close O;
```

```bash
perl create_pal2nal.pl mafftoutdir cdsdir pal2naloutdir pal2nal.sh
bash pal2nal.sh
```

---

## 2.4 将fa文件转换为axt文件

- **fa2axt.pl**
  
    ```perl
    use strict;
    use warnings;
    use Bio::SeqIO;
    use autodie;
    
    my @fas=glob("pal*");
    foreach my $fa (@fas) {
    	my $axt="$fa.axt";
    	open OUT,">","$axt";
    	print OUT "$fa\n";
    	my $obj=Bio::SeqIO->new (-format=>'fasta',-file=>"$fa");
    	while (my $seq=$obj->next_seq) {
    		my $sequence=$seq->seq;
    		print OUT "$sequence\n";
    	}
    	close OUT;
    }
    ```
    

```bash
cd pal2naloutdir
mkdir axts
mkdir fas
perl fa2axt.pl
mv *axt axts
mv pal* fas
```

---

## 2.5 计算kakscalculator计算kaks

💡 记得更改kakscalculator软件的路径

- kakscal.pl
  
    ```perl
    my @axts=glob("*axt");
    foreach my $axt (@axts) {
            ($kaks=$axt)=~s/axt/kaks/;
            `~/software/KaKs_Calculator2.0/bin/Linux/KaKs_Calculator -i $axt -o $kaks -c 1 -m YN`;
    }
    ```
    

```bash
cd axts
perl kakscal.pl
mkdir kaks
mkdir axts
mv *axt axts
mv *kaks kaks
```

---

## 2.6 统计kaks

- **kaks_sta.pl**
  
    ```perl
    my $prefix=$ARGV[0];
    my $output="$prefix.kakssta.txt";
    open OUT,">>","$output";
    print OUT "gene\tka\tks\tkaks\n";
    my @files=glob("*kaks");
    foreach my $file (@files) {
    	open IN,"<","$file";
    	while (<IN>) {
    		next if /^Sequence/;
    		chomp;
    		my ($gene,$ka,$ks,$kaks)=(split/\s+/)[0,2,3,4];
    		print OUT "$gene\t$ka\t$ks\t$kaks\n";
    	}
    }
    ```
    

```bash
cd kaks
perl kaks_sta.pl
```

---

## 2.7 ggplot2画图

- 针对单个species-pair的直方图
  
    ```r
    bgybse<-read.table(file = "bgy-bse.kakssta.txt",header=TRUE)
    bgybseplot<-ggplot(bgybse,aes(x=ks))+geom_line(stat="density")+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,by=.1))
    ggThemeAssist::ggThemeAssistGadget(bgybseplot)
    bcybseplot + theme(plot.title = element_text(size = 15, face = "bold.italic", colour = "red",hjust = 0.5, vjust = 2), panel.background = element_rect(fill = NA), plot.background = element_rect(fill = "gray92", linetype = "solid")) +labs(title = "bcy-bse")
    ```
    
- 多个species画峰峦图和直方图
    - kaks2rplot_input.pl
      
        ```perl
        use strict;
        use autodie;
        use warnings;
        
        my $out=$ARGV[0];
        open OUT,">>","$out";
        print OUT "species\tgene\tka\tks\tkaks\n";
        
        my @txts=glob("*txt");
        foreach my $txt (@txts) {
        	open IN,"<","$txt";
        	$txt=~/^(\w+\-\w+)\.kakssta\.txt/;
        	my $prefix=$1;
        	while (<IN>) {
        		next if /^gene/;
        		print OUT "$prefix\t$_";
        	}
        	close IN;
        }
        close OUT;
        ```
        
    - 峰峦图
      
        ```r
        library("ggplot2")
        library("ggridges")
        # 这个包的说明可以看一下网上的[blog](https://zhuanlan.zhihu.com/p/339605450)
        library("ggsci")
        colorplate<-pal_npg("nrc",alpha = 1)(10)
        # 配色可以看一下《R语言数据可视化之美》第一章，尤其是P70。
        try<-ggplot(kaks,aes(x=ks,y=species,fill=species))+geom_density_ridges()+scale_fill_manual(values=colorplate)+xlim(0,3)+theme(legend.position = "none")
        ggThemeAssist::ggThemeAssistGadget(try)
        try + theme(plot.title = element_text(size = 15,face = "bold.italic", colour = "red",hjust = 0.5, vjust = 1)) +labs(title = "Ks distribution")
        ```
        
    - 直方图
      
        ```r
        library("ggplot2")
        library("ggsci")
        colorplate<-pal_npg("nrc",alpha = 1)(10)
        # 配色可以看一下《R语言数据可视化之美》第一章，尤其是P70。
        ggplot(kaks,aes(x=ks,colour=species))+geom_density(size=1.2,alpha=.5)+xlim(0,3)+scale_color_manual(values=colorplate)+ theme(plot.title = element_text(size = 15,face = "bold.italic", colour = "red",hjust = 0.5, vjust = 1)) +labs(title = "Ks distribution")
        ```
        

---

## 2.8. shell and readme

- kaks.sh
  
    ```bash
    perl ../scripts/extract_collinear_pep_cds_sequence.pl cds pep collinearity cdsdir pepdir &&
    perl ../scripts/create_mafft_sh.pl pepdir mafft.sh mafftoutdir
    bash mafft.sh &&
    perl ../scripts/create_pal2nal.pl mafftoutdir cdsdir pal2naloutdir pal2nal.sh &&
    bash pal2nal.sh &&
    cd pal2naloutdir &&
    mkdir axts &&
    mkdir fas &&
    perl ../../scripts/fa2axt.pl &&
    mv *axt axts &&
    mv pal* fas &&
    cd axts &&
    perl ../../../scripts/kakscal.pl &&
    mkdir kaks &&
    mkdir axts &&
    mv *axt axts &&
    mv *kaks kaks &&
    cd kaks &&
    perl ../../../../scripts/kaks_sta.pl &&
    mv kaks_sta.txt ../../../../
    ```
    
- README
  
    ```bash
    Remember to modify the path to the software in 
    	1.create_mafft_sh.pl
    	2.create_pal2nal.pl
    	3.kakscal.pl
    and the name of collinearity file in 
    	1.kaks.sh
    
    The running directory should contain the following files:
    	1.cds file in fasta format named "cds";
    	2.protein file in fasta format named "pep";
    
    Create a copy of [kaks.sh] in the outputdir of MCScanX, and locate the scripts dir in the parent directory.
    Then run "nohup bash kaks.sh &"
    ```
    

---

# 3.WGD retention分析

## 3.1 统计每一个block中的ks中位数

<aside>
💡 一般来说最后统计出的ks中位数会分为两个明显的峰，对应两次wgd事件，一般选择近期的wgd事件进行分析。

</aside>

```bash
perl KsMedianOfCollinearBlocks.pl Prefix.kakssta.txt >Prefix.ks.median.txt
```

- KsMedianOfCollinearBlocks.pl
  
    ```perl
    use POSIX;
    my @Blocks;
    
    while (<>) {
        next if /^gene/;
        chomp;
        my ($GenePair,$ka,$ks,$kaks)=split/\s+/;
        my $Block=(split/_/,$GenePair)[0];
        $Block="\\@"."$Block";
        push @Blocks,$Block;
        push @{$Block},$ks;
    }
    
    my %num;
    map {$num{$_}++} @Blocks;
    @Blocks=sort keys %num;
    
    my %KsMedian;
    map {my $Median=&Median($_);($name=$_)=~s/\\@//;$KsMedian{$name}=$Median;} @Blocks;
    my @BlocksSortedByKsMedian=sort {$KsMedian{$a}<=>$KsMedian{$b}} keys %KsMedian;
    map {print "$_\t$KsMedian{$_}\n";} @BlocksSortedByKsMedian;
    
    sub Median {
        my @ksOfBlock=@{$_};
        @ksOfBlock=sort {$a<=>$b} @ksOfBlock;
        my $Median=($ksOfBlock[int($#ksOfBlock/2)]+$ksOfBlock[ceil($#ksOfBlock/2)])/2;
        return $Median;
    }
    ```
    

---

## 3.2 分离wgd保留和丢失的基因

<aside>
💡 1.根据ks中位数来筛选近期wgd产生的block
2.从block中根据位置信息筛选保留的基因和丢失的基因

</aside>

- Collinear2Retention_Loss.pl
  
    ```perl
    use 5.010;
    use Getopt::Long;
    use warnings;
    use autodie;
    
    my ($LowerLimit,$UpperLimit,$KsOfBlocks,$Collinear,$OutputPrefix,$GffFile);
    
    GetOptions (
            'L=f' => \$LowerLimit,
            'U=f' => \$UpperLimit,
            'K=s' => \$KsOfBlocks,
            'C=s' => \$Collinear,
            'G=s' => \$GffFile,
            'P=s' => \$OutputPrefix,
    ) or die $!;
    
    # Collinear blocks of recent WGD event 
    open KOB,"<$KsOfBlocks";
    my @BlockOfRecentWGD;
    
    while (<KOB>) {
            chomp;
            my ($NameOfBlock,$KsMedian)=split/\s+/;
            $NameOfBlock=~/(\d+)$/;
            my $NumOfBlock=$1;
            push @BlockOfRecentWGD,$NumOfBlock if ($KsMedian>=$LowerLimit && $KsMedian<=$UpperLimit);
    }
    my $NumOfBlockFromRecentWGD=@BlockOfRecentWGD+0;
    print "Number of Blcok From Recent WGD is: $NumOfBlockFromRecentWGD\n";
    
    close KOB;
    
    # Start and end positions of each gene
    open GFF,"<","$GffFile";
    my (%Chr,%StartPos,%EndPos);
    while (<GFF>) {
            chomp;
            my ($Chr,$Gene,$StartPos,$EndPos)=(split/\s+/);
            $StartPos{$Gene}=$StartPos;
            $EndPos{$Gene}=$EndPos;
            $Chr{$Gene}=$Chr;
    }
    close GFF;
    
    open PAIR,">$OutputPrefix.paired.gene";
    open COL,"<","$Collinear";
    
    my $NumOfGenePairFromRecentWGD;
    my @COL=(<COL>);
            for (my $m=1;$m<=11;$m++) {
                    shift @COL;
            }
    foreach my $line (@COL) {
            if ($line=~/^#/) {
                            $line=~/Alignment (\d+)/;
                            my $NumOfBlock=$1;
                            $line=~/(\w{2}\d+)\&(\w{2}\d+)/;
                            my $Chr1=$1;
                            my $Chr2=$2;
                            push @{$Chr1},"A$NumOfBlock" if $NumOfBlock~~@BlockOfRecentWGD;
                            push @{$Chr2},"B$NumOfBlock" if $NumOfBlock~~@BlockOfRecentWGD;
            } else {
            $line=~s/^\s+//;
            $line=~s/\s+$//;
            my ($GenePairName,$GenePair)=(split/:/,$line)[0,1];
            $GenePairName=~s/^\s+//;
            $GenePair=~s/^\s+//;
            my ($GeneA,$GeneB)=(split/\s+/,$GenePair)[0,1];
            say PAIR "$GeneA\t$GeneB";
            my $NumOfBlock=(split/-/,$GenePairName)[0];
            if ($NumOfBlock~~@BlockOfRecentWGD) {
                    my $GroupA="A$NumOfBlock";
                    my $GroupB="B$NumOfBlock";
                    $NumOfGenePairFromRecentWGD++;
                    push @{$GroupA},$GeneA;
                    push @{$GroupB},$GeneB;
                                    push @BlockGroupOfRecentWGD,$GroupA,$GroupB;
                                    push @GenesOfRecentWGD,$GeneA,$GeneB;
            }
        }
    }
    print "Number of Gene Paird From Recent WGD is: $NumOfGenePairFromRecentWGD\n";
    
    close COL;
    close PAIR;
    
    map {$elenum{$_}++} @BlockGroupOfRecentWGD;
    @BlockGroupOfRecentWGD=keys %elenum;
    
    open SINGLE,">$OutputPrefix.single.gene";
    
    my @SingleGenes;
    my (%StartPosOfBlock,%EndPosOfBlock);
    foreach my $Block (@BlockGroupOfRecentWGD) {
            my %NumOfGene;
            my (@StartPos,@EndPos);
            map {push @StartPos,$StartPos{$_};push @EndPos,$EndPos{$_};} (@{$Block});
            @StartPos=sort {$a<=>$b} @StartPos;
            @EndPos=sort {$b<=>$a} @EndPos;
            $StartPosOfBlock{$Block}=shift @StartPos;
            $EndPosOfBlock{$Block}=shift @EndPos;
    }
    
    open GFF,"<","$GffFile";
    while (<GFF>) {
            chomp;
            my %bak;
            my ($Chr,$Gene,$StartPos,$EndPos)=(split/\s+/);
            next if $Gene~~@GenesOfRecentWGD;
            foreach my $Block (@{$Chr}) {
                    $bak{$Gene}++ if ($StartPos>=$StartPosOfBlock{$Block} && $EndPos<=$EndPosOfBlock{$Block});
            }
            push @SingleGenes,$Gene if $bak{$Gene}>=1;
    }
    close GFF;
    
    map {$SingleGenes{$_}++} @SingleGenes;
    @SingleGenes=keys %SingleGenes;
    map {say SINGLE $_} @SingleGenes;
    close SINGLE;
    ```
    

```bash
perl Collinear2Retention_Loss.pl -L 0.3 -U 0.8 -K bgy.ks.median.txt -C bgy.collinearity -P bgy -G bgy.gff >STDOUT 2>STDERR

-L,   Lower limit of block ks median
-U,   Upper limit of block ks median
-K,   Ks median file from the former step
-C,   Collinearity file
-P,   Prefix of output file
-G,   Gff file used to run MCScanX
```

---

## 3.3 GO和KEGG富集分析