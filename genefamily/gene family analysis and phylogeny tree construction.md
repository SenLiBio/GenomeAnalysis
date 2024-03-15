# gene family analysis and phylogeny tree construction

## 提前配置

### 软件

|    软件     |     版本     |                             引用                             |
| :---------: | :----------: | :----------------------------------------------------------: |
| orthofinder |    2.3.12    | [Emms, D.M., Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol 20, 238 (2019)](https://doi.org/10.1186/s13059-019-1832-y)<br />[Emms, D.M., Kelly, S. OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biol 16, 157 (2015)](https://doi.org/10.1186/s13059-015-0721-2) |
| jmodeltest2 |    2.1.10    | Darriba D, Taboada GL, Doallo R, Posada D. 2012. jModelTest 2: more models, new heuristics and parallel computing. Nature Methods 9(8), 772.<br />Guindon S and Gascuel O (2003). A simple, fast and accurate method to estimate large phylogenies by maximum-likelihood". Systematic Biology 52: 696-704. |
|   pal2nal   |     v14      | Suyama, M., Torrents, D., & Bork, P. (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. *Nucleic acids research*, *34*(Web Server issue), W609–W612. https://doi.org/10.1093/nar/gkl315 |
|    phyml    | 3.3.20190909 | Guindon, S., & Gascuel, O. (2003). A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood. *Systematic biology*, *52*(5), 696–704. https://doi.org/10.1080/10635150390235520 |
|  raxml-ng   |              | Alexey M. Kozlov, Diego Darriba, Tomáš Flouri, Benoit Morel, and Alexandros Stamatakis (2019) **RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference.** *Bioinformatics, 35 (21), 4453-4455* doi:[10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305) |
|    paml     |     4.9d     | Yang, Z. 2007. PAML 4: a program package for phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24: 1586-1591 |



### Perl 模块安装

需要预先安装Bio::SeqIO模块。可以参考[网上教程](https://www.jianshu.com/p/9e90b3524fe2)，这个还是相对简单的。

---

# 2. 文件格式说明

## 2.1 命名

用于分析的文件命名为：[A-Za-z]{3}.fa

![gene%20famil%2097fb2/Untitled.png](gene%20famil%2097fb2/Untitled.png)

---

## 2.2 序列格式

用于分析的reads的名称是>[a-zA-Z]{3}|.*

![gene%20famil%2097fb2/Untitled%201.png](gene%20famil%2097fb2/Untitled%201.png)

---

---

# 3 蛋白质和cds文件预处理

## 3.1 filter protein and cds by length

```bash
perl filter_fasta_by_length_and_sort.pl LENGTH Prefix.fa Prefix.Passed.fa Prefix.unpassed.fa Prefix.length.txt

# LENGTH: 长度小于这个值的序列会被filter掉
# Prefix.length.txt：每条序列的长度，可以用于检查长度分布
```

## 3.2 add species name to sequence header

```bash
perl add_species_name_to_fa_sequences.pl Prefix.fa
```

---

# 4 获取基因家族信息+画树

## 4.1 cluster gene families by *orthofinder*

### 4.1.1 直接运行

```bash
orthofinder -t 10 -a 5 -S diamond -og -f INPUT_DIR

	-t, Number of parallel sequence search threads [Default = 16]
	-a, Number of parallel analysis threads [Default = 1]
	-S, Sequence search program [Default = blast]. Options: blast, blast_gz, diamond
	-og, Stop after inferring orthogroups
	-f, Input dir containing all protein files
```

### 4.1.2 在原有结果的基础上添加物种

### 4.1.3 在原有结果中删除物种

```bash
orthofinder -b previous_orthofinder_directory -f new_fasta_directory

# This will add each species from the 'new_fasta_directory' to existing set of species, reuse all the previous BLAST results, perform only the new BLAST searches required for the new species and recalculate the orthogroups. The 'previous_orthofinder_directory'  is the OrthoFinder 'WorkingDirectory/' containing the file 'SpeciesIDs.txt'.
```

```bash
orthofinder -b previous_orthofinder_directory

#In the 'WorkingDirectory/' from a previous analysis there is a file called 'SpeciesIDs.txt'. Comment out any species to be removed from the analysis using a '#' character and then run OrthoFinder using the command above. The 'previous_orthofinder_directory'  is the OrthoFinder 'WorkingDirectory/' containing the file 'SpeciesIDs.txt'.
```

### 4.1.4 同时进行添加和删除操作

```bash
orthofinder -b previous_orthofinder_directory -f new_fasta_directory
```

---

## 4.2 construct phylogeny tree

### 4.2.0 预准备

将所有的蛋白质文件cat到一个命名为allpepfiles的文件中。并放到工作目录下

将所有的cds文件放在一个命名为all_cds_dir的文件中

![gene%20famil%2097fb2/Untitled%202.png](gene%20famil%2097fb2/Untitled%202.png)

### 4.2.1 提取单拷贝蛋白质序列

**orthofinder会自动生成这个文件夹，这一步只是保险项。**

Orthogroups.txt来自于*Orthofinder*的输出，allpepfiles是用到的所有物种的蛋白质序列的集合，记得更改脚本里第47行的物种数目。

```bash
perl scripts/1.extract_ortholog.pep.pl ./Orthogroups/Orthogroups.txt allpepfiles Single_Copy_Orthologue_Sequences
```

---

### 4.2.2 *mafft*

注意更改脚本里面对应的单拷贝序列前缀的正则匹配。

```bash
perl scripts/2.create_sh.mafft.pl Single_Copy_Orthologue_Sequences mafft.sh mafftseqs
bash mafft.sh >mafft.log 2>>mafft.err 
mkdir mafftout
mv mafft* mafftout
```

---

---

### 4.2.3 提取蛋白质文件对应的cds序列

记得修改第52行的物种数目

```bash
perl scripts/3.extract_ortholog.cds.pl ./Orthogroups/Orthogroups.txt all_cds_dir/ singlecopycdsdir length_check.txt
```

---

### 4.2.4 cds结果排序

```bash
perl scripts/4.sort_cds_mafft.pl mafftout/mafftseqs singlecopycdsdir sortedcdsdir
```

---

### 4.2.5 *pal2nal*

注意更改脚本中的*pal2nal*软件的路径

```bash
perl scripts/5.create_pal2aln.pl mafftout/mafftseqs sortedcdsdir pal2nalfas pal2nal.sh
nohup bash pal2nal.sh >pal2nal.log 2>>pal2nal.err &
mkdir pal2naloutdir
mv pal2nal* err* pal2naloutdir

# 如果存在protein和cds不对应的情况，应当删除这些
```

---

### 4.2.6 *Gblocks*筛选保守位点（这一步可以选做）

如果不能直接调用gblocks，可以修改line10中指定的路径

```bash
perl scripts/6.create.gblocks_sh.pl pal2naloutdir/pal2nalfas gblockshtms gblocks.sh
bash gblocks.sh
mkdir gblocksoutdir
mkdir gblocksgbs
mv pal2naloutdir/pal2nalfas/*htm gblockshtms
mv pal2naloutdir/pal2nalfas/*gb gblocksgbs
mv gblocks* gblocksoutdir
```

---

### 4.2.7 merge cds

```bash
# 不做gblocks的情况
perl scripts/7.cat_cds.pl pal2naloutdir/pal2nalfas mergedcdss

# 做gblocks的情况
perl scripts/7.cat_cds.pl gblocksoutdir/gblocksgbs mergedcdss
```

---

### 4.2.8 生成phylip文件

这里有个很有趣的事情，就是有时候`convertFasta2Phylip.sh`的换行符会出现变化，然后就无法使用，所以建议直接用包里面的脚本。

```bash
mkdir finalseqs
cat mergedcdss/*fa >finalseqs/final.fa
cd finalseqs
sh ../scripts/convertFasta2Phylip.sh final.fa >final.phy
```

---

### 4.2.9 用*jmodeltest2*寻找最佳模型

```bash
mkdir tree
cd tree
ln -s ../finalseqs/final.phy .
nohup java -jar -XX:ParallelGCThreads=4 -Xmx4g ~/software/jmodeltest-2.1.10/jModelTest.jar -tr 10 -d final.phy -s 11 -f -i -g 8 -AIC -BIC -AICc -DT -p -a -w -o bestmodel &

	-tr,	number of threads to execute (default is 32)
	-d,		input data file
	-s,		number of substitution schemes (e.g., -s 11) (it has to be 3,5,7,11,203; 				default is 3)
	-f,		Include models with unequals base frecuencies
	-i,		include models with a proportion invariable sites (e.g., -i) (default is 				false)
	-g,		include models with rate variation among sites and number of categories 				(e.g., -g 8) (default is false & 4 categories)
	-AIC,	calculate the Akaike Information Criterion (e.g., -AIC) (default is 					false)
	-AICc,	calculate the corrected Akaike Information Criterion (e.g., -AICc) 						(default is false)
	-BIC,	calculate the Bayesian Information Criterion (e.g., -BIC) (default is 					false)
	-DT,	calculate the decision theory criterion (e.g., -DT) (default is false)
	-p,		calculate parameter importances (e.g., -p) (default is false)
	-a,		estimate model-averaged phylogeny for each active criterion (e.g., -a) 					(default is false)
	-w,		write PAUP block (e.g., -w) (default is false)
	-o,		set output file (e.g., -o jmodeltest.out)
```

---

### 4.2.10 *phyml*建树

其他参数是根据*jmodeltest*来的，然后-b是bootstrap值

```bash
phyml -i final.phy -d nt -n 1 -b 1000 --run_id 012345 -m GTR+I+G -f m -v e -c 4 -a e --no_memory_check -o tlr -s BEST

	-i,		seq_file_name
			seq_file_name is the name of the nucleotide or amino-acid sequence file in 				PHYLIP format.
	-d,		data_type
			data_type is 'nt' for nucleotide (default), 'aa' for amino-acid	sequences, 				or 'generic'
	-n,		nb_data_sets
			nb_data_sets is an integer corresponding to the number of data sets to  			 analyse.
	-b,		int
			int >  0: int is the number of bootstrap replicates.
			int =  0: neither approximate likelihood ratio test nor bootstrap values 				are computed.
			int = -1: approximate likelihood ratio test returning aLRT statistics.
			int = -2: approximate likelihood ratio test returning Chi2-based 						parametric branch supports.
			int = -4: SH-like branch supports alone.
			int = -5: (default) approximate Bayes branch supports.
	-m,		model
	-f,		e, m, or fA,fC,fG,fT (具体参照phyml -help)
	-v,		prop_invar
			prop_invar : proportion of invariable sites. Can be a fixed value in the 				[0,1] range or ‘e’ to get the maximum likelihood estimate.
	-c,		nb_subst_cat
			nb_subst_cat : number of relative substitution rate categories. Default : 				nb_subst_cat=4.	Must be a positive integer.
	-a,		gamma
			gamma : distribution of the gamma distribution shape parameter. Can be a 				fixed positive value or ‘e’ to get the maximum likelihood estimate.
	-o,		params
			This option focuses on specific parameter optimisation.
			t: tree topology are optimised; l: branch length are optimised;
			r: rate parameters are optimised; n: no parameter is optimised.
	-s,		move
			Tree topology search operation option. Can be either NNI (default, fast) 				or SPR (a bit slower than NNI) or BEST (best of NNI and SPR search).
```

---

### 4.2.11 *raxml-ng* 建树

```bash
raxml-ng --all -msa final.phy --model GTR+I+G --prefix final --seed 2 --threads 10 --bs-metric fbp,tbe --bs-trees 1000

		--all,			all-in-one (ML search + bootstrapping)
			-msa,			alignment file
		--model,		model specification OR partition file
		--prefix,		prefix for output files (default: MSA file name)	
		--seed,			seed for pseudo-random number generator (default: current time)
	--threads,		number of parallel threads to use (default: 24)
	--bs-metric,	branch support metric: fbp = Felsenstein bootstrap (default), tbe 						= transfer distance	
	  --bs-trees,	number of bootstraps replicates.
```

---

### 4.2.12 *iqtree*建树

```bash
iqtree -s final.phy -st CODON -m GTR+I+G -b 1000 -nt 10

	-s,		Input alignment in PHYLIP/FASTA/NEXUS/CLUSTAL/MSF format
	-st,	BIN, DNA, AA, NT2AA, CODON, MORPH (default: auto-detect)
	-m,		model specification. if '-m MFP -mtree' used, this software will 						automatically search for the best model followed by tree inference
	-bb,	Ultrafast bootstrap (>=1000) (这个是iqtree特有的快速bootstrap)
	-b,		Bootstrap + ML tree + consensus tree (>=100)
	-nt,	number of parallel threads to use
```

---

### 4.2.13 mcmctree计算分歧时间

- mcmctree.tree
  
    ```
    10 1
    (Nnu,((Tha,Ath),((Clo,Bgy)'L(0.478)',((Mes,Rco),(Spu,(Peu,Ptr))))'L(0.48)'))'B(1.196,1.2863)';
    ```
    
- mcmctree.ctl
  
    ```bash
    seed = -1
    seqfile = mcmctree.phy
    treefile = mcmctree.tree
    outfile = mcmctree.out
    ndata = 3
    usedata = 1 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
    clock = 3 * 1: global clock; 2: independent; and 3: correlated rates
    RootAge = '<1.2863' * safe constraint on root age, used if no fossil for root.
    model = 4 * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
    alpha = 0 * alpha for gamma rates at sites
    ncatG = 5 * No. categories in discrete gamma
    cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
    BDparas = 1 1 0 * birth, death, sampling
    kappa_gamma = 6 2 * gamma prior for kappa
    alpha_gamma = 1 1 * gamma prior for alpha
    rgene_gamma = 2 2 * gamma prior for rate for genes
    sigma2_gamma = 1 10 * gamma prior for sigma^2 (for clock=2 or 3)
    finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, etc.
    print = 1
    burnin = 200000
    sampfreq = 10
    nsample = 50000
    ```
    

```bash
perl scripts/phy2mcmctreephy.pl final.phy
mcmctree mcmctree.ctl
```

---

# 5. CAFE做基因家族演化分析

## 5.1 软件和自带脚本包的下载和安装

[^note]: 这里用的是CAFE4.2，下载Python包的网站有时候不可用，直接从我的github上面下载比较好。

cafe的[github主页](https://github.com/hahnlab/CAFE)、[最新版下载页面](https://github.com/hahnlab/CAFE/releases/tag/v4.2.1)、[官网](https://hahnlab.github.io/CAFE/download.html)、[自带的Python包](https://iu.app.box.com/v/cafetutorial-files)。

```bash
# 解压和安装
./configure
make

# cafe的执行文件在“release”文件夹中
```

## 5.2 输入文件和config文件准备

```bash
# "Orthogroups.GeneCount.tsv"文件在Orthofinder输出中的“Orthogroups”文件夹。

perl GeneCount2CafeInput.pl ../Orthogroups/Orthogroups.GeneCount.tsv >cafeinput.raw
python cafetutorial_clade_and_size_filter.py -s -i cafeinput.raw -o cafeinput.filtered
```

## 5.3 运行CAFE

```BASH
cafe cafe_script.sh

# 下面是cafe_script.sh文件
#! ~/software/CAFE/release/cafe
date
load -i cafeinput.filtered -p 0.01 -t 10 -l log.txt
# -t, 使用的线程数
# -p, 显著性判断阈值
#the phylogenetic tree structure with branch lengths
# 这棵树可以直接用mcmctree输出的树
tree (Nnu:124.5781,((Ath:24.2897,Tha:24.2897):94.0914,(((((Bse:13.0498,Bgy:13.0498):2.2495,Bcy:15.2993):22.3403,Rap:37.6396):15.1566,Clo:52.7963):40.9378,((Mes:67.4141,Rco:67.4141):19.6264,(Spu:22.7468,(Ptr:10.8282,Peu:10.8282):11.9186):64.2937):6.6936):24.6470):6.1971)
#search for 2 parameter model
# 也可以对不同的clade设置不同的lambda值，最终的结果好像差别不是很大
lambda -s -t (1,((1,1)1,(((((1,1)1,1)1,1)1,1)1,((1,1)1,(1,(1,1)1)1)1)1)1)
report reportout1
date
```

## 5.4 结果统计

```bash
python report_analysis.py -i reportout1.cafe -r 1 -o rapid

# -r, 选择1就是统计p值显著的基因家族，0就是所有基因家族
```

## 5.5 用拟南芥的基因做基因家族功能注释

```bash
perl Rapid_GF_Anno.pl rapid_fams.txt Orthogroups.txt ath.gene.discription
```



# 6 脚本和整合的shell

---

## 6.1 整个流程的shell

注意：

1. 脚本放在**`./scripts`**目录下
2. 记得更改`1.extract_ortholog.pep.pl`和`3.extract_ortholog.cds.pl`中的物种数目（分别在第47行和第52行）
3. 记得检查`5.create_pal2aln.pl`中pal2nal的路径（第11行）
4. 记得检查`6.create.gblocks_sh.pl`中gblocks的路径（第9行）
5. 记得更改最后的jModelTest.jar的路径
6. 由于用Windows或者某些软件打开后，换行符那里会添加一个`^M`,linux识别的时候可能会出现问题，建议先用`cat -A final.sh` 检查换行符，然后用`sed -e 's/^M//g' final.sh >1.bak | mv 1.bak final.sh` 来去除^M，注意打^M的方式是`Ctrl+M`而不是分别打出^和M。

```bash
# 1. 脚本放在**`./scripts`**目录下
# 2. 记得更改`1.extract_ortholog.pep.pl`和`3.extract_ortholog.cds.pl`中的物种数目（分别在第47行和第52行）
# 3. 记得检查`5.create_pal2aln.pl`中pal2nal的路径（第11行）
# 4. 记得检查`6.create.gblocks_sh.pl`中gblocks的路径（第9行）
# 5. 记得更改最后的jModelTest.jar的路径
# 6. 由于用Windows或者某些软件打开后，换行符那里会添加一个`^M`,linux识别的时候可能会出现问题，建议先用`cat -A [final.sh](http://final.sh)` 检查换行符，然后用`sed -e 's/^M//g' final.sh >1.bak | mv 1.bak final.sh` 来去除^M，注意打^M的方式是`Ctrl+M`而不是分别打出^和M。
# 7. 由于下载到的序列可能会有蛋白质和cds序列不匹配的情况，所以pal2nal那一步可能会有报错，导致无法往后进行，建议跑shell跑到这一步，然后手动删除文件夹中的空文件，之后手动跑后面的步骤。

# perl scripts/1.extract_ortholog.pep.pl ./Orthogroups/Orthogroups.txt allpepfiles Single_Copy_Orthologue_Sequences &&
perl scripts/1.extract_ortholog.pep.pl ./Orthogroups/Orthogroups.txt allpepfiles Single_Copy_Orthologue_Sequences 
perl scripts/2.create_sh.mafft.pl Single_Copy_Orthologue_Sequences mafft.sh mafftseqs
bash mafft.sh >mafft.log 2>>mafft.err 
mkdir mafftout
mv mafft* mafftout
perl scripts/3.extract_ortholog.cds.pl ./Orthogroups/Orthogroups.txt all_cds_dir/ singlecopycdsdir length_check.txt 
perl scripts/4.sort_cds_mafft.pl mafftout/mafftseqs singlecopycdsdir sortedcdsdir 
perl scripts/5.create_pal2aln.pl mafftout/mafftseqs sortedcdsdir pal2nalfas pal2nal.sh 
bash pal2nal.sh >pal2nal.log 2>>pal2nal.err 
mkdir pal2naloutdir 
mv pal2nal* err* pal2naloutdir 
perl scripts/6.create.gblocks_sh.pl pal2naloutdir/pal2nalfas gblockshtms gblocks.sh &&
bash gblocks.sh >gblocks.log 2>>gblocks.err
mkdir gblocksoutdir
mkdir gblocksgbs
mv pal2naloutdir/pal2nalfas/*htm gblockshtms
mv pal2naloutdir/pal2nalfas/*gb gblocksgbs
mv gblocks* gblocksoutdir
perl scripts/7.cat_cds.pl gblocksoutdir/gblocksgbs mergedcdss
mkdir finalseqs
cat mergedcdss/*fa >finalseqs/final.fa
cd finalseqs
sh ../scripts/convertFasta2Phylip.sh final.fa >final.phy
mkdir tree
cd tree
ln -s ../finalseqs/final.phy .
nohup java -jar -XX:ParallelGCThreads=4 -Xmx4g ~/software/jmodeltest-2.1.10/jModelTest.jar -tr 10 -d final.phy -s 11 -f -i -g 8 -AIC -BIC -AICc -DT -p -a -w -o bestmodel &
```

---

## 6.2 脚本

### add_species_name_to_fa_sequences.pl

```bash
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
```

---

### filter_fasta_by_length_and_sort.pl

```bash
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
```

---

### 1.extract_ortholog.pep.pl

```bash
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
```

---

### 2.create_sh.mafft.pl

```bash
use warnings;
use strict;
my $homolog_dir=shift;
my @files=<$homolog_dir/*>;
my $output_sh=shift;
open (O,">$output_sh");

my $outdir=shift;
system("mkdir $outdir");
foreach my $file (@files){
    $file=~/OG(\d+)/;
    my $sign=$1;
    my $outfile="$outdir/mafft$1";
    print O "mafft --maxiterate 1000 --localpair $file >$outfile\n";
}
close O;
```

---

### 3.extract_ortholog.cds.pl

```bash
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
```

---

### 4.sort_cds_mafft.pl

```bash
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
    my $cdsfile="$cdsdir/OG$sign";

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
```

---

### 5.create_pal2aln.pl

```bash
use warnings;
use strict;
my $pepdir=shift;
my @mafftfiles=<$pepdir/*>;
my $cdsdir=shift;
my $outdir=shift;
system ("mkdir $outdir");
my $output_sh=shift;
open(O,">$output_sh");

my $pal2alnpath="~/software/pal2nal.v14/pal2nal.pl";
foreach my $pepfile (@mafftfiles){
    $pepfile=~/mafft(\d+)/;
    my $sign=$1;
    print O "perl $pal2alnpath $pepfile $cdsdir/paltoaln$sign -output fasta >$outdir/pal2aln_$sign\n";
}
close O;
```

---

### 6.create.gblocks_sh.pl

```bash
use strict;
my $pepdir=shift;
my @files=<$pepdir/*>;
my $logdir=shift;
system ("mkdir $logdir");
my $output_sh=shift;
open(O,">$output_sh");

my $gblockspath="~/software/Gblocks_0.91b/Gblocks";
foreach my $pepfile (@files){
    print O "$gblockspath $pepfile -t=c\n";
}
close O;
```

---

### 7.cat_cds.pl

```bash
use warnings;
use strict;
use Bio::SeqIO;
my $trimmeddir=shift;
my @pepfiles=<$trimmeddir/*>;
my $pepdir=shift;
system("mkdir $pepdir");

my %pep;
my $count=0;
foreach my $file (sort @pepfiles){
    my $fa=Bio::SeqIO->new(-file=>$file,-format=>'fasta');
    while (my $seq_obj=$fa->next_seq){
        my $pepid=$seq_obj->display_name;
        my $pepseq=$seq_obj->seq;
    $pepid=~/(\w+)\|(.+)/;
    my $species=$1;
    $pep{$species}.=$pepseq;
    }
    $count++;
}
print "total files:$count\n";

my @species=sort keys %pep;
foreach my $species (@species){
    my $output="$pepdir/$species.fa";
    open (O,">$output");
    print O ">$species\n$pep{$species}\n";
    close O;
}
```

---

### convertFasta2Phylip.sh

```bash
#! /bin/sh

if [ $# != 1 ]; then
    echo "USAGE: ./script <fasta-file>"
    exit
fi

numSpec=$(grep -c  ">" $1)
tmp=$(cat $1 | sed "s/>[ ]*\(\w*\).*/;\1</"  | tr -d "\n" | tr -d ' '  | sed 's/^;//' | tr "<" " " )
length=$(($(echo $tmp | sed 's/[^ ]* \([^;]*\);.*/\1/'   | wc -m ) - 1))

echo "$numSpec $length"
echo  $tmp | tr ";" "\n"
```

---

### phy2mcmctreephy.pl

```perl
open PHY,"<$ARGV[0]";
open OUT,">","mcmctree.phy";

@lines=(<PHY>);
$header=shift @lines;
my ($species,$base)=split/\s+/,$header;
$codonbase=$base/3;

for (1..$species) {
        my $seq=shift @lines;
        my ($spe,$codon)=split/\s+/,$seq;
        my @bases=split//,$codon;
        my $count=1;
        foreach my $base (@bases) {
                $codon1{$spe}.="$base" if $count%3==1;
                $codon2{$spe}.="$base" if $count%3==2;
                $codon3{$spe}.="$base" if $count%3==0;
                $count++;
        }
}

foreach my $num (1..3) {
        print OUT "$species $codonbase\n";
        my $hashname="codon" . "$num";
        foreach my $spe (sort keys %codon1) {
                print OUT "$spe  ${$hashname}{$spe}\n";
        }
}

```

---

### mcmctree.ctl

```bash
seed = -1
seqfile = mcmctree.phy
treefile = mcmctree.tree
outfile = mcmctree.out
ndata = 3
usedata = 1 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
clock = 3 * 1: global clock; 2: independent; and 3: correlated rates
RootAge = '<1.2863' * safe constraint on root age, used if no fossil for root.
model = 4 * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 0 * alpha for gamma rates at sites
ncatG = 5 * No. categories in discrete gamma
cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
BDparas = 1 1 0 * birth, death, sampling
kappa_gamma = 6 2 * gamma prior for kappa
alpha_gamma = 1 1 * gamma prior for alpha
rgene_gamma = 2 2 * gamma prior for rate for genes
sigma2_gamma = 1 10 * gamma prior for sigma^2 (for clock=2 or 3)
finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, etc.
print = 1
burnin = 200000
sampfreq = 10
nsample = 50000
```

---

### GeneCount2CafeInput.pl

```perl
use 5.011;

while (<>) {
        if (/^Orthogroup/) {
                chomp;
                my $line=$_;
                my @eles=split/\s+/,$line;
                print "Desc\tFamily ID\t";
                shift @eles;
                pop @eles;
                my $lastone=pop @eles;
                map {print "$_\t"} @eles;
                say $lastone;
        } else {
                chomp;
                my $line=$_;
                my @eles=split/\s+/,$line;
                print "(null)\t";
                pop @eles;
                my $lastone=pop @eles;
                map {print "$_\t"} @eles;
                say $lastone;
        }
}
```

---

### clade_and_size_filter.py

```python
"""
This script defines and runs the clade_filter and size_filter functions. The former keeps only gene families with gene copies in at least two species of the specified clades; this is necessary because ancestral state reconstruction requires at least two species per clade of interest. The latter separates gene families with < 100 from > 100 gene copies in one or more species; this is necessary because big gene families cause the variance of gene copy number to be too large and lead to noninformative parameter estimates.
The output should be two CAFE input files, one for gene families with < 100 gene copies in all species, another for the remaining gene families. The first file should be used to estimate parameter values, and these values should then be used to analyse the second file.
"""

__author__ = "Fabio H. K. Mendes"

import os
import argparse

def clade_filter(mcl_dump, clade_str):
    """
    Return set of lines to print after checking if at least 2 species in specified clades have gene copies for a given gene family

    [str] mcl_dump: path to mcl's dump file
    [str] clade_str: clades of interest (separated by white spaces, with species within clades separated by comma)
                       (e.g., "ENSBTA,ENSCFA,ENSECA  ENSMUS,ENSNLE,ENSPTR,ENSPAN,ENSPPY,ENSCJA,ENSP00 ENSMMU,ENSRNO")
    """
    lines_to_keep_list = list()
    spp_idx_dict = dict()
    clades_list = list()
    if clade_str: # if clade filter was specified
        clades_list = clade_str.split(" ")

    with open(mcl_dump, "r") as input_file:
        for line_n, line in enumerate(input_file):
            line = line.rstrip()
            tokens = line.split("\t")
            spp_info = tokens[2:]

            if line.startswith("Desc"):
                spp_idx_dict = dict((sp, idx) for idx,sp in enumerate(spp_info))
                continue

            if clades_list:
                clades_ok_list = list()

                for clade in clades_list:
                    spp_list = clade.split(",")
                    clade_count = sum(1 for sp in spp_list if int(spp_info[spp_idx_dict[sp]]) >= 1)

                    if clade_count >= 2:
                        clades_ok_list.append(1)

                if sum(clades_ok_list) == len(clades_list):
                    lines_to_keep_list.append(line_n)

            # just keeping lines where >=2 species (among all of them) have gene copies
            clade_count = sum(1 for sp_count in spp_info if int(sp_count) >= 1)
            if clade_count >= 2:
                lines_to_keep_list.append(line_n)

    return set(lines_to_keep_list)

def size_filter(mcl_dump, lines_to_keep_set):
    """
    Return set of lines to print after checking if at least 2 species in specified clades have gene copies for a given gene family

    [str] mcl_dump: path to mcl's dump file
    """
    lines_to_remove_set = set()
    size_cutoff = 100
    fam_size = int()
    with open(mcl_dump, "r") as input_file:
        for line_n, line in enumerate(input_file):
            line = line.rstrip()

            if line.startswith("Desc"):
                continue

            elif line_n not in lines_to_keep_set and len(lines_to_keep_set) > 0:
                continue

            tokens = line.split("\t")
            spp_info = tokens[2:]

            for gene_count in spp_info:
                if int(gene_count) >= size_cutoff:
                    lines_to_separate_set.add(line_n)

    lines_to_keep_set -= lines_to_separate_set

    return lines_to_keep_set, lines_to_separate_set

def filter_print(mcl_dump, lines_to_keep_set, lines_to_separate_set, output_file_name):
    """
    Print two mcl input files, one with gene families having < 100 gene copies, the other with gene families having > 100 copies
    """
    if len(lines_to_keep_set) == 0 and len(lines_to_separate_set) == 0:
        exit("No filtering was done! Exiting...\n")

    with open(output_file_name, "w") as output_file:
        with open("large_"+output_file_name, "w") as output_file2:
            with open(mcl_dump, "r") as input_file:
                for line_n, line in enumerate(input_file):
                    line = line.rstrip() + "\n"

                    if line_n == 0:
                        output_file.write(line)
                        output_file2.write(line)

                    elif line_n in lines_to_keep_set and len(lines_to_keep_set) >= 1:
                        output_file.write(line)

                    elif line_n not in lines_to_separate_set and len(lines_to_keep_set) == 0:
                        output_file.write(line)

                    # has to be if, not elif
                    if line_n in lines_to_separate_set:
                        output_file2.write(line)

        # cleaning up in case size filtering was not done
        if len(lines_to_separate_set) == 0:
            os.unlink("large_"+output_file_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, prog="cafetutorial_clade_and_size_filter.py")
    parser.add_argument("-i", "--input-file", action="store", dest="input_file", required=True, type=str, help="full path to mcl's output dump file")
    parser.add_argument("-o", "--output-file", action="store", dest="output_file", required=True, type=str, help="full path to file to be written")
    parser.add_argument("-cl", "--clade-filter", action="store", dest="clade_str", default=None, required=False, type=str, help="list of clades (separated by white spaces) comprised of species identifiers (separated by comma) that must have at least two species with gene copies for a given gene family")
    parser.add_argument("-s", "--size-filter", action="store_true", dest="size_filter", required=False, help="option to perform size filtering")

    args = parser.parse_args()

    lines_to_keep_set, lines_to_separate_set = set(), set()

    # applying size filter (if no groups are specified, just the lines where just 1 species has gene counts are removed
    lines_to_keep_set = clade_filter(args.input_file, args.clade_str)

    if args.size_filter:
        lines_to_keep_set, lines_to_separate = size_filter(args.input_file, lines_to_keep_set)

    filter_print(args.input_file, lines_to_keep_set, lines_to_separate_set, args.output_file) # .add(0) to get header back
```

---

### cafe_script.sh

```
#! ~/software/CAFE/release/cafe
date
load -i cafeinput.filtered -p 0.01 -t 10 -l log.txt
# -t, 使用的线程数
# -p, 显著性判断阈值
#the phylogenetic tree structure with branch lengths
# 这棵树可以直接用mcmctree输出的树
tree (Nnu:124.5781,((Ath:24.2897,Tha:24.2897):94.0914,(((((Bse:13.0498,Bgy:13.0498):2.2495,Bcy:15.2993):22.3403,Rap:37.6396):15.1566,Clo:52.7963):40.9378,((Mes:67.4141,Rco:67.4141):19.6264,(Spu:22.7468,(Ptr:10.8282,Peu:10.8282):11.9186):64.2937):6.6936):24.6470):6.1971)
#search for 2 parameter model
# 也可以对不同的clade设置不同的lambda值，最终的结果好像差别不是很大
lambda -s -t (1,((1,1)1,(((((1,1)1,1)1,1)1,1)1,((1,1)1,(1,(1,1)1)1)1)1)1)
report reportout1
date
```

---

### Rapid_GF_Anno.pl

```perl
use strict;
use autodie;

die "Usage: perl $0 rapid_fams.txt Orthogroups.txt ath.gene.discription\n" if @ARGV != 3 ;

open RGF,"<","$ARGV[0]";
open ORTHO,"<","$ARGV[1]";
open ANNO,"<","$ARGV[2]";

my %Ath_Genes;
while (<ORTHO>) {
		s/\s+$//;
		my @eles=split/\s+/,$_;
		my $OG=shift @eles;
		$OG=~s/://;
		my @Ath_Genes;
		foreach my $ele (@eles) {
			push @Ath_Genes,$ele if $ele=~/Ath/;
		}
		$Ath_Genes{$OG}=[@Ath_Genes];
		undef @Ath_Genes;
}
close ORTHO;


my %Description;
while (<ANNO>) {
	next if /^Locus/;
	next if /^\s+$/;
	s/\s+$//;
	s/\(source\:Araport11\)//;
	s/\s+protein_coding//;
	my @eles=split/\s+/;
	my $Locus=shift @eles;
	my $Gene=shift @eles;
	my $Description=join(' ',@eles);
	$Description{$Gene}=$Description;
}

while (<RGF>) {
	next if (/^#/ || /^Overall/);
	s/\s+$//;
	my @eles=split/\s+/;
	my $spe=shift @eles;
	$spe=~s/\<//;
	$spe=~s/\>//;
	my @OGs=split/,/,(shift @eles);
	(my $out=$spe)=~s/:/\.rapid\.AthAnno/;
	open OUT,">","$out";
	foreach my $OG (@OGs) {
		$OG=~m/(\+|\-)/;
		my $Status=$1;
		$OG=~s/\[.*\]//;
		my @Ath_Genes=@{$Ath_Genes{$OG}};
		foreach my $Ath_Gene (@Ath_Genes) {
			$Ath_Gene=~s/^Ath\|//;
			print OUT "$OG$Status: $Description{$Ath_Gene}\n";
		}
	}

}
```

---

### report_analysis.py

```python
#!/usr/bin/python
########################################################################################
# The new CAFE Report Analysis script.
#
# Gregg Thomas, Spring 2016
########################################################################################

import sys, os, argparse, cafecore as cafecore

############################################
#Function Definitions
############################################

def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Analyzes a CAFE report file (.cafe)");

	parser.add_argument("-i", dest="report_input_file", help="A CAFE report file (.cafe).");
	parser.add_argument("-r", dest="rapids_out", help="1: Output file will be list of only rapidly changing families on each node. 0: Output file will be list of all changing families on each node. Default: 1", type=int, default=1);
	parser.add_argument("-l", dest="large_fam_file", help="A CAFE report file for the large families of the same set of species.");
	parser.add_argument("-o", dest="output_prefix", help="A prefix string to put on all the output files generated.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.report_input_file == None or args.output_prefix == None:
			cafecore.errorOut(1, "A CAFE report file must be specified with -i and an output file name with -o");
			optParse(1);

		if args.rapids_out not in [1,0]:
			cafecore.errorOut(2, "-r can only take values of 0 or 1");
			optParse(1);

		return args.report_input_file, args.rapids_out, args.large_fam_file, args.output_prefix;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

#######################

def formatLineParse(line):
# This function handles CAFE's weird node pair format for the p-values and node ids.

	if "=" in line:
		line = line.split("=")[1];
	if ":" in line:
		line = line.split(": ")[1].strip();
	line = line.replace("(", "").replace(")", "");
	line = line.split(" ");
	line = [f.split(",") for f in line];

	return line;

#######################

def nodeRelabel(treedict):
# Family trees are read with gene counts on the tip labels. This function removes them.

	tmp = {};
	#print treedict;

	for oldkey in treedict:
		if treedict[oldkey][3] == 'tip':
			newkey = oldkey[:oldkey.index("_")];
			tmp[newkey] = treedict[oldkey];
		else:
			tmp[oldkey] = treedict[oldkey];

	return tmp;

#######################

def nodeMap(cafetd, mytd):
# CAFE has pre-determined mappings in the tree. When I read the tree with my own script the mappings
# are different. This function creates a map from my node ids to CAFE's node ids.
# The dictionary nodemap has the following {key:value} format: {my node id:CAFE's node id}

	nodemap = {};
	# The map dictionary.

	##############
	# for node in cafetd:
	# 	if cafetd[node][3] == 'tip':
	# 		spec = node[:node.index("<")];
	# 		cafeid = node[node.index("<")+1:node.index(">")];
	# 		nodemap[cafeid] = spec;

	# while len(nodemap) != len(cafetd):
	# 	for node in cafetd:
	# 		if cafetd[node][3] == 'root':
	# 			continue;

	# 		orignode = node;
	# 		node = node[node.index("<")+1:node.index(">")];

	# 		if node in nodemap:
	# 			if cafetd[orignode][3] == 'tip':
	# 				curanc = cafetd[orignode][1];
	# 				mapanc = mytd[nodemap[node]][1];
	# 			else:
	# 				curanc = cafetd[orignode][1];
	# 				mapanc = mytd["<" + nodemap[node] + ">"][1];

	# 			nodemap[curanc.replace("<","").replace(">","")] = mapanc.replace("<","").replace(">","");
	##############
	# The above formats nodemap with the reverse {key:value} format: {CAFE's node id:my node id}

	for node in cafetd:
		if cafetd[node][3] == 'tip':
			spec = node[:node.index("<")];
			cafeid = node[node.index("<"):];
			nodemap[spec] = node;
	# First map the tips by their unique species labels.

	while len(nodemap) != len(mytd):
		for node in mytd:
			if mytd[node][3] == 'root':
				continue;

			if node in nodemap:
				curanc = mytd[node][1];
				mapanc = cafetd[nodemap[node]][1];

				nodemap[curanc] = mapanc;
	# Then do a post-order traversal and map the current node's ancestor to it's map's ancestor.

	return nodemap;

#######################
def cra(inlines, results, node_fams, linestart, afilename, s_nodes, v):
	numbars = 0;
	donepercent = [];
	i = 0;
	acount = 0;
	afile = open(afilename, "a");

	for inline in inlines:
	# Each line of the report file is read.
		if v == 1:
			numbars, donepercent = cafecore.loadingBar(i, len(inlines), donepercent, numbars);
		i = i + 1;

		if i <= linestart:
			continue;
		# If the line is a CAFE info line, skip it.

		inline = inline.strip().split("\t");
		famid = inline[0];
		famtree = inline[1];
		nodeformat = inline[3].replace("),(", ") (");
		# Parsing the information for the current family.

		outline = famid + "\t";
		outlist = [0 for n in s_nodes];
		# Prep for the anc states file.

		nodes = formatLineParse(nodeformat);

		tlinfo, newfamtree = cafecore.treeParseNew(famtree, 1);
		# Reading the tree and adding my node labels.

		for tlnode in tlinfo:
			if tlinfo[tlnode][3] == 'root':
				tlinfo[tlnode].append(famtree[famtree.rfind("_")+1:]);
			elif tlinfo[tlnode][3] == 'tip':
				tlinfo[tlnode].append(tlnode[tlnode.index("_")+1:]);
			else:
				tlinfo[tlnode][4] = tlinfo[tlnode][4][1:];
		# Gene counts for each node are read as support values for internal nodes, but must
		# have the underscore removed. Tip and root node counts are added here as well.

		tlinfo = nodeRelabel(tlinfo);
		# Removes the gene counts from the tip node labels.

		if i == (linestart + 1):
			maps = nodeMap(tinfo, tlinfo);
		# If this is the first family, we need to build our maps from my node ids to CAFE's.

		for tlnode in tlinfo:
		# For each node in the current gene family tree, we make our counts.

			if tlinfo[tlnode][3] == 'root':
				continue;
			# No count is made at the root of the tree.

			curanc = tlinfo[tlnode][1];
			curmap = maps[tlnode];
			# Get the ancestor and the map of the current node.

			curcount = int(tlinfo[tlnode][4]);
			anccount = int(tlinfo[curanc][4]);
			# Get the gene counts of the current node and the ancestor.

			outlist[s_nodes.index(curmap)] = str(curcount);
			# Save the count of the current node to be sent to the anc states file.

			diff = curcount - anccount;
			# Calculate the difference in gene count.

			typeflag = 0;
			# typeflag tells whether an expansion or contraction has occurred.

			if curcount > anccount:
				typeflag = 1;
				results[curmap][0] += 1;
				results[curmap][1] += diff;

				if r_opt == 0:
					node_fams[curmap][0].append(famid + "[+" + str(diff) + "]");

					if famid not in node_fams['total']:
						node_fams['total'].append(famid);
			# If the difference in gene count between the current node and the ancestor is positive, an
			# expansion has occurred. This makes the appropriate counts.

			elif curcount < anccount:
				typeflag = 2
				results[curmap][3] += 1;
				results[curmap][4] += abs(diff);

				if curcount == 0 and anccount != 0:
					results[curmap][5] += 1;

				if r_opt == 0:
					node_fams[curmap][1].append(famid + "[" + str(diff) + "]");

					if famid not in node_fams['total']:
						node_fams['total'].append(famid);
			# If the difference in gene count between the current node and the ancestor is negative, a
			# contraction has occurred. This makes the appropriate counts. It also checks for family losses
			# along that branch by seeing if the current node has transitioned to a count of 0.

			elif curcount == anccount:
				results[curmap][2] += 1;
			# Otherwise, the counts at the current node and the ancestor are the same and no change has occurred.

			if float(inline[2]) < 0.01:
			# If the family p-value is below a threshold, the family is rapidly evolving.

				if r_opt == 1:
					if famid not in node_fams['total']:
						node_fams['total'].append(famid);
				# Add the family id to the 'total' key of node_fams. This also parses the nodeformat line which is
				# in the paired CAFE node format.

				pairnodeid = curmap[curmap.index("<")+1:curmap.index(">")];
				# Since the paired format does not include the brackets that the other node labels do, I have to
				# remove them to check against that format.

				for j in range(len(nodes)):
					for k in range(len(nodes[j])):
						if formatline[j][k] == pairnodeid and float(nodes[j][k]) < 0.01:
							if typeflag == 1:
								results[curmap][7] += 1;
								if r_opt == 1:
									node_fams[curmap][0].append(famid + "[+" + str(diff) + "*]");
								elif r_opt == 0:
									node_fams[curmap][0].pop();
									node_fams[curmap][0].append(famid + "[+" + str(diff) + "*]");
							elif typeflag == 2:
								results[curmap][8] += 1;
								if r_opt == 1:
									node_fams[curmap][1].append(famid + "[" + str(diff) + "*]");
								elif r_opt == 0:
									node_fams[curmap][1].pop();
									node_fams[curmap][1].append(famid + "[" + str(diff) + "*]");
							results[curmap][9] += 1;
				# Runs through the paired format as a list of lists. If the p-value of that node is less than a threshold
				# that branch is rapidly evolving. Based on typeflag, the appropriate counts are made. The family id is
				# also added to the current node in rapids.

		outline += "\t".join(outlist) + "\n";
		afile.write(outline);
		# Write the states of the current family to the anc states file

	if v == 1:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	afile.close();
	return results, node_fams;

############################################
#Main block
############################################

infilename, r_opt, largefilename, outprefix = optParse(0);
#Get the input parameters.

famfilename = outprefix + "_fams.txt";
nodefilename = outprefix + "_node.txt";
pubfilename = outprefix + "_pub.txt";
ancfilename = outprefix + "_anc.txt";

print "=======================================================================";
print "\t\tCAFE Report File Analysis"
print "\t\t" + cafecore.getDateTime();
print "---------";
print "Parsing format information...\n";

infile = open(infilename, "r");
inlines_main = infile.readlines();
infile.close();
# Reads the input report file.

if inlines_main[2].find("Lambda tree:") != -1:
	treeline = inlines_main[3];
	formatline = inlines_main[4];
	avgline = inlines_main[6];
	linestart_main = 11;
else:
	treeline = inlines_main[2]
	formatline = inlines_main[3];
	avgline = inlines_main[5];
	linestart_main = 10;
# If CAFE was run with a lambda tree structure as input, the report file places that on the third
# line. This shifts all the other relevant lines down by 1. This if/else accounts for that.

labeled_tree = treeline[treeline.index(":")+1:].strip();
tinfo, newtree = cafecore.treeParseNew(labeled_tree,2);
# This reads the CAFE tree with its node labels.

formatline = formatLineParse(formatline);
# formatline is CAFE's line with its paired node format with node ids. The formatLineParse function
# reads that format and returns it as a list of lists.

avgline = avgline.split(":\t")[1].strip().replace("\t", " ");
avgline = formatLineParse(avgline);
# The line of average expansions for each node, in the paired node format. Again passed to formatLineParse
# to make it interpretable.

print "---------";
print "Initializing output structures...\n";

node_fams_main = {"total" : []};
results_main = {};

sorted_nodes = [];
ancfile = open(ancfilename, "w");
header = "Family ID\t";

for node in tinfo:
	if tinfo[node][3] == 'root':
		continue;
	node_fams_main[node] = [[],[]];
	results_main[node] = [0,0,0,0,0,0,0,0,0,0];

	sorted_nodes.append(node);
	header += node + "\t";

ancfile.write(header[:-1] + "\n");
ancfile.close();
# [expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigcontract,total sig changes]
# node_fams and results are the two main dictionaries to store CAFE's results.
# node_fams {key:value} format: {node:list of two lists containing family ids for (rapid) expansions and (rapid) contractions, respectively}
# results {key:value} format: {node:[expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigcontract,total sig changes]}
# This loop also does the prepping of the header and sorted nodes for the anc count file

for j in range(len(formatline)):
	for k in range(len(formatline[j])):
		n = "<" + formatline[j][k] + ">";
		for r in results_main:
			if n in r:
				results_main[r][6] = avgline[j][k];
# Setting average expansion in results as read from avgline

print "---------";
print "Counting changes per branch...\n";

results_main, node_fams_main = cra(inlines_main, results_main, node_fams_main, linestart_main, ancfilename, sorted_nodes, 1);

if largefilename != None:
	print "\n\n---------";
	print "Parsing large families...\n";

	lfile = open(largefilename, "r");
	llines_main = lfile.readlines();
	lfile.close();
	# Reads the input report file.

	if llines_main[2].find("Lambda tree:") != -1:
		linestart_main = 11;
	else:
		linestart_main = 10;

	results_main, node_fams_main = cra(llines_main, results_main, node_fams_main, linestart_main, ancfilename, sorted_nodes, 0);

print "\nDone!";
print "=======================================================================";

print "Writing output files...";

## Begin fam output block.
outfile = open(famfilename, "w");
outfile.write("");
# Initialize the output file.
outfile.write("# The labeled CAFE tree:\t" + labeled_tree + "\n");

if r_opt == 0:
	desc_str = " ";
elif r_opt == 1:
	desc_str = " rapid ";

outline = "Overall" + desc_str + ":\t"
for f in node_fams_main['total']:
	outline = outline + f + ",";
outline = outline[:-1] + "\n";
outfile.write(outline);

for spec in node_fams_main:
	if spec == 'total':
		continue;

	# for f in range(len(node_fams_main[spec])):
	# 	if f == 0:
	# 		outline = spec + desc_str + "expansions:\t";
	# 	elif f == 1:
	# 		outline = spec + desc_str + "contractions:\t";

	# 	for rapid_f in node_fams_main[spec][f]:
	# 		outline = outline + rapid_f + ",";
	# 	outline = outline[:-1] + "\n";
	# 	outfile.write(outline);
	# For output on separate lines for expansions and contractions per species.

	outline = spec + ":\t";
	outline += ",".join(node_fams_main[spec][0] + node_fams_main[spec][1]) + "\n";
	outfile.write(outline);
	# For output on a single line per species.

outfile.close();
## End fam output block


## Begin node and pub output block
nodefile = open(nodefilename, "w");
nodefile.write("Node\tExpansions\tContractions\tRapidly evolving families\n");

pubfile = open(pubfilename, "w");
pubfile.write("Species\tExpanded fams\tGenes gained\tgenes/expansion\tContracted fams\tGenes lost\tgenes/contraction\tNo change\tAvg. Expansion\n");

for node in results_main:
	outline = node + "\t" + str(results_main[node][0]) + "\t" + str(results_main[node][3]) + "\t" + str(results_main[node][9]) + "\n";
	nodefile.write(outline);

	if node.replace("<","").replace(">","").isdigit():
		continue;

	spec = node.title()[:node.index("<")];
	exp = results_main[node][0];
	con = results_main[node][3];
	outline = spec + "\t" + str(results_main[node][0]) + " (" + str(results_main[node][7]) + ")\t" + str(results_main[node][1]) + "\t";
	if exp != 0:
		outline += str(round(float(results_main[node][1])/float(exp),2)) + "\t";
	else:
		outline += '0' + "\t";
	outline += str(results_main[node][3]) + " (" + str(results_main[node][8]) + ")\t" + str(results_main[node][4]) + "\t";
	if con != 0:
		outline += str(round(float(results_main[node][4])/float(con),2)) + "\t";
	else:
		outline += '0' + "\t";
	outline += str(results_main[node][2]) + "\t" + str(results_main[node][6]) + "\n";
	pubfile.write(outline);

nodefile.close();
pubfile.close();
## End node and oub output block

## Begin plot block
# print "Generating plots...";

# x_nodes = [];
# y_rapids = [];
# y_changes = [];
# y_exp = [];
# y_rexp = [];
# y_con = [];
# y_rcon = [];
# for node in results_main:
# 	y_rapids.append(results_main[node][9]);
# 	y_changes.append((results_main[node][0]+results_main[node][3])-results_main[node][9]);

# 	y_rexp.append(results_main[node][7]);
# 	y_exp.append(results_main[node][0]-results_main[node][7]);

# 	y_rcon.append(results_main[node][8]);
# 	y_con.append(results_main[node][3]-results_main[node][8]);

# 	if node[0] != "<":
# 		node = node.title()[:node.index("<")];

# 	x_nodes.append(node);

# #barcols = ['#ffef52','#5b5bd7'];
# barcols = ['#e5653c', '#2aa064'];

# y_data = [y_rapids, y_changes];
# y_names = ['total rapids', 'changes'];
# crplot.barPlotStack(x_nodes,y_data,y_names,"","# changing families","# of changing families",outprefix+"_change.html",barcols,w=1200);
# # Total plot

# y_data = [y_rexp, y_exp];
# y_names = ['rapid expansions', 'expansions'];
# crplot.barPlotStack(x_nodes,y_data,y_names,"","# expanding families","# of expanding families",outprefix+"_expand.html",barcols,w=1200);
# # Expansion plot

# y_data = [y_rcon, y_con];
# y_names = ['rapid contractions', 'contractions'];
# crplot.barPlotStack(x_nodes,y_data,y_names,"","# contracting families","# of contracting families",outprefix+"_contract.html",barcols,w=1200);
# # Contraction plot
## End plot block

# node_fams {key:value} format: {node:list of two lists containing family ids for (rapid) expansions and (rapid) contractions, respectively}
# results {key:value} format: {node:[expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigcontract,total sig changes]}
print "RESULTS TABLE -- tab delimted for easy copy/pasting into your favorite spreadsheet program"
print "\tExpansions\tGenes Gained\tEqual\tContractions\tGenes Lost\tFamilies Lost\tAverage Expansion\tSig Expansions\tSig Contractions\tTotal Sig Changes";
for species in results_main:
	outline = species + "\t";
	for col in results_main[species]:
		outline = outline + str(col) + "\t";
	print outline;

print
print "CAFE labeled tree:\t" + labeled_tree;
# This block simply prints the information stored in results to the screen.
print "=======================================================================";
```

