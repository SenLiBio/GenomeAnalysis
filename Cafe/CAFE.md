# CAFE - genefamily

# 1. 软件和自带脚本包的下载

cafe的[github主页](https://github.com/hahnlab/CAFE)、[最新版下载页面](https://github.com/hahnlab/CAFE/releases/tag/v4.2.1)、[官网](https://hahnlab.github.io/CAFE/download.html)、[自带的Python包](https://iu.app.box.com/v/cafetutorial-files)。


💡 有时候需要翻墙，然后Python包暂时看起来不可用，直接从这个页面的附件下载比较好。

---

# 2. 输入文件和conf文件准备

## 2.1 输入文件准备


💡 `Orthogroups.GeneCount.tsv`在*Orthofinder输出中的`Orthogroups`文件夹。

```bash
perl GeneCount2CafeInput.pl ../Orthogroups/Orthogroups.GeneCount.tsv >cafeinput.raw
python clade_and_size_filter.py -s -i cafeinput.raw -o cafeinput.filtered
```

- GeneCount2CafeInput.pl
  
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
    
- clade_and_size_filter.py
  
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

## 2.2 运行cafe

```bash
cafe cafe_script.sh
```

- cafe_script.sh
  
    ```bash
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

## 2.3 结果统计

```bash
python /path2/python_scripts/cafetutorial_report_analysis.py -i reportout1.cafe -r 1 -o rapid
```

---

## 2.4 用拟南芥的基因做基因家族功能注释

```perl
perl Rapid_GF_Anno.pl rapid_fams.txt Orthogroups.txt ath.gene.discription
```

-   Rapid_GF_Anno.pl

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

    

