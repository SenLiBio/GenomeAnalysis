# RepeatModeler+RepeatMasker



# 1. 安装与配置

这个软件的安装还是比较麻烦的，建议参考[教程1](https://www.jianshu.com/p/8c20f7922f90)和[教程2](https://www.jianshu.com/p/ffdbedae80fa)，当然后续的使用也可以参考这个文件。需要注意的一点是，在选择repeatmasker和repeatmodeler调用的其他软件路径的时候，最好是选择该软件自己的安装路径的bin文件夹下的运行文件，因为这个软件调用的某些模块可能是从自己的安装路径找，如果简单的选`/pathto/miniconda2/bin`可能会出现运行问题。

---

# 2. RepeatModeler构建重复序列的library

具体命令可以参考[官网](http://www.repeatmasker.org/RepeatModeler/)。

## 2.1 构建数据库

```bash
BuildDatabase -name 001brgy 001brgy.fa 2>>bgy.repeatmodeler.err
# 可以用“-engine”参数指定blast软件，结果相差不多，具体可以参考官网说明，一般可以用“-engine ncbi”
```

## 2.2 构建library

```bash
RepeatModeler -pa 10 -database 001brgy -LTRStruct >& run.out
# -pa,  用到的线程数
# -LTRStruck, 同时注释ltr，这一步要求安装LTRHarvest等软件
# -engine, blast用的软件
```

---

# 3. TEclass鉴定unknown的重复元件

## 3.1 安装和配置

具体可以参考[网上教程](https://www.jianshu.com/p/8d20a35330db)，这个软件真的很老旧了。。遇到过跟最新的repbase不兼容的问题，跟作者交流后也没有好的解决方法，不过对于自己构建的library还是可以用的。

## 3.2 提取Unknown序列进行分类

记得现在目录下创建相应的perl脚本。

- teclass.sh
  
    ```bash
    # teclass.sh
    #USAGE: sh this_script custom
    #custom must be 123wxyz format (3num4letters)
    mkdir 0out_rmd 1known_rmd 2unknown_rmd 3out_tecls 4transOut_tecls 5ultRslt
    ln -s ../$1-families.fa 0out_rmd/
    for i in ./0out_rmd/*.fa;do echo $i;perl seperate_denoPrdcTEconsn.pl $i;done 
    
    cd 2unknown_rmd/
    for i in *.fa; do numbAbbr=`echo $i|sed 's/.*\([0-9]\{3\}....\).*\.fa/\1/'`; TEclassTest -o ../3out_tecls/$numbAbbr $i>outstd.teclass_$numbAbbr 2>errstd.teclass_$numbAbbr; done
    
    cd ../
    
    for i in 3out_tecls/*/*-families.fa.lib;do echo $i;perl teclass_title_transition.pl $i;done
    
    for i in ./1known_rmd/*.fa; do file=`echo $i|sed 's/.*\([0-9]\{3\}.*\.fa\)/\1/'`; cat ./1known_rmd/$file ./4transOut_tecls/$file >./5ultRslt/$file; done
    
    #for fm in 5ultRslt/*-families.fa; do nw=`echo $fm|sed 's/.*[0-9]\{3\}\(....\).*ies\.fa/\1/'`;echo $nw; sed -i "s/^>rnd-/>$nw/" $fm;sed -i 's/family-/fmly/' $fm;done
    
    ```
    
- seperate_denoPrdcTEconsn.pl
  
    ```perl
    # seperate_denoPrdcTEconsn.pl
    #!/usr/bin/perl -w
    
    #used to seperate known and unknown TE consensus sequence 
    #output of RepeatModeler.
    #two directoried must exist when run this script!!!!!!!!!
    #they are 1known_rmd & 2unknown_rmd 
    
    #USAGE: perl this_script file.fasta
    
    use strict;
    my %id2seq=&hash_fasta($ARGV[0]);
    $ARGV[0]=~m/\d{3}\w{4}.*-families\.fa/;
    open KNN, ">1known_rmd/$&" or die "$!";
    open UNKNN, ">2unknown_rmd/$&" or die "$!";
    foreach(sort{$a cmp$b}keys%id2seq){
            if(/#Unknown/){
                    print UNKNN ">$_"."$id2seq{$_}";
            }else{
                    print KNN ">$_"."$id2seq{$_}";
            }
    }
    
    sub hash_fasta{#scaffold id to sequence hash generator
            open FA, "$_[0]" or die "$!";
            my %fa_hash=(); 
            my $id=""; 
            my @ids=();     
            while (<FA>){
                    if (/^>/){
                            s/^>//;
                            my@_1=split "\\s+", $_;
                            #$id=$_1[0];
                            $id=$_;
                            push @ids, $_1[0];
                    }else{
                            #s/\s+$//g;
                            $fa_hash{$id}.=$_;
                    }
            }close FA;
            my @dup=grep{$_{$_}++}@ids; %_=();
            die "there are duplicated sequence ID(@dup) in $_[0] file" if @dup>0;   
            %fa_hash;
    }
    ```
    
- teclass_title_transition.pl
  
    ```perl
    # teclass_title_transition.pl
    #!/usr/bin/perl -w
    #used to trans id line of teclass out file.
    #Unknown will be substibuted to proper 
    #TE type(one of 'DNA LINE LTR or SINE') if classification done
    #make sure directory 4transOut_tecls exist, see line 14 which is
    #'open OUT, ">./4transOut_tecls/$nfn";'
    
    #USAGE: perl this_script TEclass_out.lib 
    use strict;
    my@te=qw/DNA LINE LTR SINE/;
    open TECLSLIB, "<$ARGV[0]"or die "$!";
    $ARGV[0]=~m/.+\/(\d{3}\w{4}.*?\.fa)\.lib/;
    my$nfn=$1;#new file name
    open OUT, ">./4transOut_tecls/$nfn";
    my$sum=0;
    my%te_counter=();
    while(<TECLSLIB>){
            my $new_idline='';
            my $tecls='';
            if(m/(.+?)\|TEclass result: (.+?)\|/){
                    $new_idline=$1;
                    $tecls=$2;
                    if($tecls~~@te){
                            $new_idline=~s/#Unknown/#$tecls/;
                            print OUT "$new_idline\n";
                    }else{print OUT "$new_idline\n";
                    }
            }else{print OUT $_;
            }
    }
    close OUT;
    close TECLSLIB;
    ```
    

```bash
mkdir teclass
cd teclass
bash teclass.sh 001brgy
cd ..
```

---

# 4. RepeatMasker做鉴定和汇总

对结果的解读可以参考[官网](http://repeatmasker.org/webrepeatmaskerhelp.html)

```bash
ln -s ./teclass/5ultRslt/001brgy-families.fa 001brgy-families_tecls.fa
mkdir ./repeatmasker
RepeatMasker -pa 4 -html -gff -x -poly -dir ./repeatmasker -lib ./"001brgy-families_tecls.fa" 001brgy.fa >./outstd_RMs_001brgy 2>./errstd_RMs_001brgy

# -pa, 线程
# -html, Creates an additional output file in xhtml format
# -gff, Creates an additional Gene Feature Finding format output
# -x, Returns repetitive regions masked with Xs rather than Ns
# -poly, Reports simple repeats that may be polymorphic (in file.poly)
```

---

# 5. shell scripts

记得在目录下创建my_scripts文件夹，并把`teclass.sh`和对应的两个perl脚本放进去。后续可以通过`customLib_defSpd_$1/001brcy.fa.masked`文件里面的X的数目来统计不同step window中的重复序列比例来画circos图。

```bash
#USAGGE: sh this_script custom genome.fasta
#the parallel version of repeatm.sh 
#custom must be 123wxyz format (3num4letters) 
BuildDatabase -name $1 -engine ncbi $2 >./outstd_BD_$1 2>./errstd_BD_$1
RepeatModeler engine ncbi -pa 4 -database $1 >./outstd_RMd_$1 2>./errstd_RMd_$1
mkdir ./teclass
cd ./teclass
sh ./my_script/teclass.sh  $1
cd ..
ln -s ./teclass/5ultRslt/$1-families.fa $1-families_tecls.fa
mkdir ./customLib_defSpd_$1
RepeatMasker -html -gff -x -poly -pa 4 -dir ./customLib_defSpd_$1 -lib ./"$1-families_tecls.fa" $2 >./outstd_RMs_$1 2>./errstd_RMs_$1
```