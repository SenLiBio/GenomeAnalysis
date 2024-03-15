# psmc

```markdown
samtools mpileup -C50 -uf ref.fa aln.bam | bcftools view -c - \   | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > diploid.fq.gz
utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa
utils/psmc_plot.pl -R diploid diploid1.psmc diploid2.psmc  
<!--（其中diploid是接下来生成文件的前缀，diploid.psmc可以将所有生成的psmc文件都列上，结果会生成相应的txt、eps、gp文件，-g参数更改多少年一代，-u参数更改突变率，-w设置线的粗细）-->
vim diploid.gp<!--（根据说明选择相应的线型和颜色）-->           
gnuplot diploid.gp<!--（将生成相应的eps文件）-->   
epstopdf diploid.eps<!--（生成pdf）-->
```

