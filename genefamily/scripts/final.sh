# 1. 脚本放在**`./scripts`**目录下
# 2. 记得更改`1.extract_ortholog.pep.pl`和`3.extract_ortholog.cds.pl`中的物种数目（分别在第47行和第52行）
# 3. 记得检查`5.create_pal2aln.pl`中pal2nal的路径（第11行）
# 4. 记得检查`6.create.gblocks_sh.pl`中gblocks的路径（第9行）
# 5. 记得更改最后的jModelTest.jar的路径
# 6. 由于用Windows或者某些软件打开后，换行符那里会添加一个`^M`,linux识别的时候可能会出现问题，建议先用`cat -A [final.sh](http://final.sh)` 检查换行符，然后用`sed -e 's/^M//g' final.sh >1.bak | mv 1.bak final.sh` 来去除^M，注意打^M的方式是`Ctrl+M`而不是分别打出^和M。



# perl scripts/1.extract_ortholog.pep.pl ./Orthogroups/Orthogroups.txt allpepfiles Single_Copy_Orthologue_Sequences &&
perl scripts/2.create_sh.mafft.pl Single_Copy_Orthologue_Sequences mafft.sh mafft_output_dir &&
bash mafft.sh &&
perl scripts/3.extract_ortholog.cds.pl ./Orthogroups/Orthogroups.txt all_cds_dir Single_Copy_cds_dir length_check.txt &&
perl scripts/4.sort_cds_mafft.pl mafft_output_dir Single_Copy_cds_dir Single_Copy_cds_sorted_dir &&
perl scripts/5.create_pal2aln.pl mafft_output_dir Single_Copy_cds_sorted_dir pal2nal_out_dir pal2nal.sh &&
bash pal2nal.sh &&
perl scripts/6.create.gblocks_sh.pl pal2nal_out_dir gblocks_log_dir gblocks.sh &&
sh gblocks.sh && 
mkdir gblocks_out_dir &&
mv pal2nal_out_dir/*htm gblocks_log_dir &&
mv pal2nal_out_dir/*gb gblocks_out_dir &&
perl scripts/7.cat_cds.pl gblocks_out_dir merged_aligned_cds_dir &&
cat merged_aligned_cds_dir/*fa >final_mafft_gblocks_aligned.fa &&
sh scripts/convertFasta2Phylip.sh final_mafft_gblocks_aligned.fa >final_mafft_gblocks_aligned.phy &&
java -jar -XX:ParallelGCThreads=4 -Xmx4g ~/software/jmodeltest-2.1.10/jModelTest.jar -tr 10 -d final_mafft_gblocks_aligned.phy -s 11 -f -i -g 8 -AIC -BIC -AICc -DT -p -a -w -o bestmodel