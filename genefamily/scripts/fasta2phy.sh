#!/usr/bin/perl
use strict;
use warnings;

# 检查输入参数
if (@ARGV != 1) {
    die "USAGE: ./script.pl <fasta-file>\n";
}

# 打开输入文件
my $file = $ARGV[0];
open(my $fh, '<', $file) or die "Could not open file '$file': $!";

my @lines = <$fh>;
close($fh);

# 初始化变量
my $num_species = 0;
my $all_sequences = "";
my @output_lines;

foreach my $line (@lines) {
    chomp $line;

    if ($line =~ /^>(\S+)/) { # 提取种名
        my $species_name = $1;
        $num_species++;
        push @output_lines, "$species_name  "; # 添加种名，后面两个空格
    } else { # 处理序列
        $output_lines[-1] .= $line; # 将序列附加到当前种名后
        $all_sequences .= $line;   # 统计所有序列
    }
}

# 计算序列总长度
my $sequence_length = length($all_sequences);

# 输出结果
print "$num_species $sequence_length\n";
print join("\n", @output_lines), "\n";
