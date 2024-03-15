# Documentary

# 1. 一些circos的基本知识：

1. circos读入的主体是circos.conf文件，这个文件里面指定了一些基础的画图信息，以及从哪些文件读取进一步的画图信息。建议把所有用到的文件放在“./etc”路径下。
2. 要搞清楚ideogram和chromosome的区别，ideogram说白了就是画出来的chromosome的部分，基础用法当然是一个ideogram对应一个chromosome，但是也可以通过命令只画部分chromosome，或者将一条chromosome拆分成几个ideogram
3. 有些参数后面会带有一个“*”，这个符号表示这里设置的参数值会覆盖掉任何之前设置的同样的参数。同样的，“**”会覆盖掉“*”
4. 如果想画某些类型的图，可以去看下官网的[documentary](http://www.circos.ca/documentation/tutorials/quick_start/)

---

# 2. 安装和使用

## 2.1 Installation and run

```bash
# 最方便的方法就是通过conda来安装
conda create -n circos
conda install -c bioconda circos
circos
```

---

## 2.2 circos.conf

## 2.2.1 基本设置和karyotype.txt

这里的是circos.conf文件里面的基础设置，并且指定了其他conf文件的位置。粉色底的部分放在文件开头，黄色底的部分放在文件结尾。karyotype.txt文件的例子可以参考[karyotype.txt](../circos%20be3c7.md)

```
# the karyotype parameter specifies the file which defines the
# size and name of each chromosome for the figure
karyotype = /path_to_the_file/karyotype.txt
# 同时画两个基因组的染色体
# karyotype = /path_to_the_file/karyotype1.txt,/path_to_the_file/karyotype2.txt

***# Always use “;” as the delimiter between records in the chromosomes parameter***
# unit of length for the chromosomes – this is used
# in other parts of the file where a position is referenced
chromosomes_units = 1000000

# toggle to display all of the chromosomes in the
# karyotype file, in the order of appearance
chromosomes_display_default = yes

# spacing, position, thickness, color/outline format of the ideograms. see [2.3.1 ideogram.conf](../circos%20be3c7.md)
# from same directory as circos.conf
<<include ideogram.conf>>

# position, labels, spacing and grids for tick marks
# from shared session directory. see [2.3.3 ticks.conf](../circos%20be3c7.md)
<<include ../etc/ticks.conf>>

# size, format and location of the output image, shared by
# all the sessions in this course
<<include ../../../etc/image.conf>>

# the files below are included from the Circos distribution
# DO NOT REMOVE THESE
# defines colors, fonts and fill patterns
<<include etc/colors_fonts_patterns.conf>>

# the housekeeping file contains system parameters that
# define file delimiters, debugging options and other global settings
<<include etc/housekeeping.conf>>

```

```bash
***# karyotype.txt***
# A simple karyotype with 5 chromosomes:
#
# chr1 5Mb
# chr2 10Mb
# chr3 20Mb
# chr4 50Mb
# chr5 100Mb
#
# The format of this file is
#
# chr - CHRNAME CHRLABEL START END COLOR
#
# In data files, chromosomes are referred to by CHRNAME.
# On the image, they are labeled by CHRLABEL
#
# Colors are taken from the spectral Brewer palette.
# To learn about Brewer palettes, see
#
# www.colorbrewer2.org
# mkweb.bcgsc.ca/brewer
# 色板里面大多数最多只有12种颜色可选，也可以用UCSC人类基因组的配色，颜色直接设置为chr1，chr2...chr22，chrY，chrX，chrM，chrNA，chrUn。
chr - chr1 1 0 5000000 spectral-5-div-1
chr - chr2 2 0 10000000 spectral-5-div-2
chr - chr3 3 0 20000000 spectral-5-div-3
chr - chr4 4 0 50000000 spectral-5-div-4
chr - chr5 5 0 100000000 spectral-5-div-5

***# Add cytogenetic bands to the ideogram***
# the first field must be “band”
# the second field is the chromosome with which the band is associated
# next two fields are the name & label of the band (not used, but must be present)
# the last three fields are the start, end and color of the band. You can use any color for the bands, but conventionally colors gposNNN (e.g., gpos25, gpos50, gpos100),
gneg, gvar, acen and stalk are used.
band chr1 band1 band1 0 2500000 gneg
band chr1 band2 band2 2500000 5000000 gpos25
band chr2 band1 band1 0 2500000 gneg
band chr2 band2 band2 2500000 5000000 gpos25
band chr2 band3 band3 5000000 7500000 gpos100
band chr2 band4 band4 7500000 10000000 gvar
band chr3 band1 band1 0 1000000 stalk
band chr3 band2 band2 1000000 5000000 gpos50
band chr3 band3 band3 5000000 7500000 gpos100
band chr3 band4 band4 7500000 10000000 gvar
band chr3 band5 band5 10000000 15000000 acen
band chr3 band7 band7 15000000 19000000 gneg

***TRANSPARENT COLORS***
# By changing band_transparency in the [<ideogram> block](../circos%20be3c7.md), the band pattern can be made semitransparent to allow the color of the ideogram to show through. Using this parameter you can combine the ideogram color and the band pattern
```

---

## 2.2.2 chromosome和ideogram设置

### 2.2.2.1 选择要展示哪些chromosome的位置

```
***## 1. filter chromosomes***
# show only chromosomes 1, 2 and 3
#chromosomes_display_default = no
#chromosomes = chr1;chr2;chr3

# remove some chromosomes, you can use regular expression instead of ideogram name
#chromosomes_display_default = yes
#chromosomes = -hs1;-hs2;-mm1;-mm2
#chromosomes = -/[XY]/
#chromosomes = -hs1;-mm1;-/Y/;-/\d\d/
```

### 2.2.2.2 染色体颜色设置

```
***## 2. color of chromosomes*
# Both “=” and “:” can be used as the assignment operator (e.g. /hs/:reds-5-seq-5 and /hs/=reds-5-seq-5).
****chromosomes_display_default = yes
chromosomes_color = /hs/=reds-5-seq-5,/mm/=blues-5-seq-5
# chromosomes_color = hs1=rdylbu-11-div-2,hs2=rdylbu-11-div-3,mm1=rdylbu-11-div-10,mm2=rdylbu-11-div-9
```

```
***## 3.1 Order of chromosomes***
# explicitly define order
#chromosomes_order = chr1,chr2,chr5,chr4,chr3
# relative order. Only the order of chromosomes 3 and 5 is affected. In this case, when the value of the parameter lists a subset of chromosomes, the first chromosome in the list acts as an anchor (its order is not affected – it appears in the same position as it would by default) and any subsequent chromosomes listed in the value appear next.
#chromosomes_order = chr3,chr5
# relative order，In this case there are two subsets that are reordered: chromosomes 1,4 and chromosomes 2. In the first subset chromosome 1 acts as the anchor (it is placed first, as it would by default) and followed by chromosome 4. The last chromosome is chromosome 2 and the two “-“ correspond to the remaining chromosomes whose order has not been explicitly stated (chromosomes 3 and 5).
#chromosomes_order = chr1,chr4,-,-,chr2

***## 3.2 Ordering ideogram regions***
# ***Tags*** are delimited by ***[]*** and are defined in the chromosomes parameter. In the example below, chromosome 4 is associated with two ideograms. The first ideogram, named ‘a’, represents the region 0-11 Mb. The second ideogram, named ‘b’, represents 19-31 Mb.
chromosomes = chr3:0-11;
							chr4[4a]:0-11;
							chr4[4b]:19-31;
							chr5[5a]:9-21;
							chr5[5b]:29-41;
# Once ideograms of a chromosome are assigned tags, the tags can be used to reorder the ideograms, ***^*** character represents***the start of the ideogram circle***.
chromosomes_order = ^,chr1,4a,chr3,4b,-5b,5a
```

### 2.2.2.4 染色体的尺度调整

```
***## 4. scale of chromosomes***
# Ideograms formed by chromosomes that have been broken into pieces can be individually scaled by using the tag name of the ideogram.
chromosomes_scale = 5a=2,5b=2

chromosomes_scale = chr1=0.2,chr2=2,chr3=10
# which would display chromosome 1 at 0.2× magnification (5× reduction), chromosome 2 at 2× magnification and chromosome 3 at 10× magnification, respectively.

#chromosomes_scale = chr1=0.5r
# chr1 occupies 50% of figure

#chromosomes_scale = /chr[123]/=0.5rn
# This will apply scale to chr1, chr2 and chr3, which are the chromosome names that match the regular expression. The special rn suffix indicates normalized relative scale. Half of the image (0.5r) will be assigned to the m chromosomes that match the regular expression, and each chromosome will occupy 0.5/m of the image

#chromosomes_scale = /./=1rn
# all the chromosomes to appear to be the same size

# use regular expression to set the scale of chromosomes
#chromosomes_scale = /hs/=1rn
#chromosomes_scale = /hs/:0.5rn,/mm/:0.5rn

# reverse the scale of mm1 and mm2 chromosomes
#chromosomes_reverse = mm1,mm2
```

### 2.2.2.5 将染色体分为不同的ideogram

```
***## 5.Drawyhning ideogram regions***
# Notice the region of chromosome 5 is expressed in units of
chromosomes_units, to avoid having many non-significant zeroes in the configuration file. Since the unit is 1 Mb, the range is 9-21 Mb
chromosomes = chr3=8-12;chr4=4-11;chr5=9-21

***## 6.Chromosome breaks** 
# The format of the axis break is controlled in the <ideogram><spacing> block. Add one line "**<<include break.conf>>**" in the **ideogram.conf** file and see session **[break.conf](../circos%20be3c7.md)** for detail.
*****# to remove the 11-19 Mb region of chromosome 5 from the display, set
chromosomes_breaks = -chr5:11-19
# move multiple breaks
# chromosomes_breaks = -chr3:13-17;-chr4:(-9;-chr4:41-);-chr5:11-19
# you can use *“(“ and “)”* to represent the ***start (or end)*** of the chromosome. Thus 0-9 is the same as (-9. This notation is helpful because you do not need to remember the start (or end) value of the chromosome to define a break at its edge.
# multiple breaks can be defined for a given chromosome. In this example, chromosome 4 has a break at its start (-9, which removes the first 9 Mb and at its end 41-), which removes
the last 9 Mb.
```

### 2.2.2.6 设置ideogram的高亮区域

```
***## 7.highlights***
# The data file defines genome regions (which may overlap) for which a colored segment is drawn in the figure. The highlight segments are drawn underneath any ticks and grids and are useful for annotating parts of the figure with colors (for lookup or focus).
# see [highlight.txt](../circos%20be3c7.md)
<highlights>
<highlight>
file = ../data/highlight.txt
r0 = 1r+40p
r1 = 1r+45p
</highlight>
</highlights>
# All regions found in the file 2/data/highlight.txt will be drawn as segments between radial positions r0 = 1r+40p and r1 = 1r+45p. This radial position format expands on expressions you’ve already seen (e.g. 1.125r) by adding an absolute value offset. The expression 1r+40p means a radial position obtained by adding 40 pixels to the outer rim of the ideogram circle. Similarly, 1r+45p means a position obtained by adding 45 pixels. 
# You can have any number of <highlight> blocks. The blocks may refer to highlight regions that overlap both radially and with respect to their genomic coordinates.
```

```
***# highlight.txt***
# chr start end options
hs1 0 249250621 fill_color=rdylbu-11-div-2
hs2 0 243199373 fill_color=rdylbu-11-div-3
mm1 0 197195432 fill_color=rdylbu-11-div-10
mm2 0 181748087 fill_color=rdylbu-11-div-9
```

---

## 2.2.3 data tracks

In this session we’ll learn how to add tracks in the image. 

**注意：这里的每一种track的添加方法都是在circos.conf文件中添加相应的block**

### 2.2.3.1 histograms

```

***## 8. histograms*
# note: 每一个<plots> block对应的是图里面的一个圈，之后的每一个<plot>对应的是这个圈里面的一条线，通过设置show = no可以选择不展示这条线，然后根据z的值可以控制图层的优先级，这样就可以展示出堆叠的效果,例如[3.1 histogram](../circos%20be3c7.md)
# 数据格式参考[2.4.3 histogram.txt](../circos%20be3c7.md)**
<plots>
# global settings for each <plot> block
type = histogram
thickness = 1p
color = black
#color = white
# min and max range of the data
min = 0
max = 1
# The histograms are all placed between inner radius r0 = 0.85r and outer radius r1 = 0.975r (relative to the ideogram
circle).
r0 = 0.85r
r1 = 0.975r

# All track types have sufficient parameters set to default values to be drawn by only specifying the file name and track type. This gets you started quickly. whose remaining necessary parameters are set by the /path_to_the_circos_software/etc/track/histogram.conf file in the Circos distribution directory.
<plot>
type = histogram
file = ../data/both.cons.2e6.avg.txt
</plot>

# plot track definition and overriding settings, the example of the data file is shown in [2.4.3 histogram.txt](../circos%20be3c7.md)
<plot>
file = ../data/both.cons.2e6.max.txt
fill_color = spectral-5-div-3 # yellow
</plot>

# another plot track and overriding settings
<plot>
# disable the track
show = no
file = ../data/both.cons.2e6.avg.txt
fill_color = spectral-5-div-4 # green
#thickness = 2p
# The order in which tracks are drawn is controlled by the z (think of it as depth) parameter. By default, each track has z = 0. Tracks are drawn in order of their z value; thus, tracks with higher z are drawn on top of tracks with a lower z. Individual data points can have their own z value as well.
z = 5
</plot>

<plot>
show = no
file = ../data/both.cons.2e6.min.txt
fill_color = spectral-5-div-5 # blue
#fill_color = white
z = 10
</plot>
</plots>

```

```
# chr start end value
mm1 10000000 11999999 0.439079022807017 color=blue
mm1 12000000 13999999 0.44630188700565 color=blue,thickness=3p
mm1 14000000 15999999 0.452856295597484 color=green,id=abc
```

![Documentar%2086d1d/Untitled.png](Documentar%2086d1d/Untitled.png)

![Documentar%2086d1d/Untitled%201.png](Documentar%2086d1d/Untitled%201.png)

![Documentar%2086d1d/Untitled%202.png](Documentar%2086d1d/Untitled%202.png)

![Documentar%2086d1d/Untitled%203.png](Documentar%2086d1d/Untitled%203.png)

---

### 2.2.3.2 heat maps

```
# 由于热图并不像曲线一样可以叠加,所以不存在<plots>这种组合的block，每一个热图的添加都是由一个单独的<plot>来指定。
# see [2.4.4 heatmap.txt](../circos%20be3c7.md) for the data format.
<plot>
type = heatmap
file = ../data/both.cons.2e6.rhe.avg.txt

# min and max range of the data
min = 0.1
max = 0.9

# The histograms are all placed between inner radius r0 = 0.85r and outer radius r1 = 0.975r (relative to the ideogram
circle).
r0 = 0.73r
r1 = 0.75r

# The colors of the heatmap can be defined by using an explicit list of colors (tedious).
colors = spectral-11-div-1,spectral-11-div-2,spectral-11-div-3,spectral-11-div-4,spectral-11-div-5,spectral-11-div-6,spectral-11-div-7,spectral-11-div-8,spectral-11-div-9,spectral-11-div-10,spectral-11-div-11
# or by using a color list (convenient)
color = spectral-11-div
# Reversing the colors in a heatmap is very easy, if a list was used. Just append –rev to the color name.
color = spectral-11-div-rev
# You can combine color lists to create longer lists. In this case, an 18-color red-blue color palette was generated.
color = reds-9-seq-rev,blues-9-seq
#You can interchangeably mix colors and lists together. Below, the spectral color list will be flanked by black. Smallest and largest values will be black
color = black,spectral-11-div,black

# By default, color mapping is done linearly within the [min,max] range of the heat map. In other words, the range is divided by the number of colors, and each color is assigned the same interval. When your input values data are highly skewed, a linear mapping results in most of the values assigned to a small number of colors. Here, a logarithmic mapping is helpful, because it distributes values among colors more uniformly (Figure 15). To apply logarithmic mapping, use the scale_log_base parameter. For scale_log_base < 1 the color index grows quickly at the beginning (emphasizing differences in small values) and slowly at the end (attenuating differences in large values). The opposite happens for scale_log_base > 1.
scale_log_base = 0.5
</plot>
```

```
# chr start end value
mm1 10000000 11999999 0.439079022807017
mm1 12000000 13999999 0.44630188700565
mm1 14000000 15999999 0.452856295597484
```

![Documentar%2086d1d/Untitled%204.png](Documentar%2086d1d/Untitled%204.png)

![Documentar%2086d1d/Untitled%205.png](Documentar%2086d1d/Untitled%205.png)

![Documentar%2086d1d/Untitled%206.png](Documentar%2086d1d/Untitled%206.png)

![Documentar%2086d1d/Untitled%207.png](Documentar%2086d1d/Untitled%207.png)

---

### 2.2.3.3 links

btw，circos提供了一个比较方便的工具来生成区域内link终点的密度分布曲线。详情可以参考session2 lesson5，或者[在线的document](http://circos.ca/documentation/tutorials/utilities/density_tracks/)。

```
<links>
<link>
# see [2.4.5 link.txt](../circos%20be3c7.md) for the data format
file = ../data/links.txt
# bezier_radius控制的是曲线的弧度，简单的来说可以认为这个值越大，圆心附近的位置越空旷。具体设置参考这篇教程：[http://www.circos.ca/documentation/tutorials/links/geometry/images](http://www.circos.ca/documentation/tutorials/links/geometry/images)
bezier_radius = 0r

# radius参数控制的是每一根线的起始和终止位置，比如0.5r就是指画links的圆的半径是ideogram的圆的半径的一半
radius = 0.5r

thickness = 1p
color = black
# Recall that suffixing a color name with _a + digit creates a transparent color. The degree of transparency depends on auto_alpha_steps. Smaller values mean that the color is more opaque. Notice that red and red_a0 are both the same color—a fully opaque red. When set auto_alpha_steps = 5, red_a0~a5 represent 0% ~ 83% transparency.
#color = black_a5
</link>
</links>
```

```
# links data format
# chr1 start1 end1 chr2 start2 end2 value
hs1 120788758 120834144 mm1 63903740 63957877 id=abc,color=blue
hs2 190247704 190320005 mm1 124448312 124506291 id=def,color=red
```

![Documentar%2086d1d/Untitled%208.png](Documentar%2086d1d/Untitled%208.png)

![Documentar%2086d1d/Untitled%209.png](Documentar%2086d1d/Untitled%209.png)

---

### 2.2.3.4 tiles

- Tiles (or cover elements) are an important data type in genomics. They are useful to visualize data that is interpreted as sampling or coverage process (reads, clones, alignments, genes, etc).
- For tiles, elements are stacked in layers to avoid overlap. You do not have direct control in which layer within the track the element will appear. Several parameters control element placement, including element padding and margin and how elements that cannot fit within the track are to be handled.

```
<plot>
type = tile
file = ../data/tiles.txt
r0 = 1r+2p
r1 = 1r+40p
# The number of layers, Elements in the same layer cannot be closer than margin distance.
layers = 7
layers_overflow = hide
layers_overflow_color = red
# the margin distance
margin = 1u
# The radial size of each tile
thickness = 3
# the distance between
layers
padding = 2
# Tiles can stack inward (orientation=in) or outward (orientation=out).
orientation = out
color = black
stroke_thickness = 1
stroke_color = vdgrey
<backgrounds>
<background>
y0 = 0.75r
color = grey_a1
</background>
<background>
y0 = 0.25r
y1 = 0.75r
color = grey_a3
</background>
<background>
y1 = 0.25r
color = grey_a5
</background>
</backgrounds>
</plot>
```

```
# pile data format
# chr start end
mm1 10000000 11999999
mm1 12000000 13999999
mm1 14000000 15999999
```

![Documentar%2086d1d/Untitled%2010.png](Documentar%2086d1d/Untitled%2010.png)

![Documentar%2086d1d/Untitled%2011.png](Documentar%2086d1d/Untitled%2011.png)

---

### 2.2.3.5 highlights

与2.2.2.6里面的highlight设置不同，这里设置的是以track形式存在的高亮区域。

```
<plot>
type = highlight
file = ../data/highlight.max.top20.txt
r0 = 0.975r
r1 = 0.995r
# reference the inner and outer radius of the ideograms to draw the highlights on top of ideograms
#r0 = dims(ideogram,radius_inner)+5p
#r1 = dims(ideogram,radius_outer)-5p
# 填充的颜色，如果不设置填充的颜色，直接把highlight位置设置和其它track一样，那么可以做到直接在别的track上画高亮框
fill_color = spectral-5-div-3
# 边框的粗细
stroke_thickness = 1p
# 边框的颜色
stroke_color = red
z = 15
</plot>
<plot>
type = highlight
file = ../data/highlight.min.top20.txt
r0 = 0.835r
r1 = 0.855r
#r0 = dims(ideogram,radius_inner)+5p
#r1 = dims(ideogram,radius_outer)-5p
fill_color = spectral-5-div-4
stroke_thickness = 1p
stroke_color = red
z = 15
</plot>
```

```
# highlight data format
# chr start end
mm1 10000000 11999999
mm1 12000000 13999999
mm1 14000000 15999999
```

![Documentar%2086d1d/Untitled%2012.png](Documentar%2086d1d/Untitled%2012.png)

![Documentar%2086d1d/Untitled%2013.png](Documentar%2086d1d/Untitled%2013.png)

注意高亮度是在ideogram内测的trac附近的那些红框小块

---

### 2.2.3.6 scatter plots

Scatter plots use the same input format as histograms, but draw the data using glyphs. Glyph
size and color can be controlled, as can glyph shape. You can probably guess that we’ll be using
rules to adjust glyphs!

```
<plot>
type = scatter
file = ../data/scatter.cons.txt
min = 0.39
max = 0.55
r0 = 0.80r
r1 = 0.90r
glyph = square
# 形状，也可以是circle等
glyph_size = 3
#glyph_size = eval(remap_round(abs(var(value) - 0.455),0,0.1,5,20)
# 大小
fill_color = grey
```

![Documentar%2086d1d/Untitled%2014.png](Documentar%2086d1d/Untitled%2014.png)

## 2.2.4 Introduction to Rules

**注意：rules-conditions的组合相当于是circos的判断和转换语句。**

### 2.2.4.1 基础的rules

```
# 直接赋值
thickness = 3p
# perl code里面的eval，min这种参数
thickness = eval( min(5,int(var(size1)/10e3) ) )
```

### 2.2.4.2 单次判断

```
<links>
<link>
...
<rules>
<rule>
# each link is tested with a rule
# if the condition passes...
condition = on(hs1)
# the color will be changed
color = rdylbu-11-div-2_a3
# as well as the z value
z = 10
</rule>
</rules>
</link>
</links>
```

The rule shown above will change the color all links that start on human chromosome 1 (hs1).
The assigned color will be rdylbu-11-div-2, which is the color used for the chromosomes
color index.

![Documentar%2086d1d/Untitled%2015.png](Documentar%2086d1d/Untitled%2015.png)

---

### 2.2.4.3 rule chain

可以设置多个rule，这些rule会按照设置的先后顺序进行判断，一旦某个为真，就会自动短路，有点像for-elsif-else循环的感觉。也可以在判断成功后用flow=continue这个命令来继续进行接下来的判断。

```
<rule>
condition = on(hs1)
color = rdylbu-11-div-2_a3
z = 10
</rule>
# any link that passed the above rule is not tested further
<rule>
# this rule always matches, since its condition is always true
# any link that filed the previous rule is matched by this one
condition = 1
# the link will be hidden
show = no
</rule>
```

![Documentar%2086d1d/Untitled%2016.png](Documentar%2086d1d/Untitled%2016.png)

也可以在一个rule里面设置多个condition，来进行判断。

```
condition = between(hs1,mm1)
condition = var(start2) < 40mb || var(start2) > 160mb
# gb，mb，kb后缀分别代表十的九、六、三次方
```

![Documentar%2086d1d/Untitled%2017.png](Documentar%2086d1d/Untitled%2017.png)

---

### 2.2.4.4 动态取值

**注意：动态取值为设置color、thickness等参数提供了一个非常好的方法，这样就可以根据data来进行参数的调整，比如为起点在染色体的不同位置的线设置不同颜色。更多的关于rule的教学可以参考在线的[tutorials](http://www.circos.ca/documentation/tutorials/links/)。**

You can write rules that use the properties of the data point to set parameter values. These are dynamic rules—their behavior changes depending on the data.

Here we’ll write a rule that changes the thickness, color and order of a link. To use Perl code in
a rule, you must delimit the code using eval(). For example,

```
thickness = eval(2+2 . “p”)
# is the same as
thickness = 4p
```

Let’s make the thickness a function of the size of the start of the link. We’ll want to map the thickness to a range between 1 and 5. The following will do the job, and although it seems complicated at first I will walk you through the expression.

```
thickness = eval( sprintf(“%dp”,remap_round(var(size1),0,50000,1,5) )
```

The sprintf() function takes a format string and a list of values which are interpolated into the format string. Special tokens (e.g. %d) in the format string correspond to different ways in which numbers (and strings in general) can be formatted. For example %d means integer, created by truncating any fractional component, so that sprintf(“%d”,2.2) will return 2. Other useful tokens are %s (string) and %f (float).

In the above rule, the format string is %dp which produces an integer followed by a “p”, to express the value in units of pixels. Note that the “p” is not part of the special format token. The value to be formatted is provided by the remap_round function, available to you in all rules. This function takes a value (first argument), and linearly remaps it from (min,max) (second and third arguments) onto the range (remap_min,remap_max) (fourth and fifth arguments)

```
remap_round(value,min,max,remap_min,remap_max)
```

Thus, the size of the start of the link size1 will be remapped from 0-50kb to 1-5. Any values of
size1 > 50kb will be remapped to 5.

The rules for color and z work similarly. We will remap size onto z and position onto colors rdylbu-11-div-{1...11}_a3.

```
<rule>
condition = fromto(hs1,mm1)
condition = var(start2) < 40mb || var(start2) > 160mb
thickness = eval(sprintf("%dp",remap_round(var(size1),0,50000,1,5)))
z = eval(sprintf("%dp",remap_round(var(size1),0,50000,1,5)))
color = eval(sprintf("rdylbu-11-div-%d_a3",
remap_round(var(size1),0,250e6,1,11)))
</rule>
```

You can compare how hs1/mm1 and hs2/mm2
are related by applying the same rule to links between hs1 and mm2.

```
<rule>
# you can toggle whether a rule is considered using ‘use’
condition = fromto(hs1,mm1) || fromto(hs2,mm2)
condition = var(start2) < 40mb || var(start2) > 160mb
thickness = eval(sprintf("%dp",remap_round(var(size1),0,50000,1,5)))
z = eval(sprintf("%dp",remap_round(var(size1),0,50000,1,5)))
color = eval(sprintf("rdylbu-11-div-%d_a3",
remap_round(var(size1),0,250e6,1,11)))
</rule>
```

![Documentar%2086d1d/Untitled%2018.png](Documentar%2086d1d/Untitled%2018.png)

![Documentar%2086d1d/Untitled%2019.png](Documentar%2086d1d/Untitled%2019.png)

---

# 2.3 其它conf文件设置

### 2.3.1 ideogram.conf

```bash
<ideogram>
<spacing>
default = 2u 
#default = 10u
#default = 0.1r
# 直接设置default值可以设置默认的chromosomes的block之间的间距
# u设置的是绝对距离，1u就是conf文件中的1 chromosomes_unit
# 0.1r中的r表示相对距离，总的ideogram size的10%（如基因组总长185Mb，这个距离就是18.5Mb）

#<pairwise chr1>
#spacing = 5u
#spacing = 2r
#</pairwise>
# 调整chr1与隔壁两个block之间的间距

#<pairwise chr2 chr3>
#spacing = 15u
#</pairwise>
#<pairwise chr3 chr4>
#spacing = 25u
#</pairwise>
# 分别是调整chr2与chr3以及chr3与chr4之间的间距

# 如果在circos.conf文件中设置了chromosomes_break,并且想调整break的格式，加上下面这句并在[2.5 break.conf](../circos%20be3c7.md)文件中进行设置
  break   = 2u
# when relative, size of break is computed
# relative to default spacing
# break = 0.25r
<<include break.conf>>
	

</spacing>
<<incluce ../../etc/ideogram.conf>>

# the transparent colors of the band patterns. the value should be 1..N, where N is defined by auto_alpha_steps in the <image> block.
band_transparency* = 4

# In the karyotypes the labels of the chromosomes do not contain the species prefix. You can change the label string using label_format in ideogram.conf. This change the label to be the same as the chromosome name, which can be referenced with var(chr).
label_format = eval(sprintf(“%s”,var(chr)))

# ideogram的半径
radius*           = 0.5r
# 标签的半径
label_radius*     = 1.9r
# 边框的宽度
stroke_thickness* = 1
# 边框的颜色
stroke_color*     = vdgrey

</ideogram>
```

---

### 2.3.2 break.conf

```
***# break.conf***
	axis_break         = yes
  axis_break_style   = 2

# If an ideogram does not begin at its chromosome start or end, you can choose to place an axis break at the edge with
  axis_break_at_edge = yes

# There are two axis break styles, with style properties defined in break_style blocks. The size of the axis break is controlled by the break parameter and can be defined in absolute units (u) or relative to the default spacing (r).
  <break_style 1>
    stroke_color = black
    fill_color   = blue
    thickness    = 0.25r
    stroke_thickness = 2
  </break_style>

  <break_style 2>
    stroke_color     = black
    stroke_thickness = 3
    thickness        = 1.5r
  </break_style>

  </spacing>
</ideogram>

```

---

### 2.3.3 ticks.conf

ticks 会被分成多个组，每一组都是用<tick>和</tick>来定义的，同一组中的ticks会有相同的间隔，不同组可以创建多种ticks形式。以下code分别是先在1.1 20Mb为间隔的地方（20Mb, 40Mb, 60Mb..)画了ticks,之后是1.2.每10Mb为间隔的地方，但是之前已经画过的地方不会被替换,所以画的地方为10Mb，30Mb，50Mb....接下来是1.3.每隔2Mb的间隔区域。当然也可以加上栅格（grids）。

```
# 是否要画ticks
show_ticks = yes
show_tick_labels = yes
show_grid = yes

<ticks>
# 字体
tick_label_font = light
# 位置，这里是在ideogram外侧
radius = dims(ideogram,radius_outer) + 45p

label_offset = 5p
label_size = 8p
multiplier = 1e-6
color = black

# 1.1 absolute tick, every 20 Mb, with label and grid
<tick>
spacing = 20u
size = 12p
thickness = 2p
show_label = yes
label_size = 10p
format = %d
grid_start = 1r
grid_end = 1r+45p
grid_color = vdgrey
grid_thickness = 1p
grid = yes
</tick>

# 1.2 absolute tick, every 10 Mb, with grid
<tick>
spacing = 10u
size = 7p
thickness = 2p
show_label = no
grid_start = 1r
grid_end = 1r+45p
grid_color = grey
grid_thickness = 1p
grid = yes
</tick>

# 1.3 absolute tick, every 2 Mb
<tick>
spacing = 2u
size = 3p
thickness = 2p
show_label = no
</tick>

# 2.1 relative tick, every 2%
<tick>
# 位置信息，radius<1表示在ideogram的环内侧
radius = 0.75r
spacing_type = relative
rspacing = 0.02
size = 3p
thickness = 1p
show_label = no
</tick>

# 2.2 relative tick, every 10%, with label and grid
<tick>
radius = 0.75r
spacing_type = relative
rspacing = 0.10
size = 6p
thickness = 1p
show_label = yes
label_relative = yes
rmultiplier = 100
format = %d
suffix = %

grid_start = 0.5r
grid_end = 0.75r
grid_color = grey
grid_thickness = 1p
grid = yes

</tick>

</ticks>
```

---

# 3. 实战

## 3.1 根据位置来给links上色

本节以人类和小鼠的染色体为例，示例文件见circos-course/2/data/ucsc。其中人类和小鼠的染色体编号分别为hgN和mmN的格式（N为数字）。

### 3.1.1 提取染色体的编号

```
substr(var(chr),2)
# 转换为perl的语句类似于“foreach my $chr （@chrs) {$num=substr($chr,2)"
```

### 3.1.2 格式化染色体编号

```perl
sprintf(“chr%s_a4”,substr(var(chr),2))
# 在染色体编号前加上了前缀“chr”和后缀“_a4”。mm4会被转换为“chr4_a4”。
```

### 3.1.3 根据染色体编号设置颜色

```
<rule>
# 该判断一定成立
condition = 1
# 根据link的起始染色体设置颜色
color = eval(sprintf("chr%s_a4",substr(var(chr2),2)))
...
</rule>
```

![Documentar%2086d1d/Untitled%2020.png](Documentar%2086d1d/Untitled%2020.png)

当然也可以设置起始于某条特定的染色体上的link的颜色以及图层优先级

```
<rule>
# Instead of on(), you can use from() and to() to test whether a link starts and ends at a given
chromosome.
condition = on(mm11)
color = chr11_a3
z = 10
thickness = 2p
</rule>

```

![Documentar%2086d1d/Untitled%2021.png](Documentar%2086d1d/Untitled%2021.png)

进一步设置：

Below is an extension of the rule above. The condition of the second rule requires that a data
point have color chr11_a3 (i.e. as assigned by the previous rule, thus this condition tests whether the previous rule matched) as well as requires that the start and end position be within 20-50Mb on hs1 (all links have their start assigned to hs1). When more than one condition is defined, all need to match. The rule colors these links red (with transparency, _a3), sets their thickness to 3 pixels and their z value to a high value, to make these links drawn on top.

```
<rule>
condition = on(mm11)
color = chr11_a3
z = 10
thickness = 2p
# links that pass are tested by remaining rules
flow = continue
</rule>
<rule>
condition = var(color) eq “chr11_a3”
condition = var(start1) > 20e6
condition = var(end1) < 50e6
# 以上的condition也可以设置为：
#condition = on(m11)
#condition = var(start1) > 20e6 && var(end1) < 50e6
color = red_a3
z = 20
thickness = 3p
</rule>
```

![Documentar%2086d1d/Untitled%2022.png](Documentar%2086d1d/Untitled%2022.png)

### 3.1.4 把相邻的links转换为bundle

这个可以通过circos-tools里面自带的bundlelinks功能来实现，详细的说明可以参考[在线教程](http://www.circos.ca/documentation/tutorials/utilities/bundling_links/)。

```
$bundlelinks -links links.txt -max_gap_1 3e6 -min_bundle_membership 3 > bundles.txt$
# -max_gap: the distance between link starts or end (distance( start(L1), start(L2) ) <= MAX_GAP,distance( end(L1), end(L2) ) <= MAX_GAP)
# -map_gap_1: the distance between link starts
# -min_bundle_membership the minimum number of links required in a bundle for the bundle to be accepted
# -strict equivalent to -min_bundle_membership 2
# -min_bundle_size accept bundles based on the size of the links in the bundle. The minimum size parameter is applied independently to both ends of all links in a bundle.
# -min_bundle_identity filters bundles based on the bundle identity (identity = size(merged links) / extent(merged links))
```

```
<links>
<link>
ribbon = yes
file = ../data/bundles.txt
bezier_radius = 0r
radius = 0.85r
thickness = 0p
color = grey_a5
</link>
```

First, notice that ribbon=yes is set. This turns the links into ribbons—curved regions whose ends reflect the size of the ends of the link in the data file. Because bundles have ends that are typically much larger than the ends of the links from  which they are formed, using ribbons to draw bundles is a good idea.

![Documentar%2086d1d/Untitled%2023.png](Documentar%2086d1d/Untitled%2023.png)

Let’s focus the figure on links from mm5, mm8 and mm11. The first rule will change the color of all
links grey using the 9-color greys Brewer palette. The second rule will trigger only for links on
chromosomes mm5, mm8 or mm11 and change the color of these links to chr5, chr8 and
chr11.

```
<rules>
<rule>
condition = 1
color = eval(sprintf("greys-9-seq-%d_a5",remap_round(var(size1),0,1e6,2,9)))
z = eval(-var(size2))
flow = continue
</rule>
<rule>
condition = on(mm(5|8|11))
color = eval(sprintf("chr%s_a1",substr(var(chr2),2)))
z = 10
</rule>
</rules>
```

![Documentar%2086d1d/Untitled%2024.png](Documentar%2086d1d/Untitled%2024.png)

---

## 3.2 track automation

可以用track counters来自动生成和叠放tracks。可以参考网上教程：[automating tracks](http://circos.ca/documentation/tutorials/recipes/automating_tracks/)和[automating heatmaps](http://circos.ca/documentation/tutorials/recipes/automating_heatmaps/)

```
h0 = 0.75 # heatmap start
hs = 0.07 # heatmap step
hw = 0.06 # heatmap width
<plots>
<plot>
# 初始的counter（h）值设置为1
init_counter = h:1
type = heatmap
# 读入的文件名也跟h值有关
file = ../data/heatmap.counter(h).txt
# 相当于r1 = 0.81r # 0.75 + (1-1) * 0.07 + 0.06
r1 = eval(sprintf("%fr",conf(h0) + (counter(h)-1) * conf(hs) + conf(hw) ))
# 相当于r0 = 0.75r # 0.75 + (1-1) * 0.07
r0 = eval(sprintf("%fr",conf(h0) + (counter(h)-1) * conf(hs)))
# 相当于color = chr1_a3
color = eval(sprintf("chr%s_a3",counter(h)))
</plot>
<plot>
# counter值相当于前一个counter的增量为1
pre_increment_counter = h:1
type = heatmap
file = ../data/heatmap.counter(h).txt
r1 = eval(sprintf("%fr",conf(h0) + (counter(h)-1) * conf(hs) + conf(hw) ))
r0 = eval(sprintf("%fr",conf(h0) + (counter(h)-1) * conf(hs)))
color = eval(sprintf("chr%s_a3",counter(h)))
</plot>
<plot>
pre_increment_counter = h:1
type = heatmap
file = ../data/heatmap.counter(h).txt
r1 = eval(sprintf("%fr",conf(h0) + (counter(h)-1) * conf(hs) + conf(hw) ))
r0 = eval(sprintf("%fr",conf(h0) + (counter(h)-1) * conf(hs)))
color = eval(sprintf("chr%s_a3",counter(h)))
</plot>
</plots>
```

![Documentar%2086d1d/Untitled%2025.png](Documentar%2086d1d/Untitled%2025.png)

Counters are an advanced automation topic, so I’ve left them until the end of this session. Take
a look at the heatmap blocks below. The first one initializes a counter named h with value 1. The
value of this counter is available via counter(h). Subsequent blocks increment the value of the
counter before the block is processed using pre_increment_counter.

The position (r0, r1) of each track is defined to be a function of the counter and the global variables that set the start (h0), width (hw) and spacing (hs) of each track.

The result will be three heatmaps, each sourced from a different file (heatmap.1.txt, heatmap.2.txt, heatmap.3.txt), each at different r0/r1 and using a different color. I emphasize that the block definition for each heat map is identical (pre_increment_counter acts the same
as init_counter when called on an undefined counter name). The difference between the
blocks is due to the counter changing values.