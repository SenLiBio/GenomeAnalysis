my @phys=glob("*phy");

foreach my $phy (@phys) {

	(my $prefix=$phy)=~s/\.phy//;
	open MA,">","$prefix.ma.ctl";
	print MA "seqfile = $phy;
treefile = speciestree.tree
outfile = $phy.ma
noisy = 3
verbose = 0
runmode = 0
seqtype = 1 
CodonFreq = 2
clock = 0 
model = 2
NSsites = 2
icode = 0
Mgene = 0
fix_kappa = 0 
kappa = 1.234567
fix_omega = 0
omega = 1.414
fix_alpha = 1
alpha = 0
ncatG = 3
getSE = 1
Small_Diff = .5e-6
*    cleandata = 0
*        ndata = 10
*  fix_blength = 0
method = 0";

	open M0,">","$prefix.m0.ctl";
	print M0 "seqfile = $phy
treefile = speciestree.tree
outfile = $phy.m0
noisy = 3
verbose = 0
runmode = 0
seqtype = 1 
CodonFreq = 2
clock = 0 
model = 2
NSsites = 2
icode = 0
Mgene = 0
fix_kappa = 0 
kappa = 1.234567
fix_omega = 1
omega = 1
fix_alpha = 1
alpha = 0
ncatG = 3
getSE = 1
Small_Diff = .5e-6
*    cleandata = 0
*        ndata = 10
*  fix_blength = 0
method = 0";

open SH0,">>","codeml.m0.sh";
open SHA,">>","codeml.ma.sh";
print SH0 "codeml $prefix.m0.ctl\n";
print SHA "codeml $prefix.ma.ctl\n";
}

`split -l 500 -d -a 2 codeml.m0.sh codeml.m0.sh_`;

my @m0shs=glob("codeml.m0.sh_*");
foreach my $m0sh (@m0shs) {
        `nohup bash $m0sh &`;
}

`split -l 500 -d -a 2 codeml.ma.sh codeml.ma.sh_`;

my @mashs=glob("codeml.ma.sh_*");
foreach my $mash (@mashs) {
        `nohup bash $mash &`;
}