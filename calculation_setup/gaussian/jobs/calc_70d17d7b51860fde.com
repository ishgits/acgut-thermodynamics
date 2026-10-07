%nprocshared=8
%mem=6GB
%chk=calc_70d17d7b51860fde.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

calc_70d17d7b51860fde; regenerated study input

0 1
P -0.055664 0.039098 -0.109695
O -0.957850 -0.213613 1.181685
O 0.773345 -1.313290 -0.197329
O 1.065957 1.085814 0.292517
O -0.784517 0.469716 -1.310652
H -1.723592 0.371638 1.239133
H 0.253504 -2.088963 -0.442691
H 1.529569 0.901851 1.119218

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_70d17d7b51860fde.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read

