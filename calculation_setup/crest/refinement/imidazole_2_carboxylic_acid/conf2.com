%nprocshared=8
%mem=6GB
%chk=conf2.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf2; regenerated study input

0 1
O -1.813623 1.208492 -0.028086
O -2.119748 -1.005405 0.026402
N 0.857941 1.015796 -0.018001
N 0.832235 -1.197723 0.021434
C 0.084547 -0.107686 0.002825
C 2.144913 0.610471 -0.012788
C 2.106834 -0.766102 0.011558
C -1.370660 -0.050325 0.002350
H 0.510157 1.963236 -0.033600
H 2.970154 1.290574 -0.025812
H 2.918815 -1.463673 0.022008
H -2.789430 1.251695 -0.030594

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf2.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

