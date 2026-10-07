%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O 1.779662 1.867492 0.000034
N -1.295289 1.351722 -0.000039
N 1.906797 -0.412367 0.000075
N -2.075210 -0.733854 -0.000020
N 0.069039 -1.883146 -0.000155
C -0.204084 0.530346 -0.000054
C -0.718926 -0.778011 -0.000040
C 1.199821 0.783627 -0.000174
C -2.375230 0.540219 0.000155
C 1.332282 -1.644980 0.000113
H -1.278033 2.360543 -0.000093
H 2.918117 -0.328331 0.000366
H -3.382099 0.906032 0.000539
H 2.002884 -2.490600 0.000568

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

