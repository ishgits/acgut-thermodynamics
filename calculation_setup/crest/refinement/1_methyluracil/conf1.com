%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O 1.464034 1.754196 -0.000594
O -2.823707 0.304918 0.001609
N 1.120424 -0.501730 -0.001992
N -0.663514 0.978716 -0.000322
C 0.696292 0.805576 -0.000870
C 0.231061 -1.527513 -0.001244
C 2.543182 -0.784382 0.002911
C -1.104514 -1.343301 -0.000780
C -1.639909 -0.006849 0.000314
H 0.667388 -2.515119 -0.001718
H 2.812089 -1.326043 0.908169
H 2.808180 -1.374897 -0.872025
H 3.074205 0.164757 -0.023203
H -0.984983 1.941853 0.001561
H -1.796106 -2.165977 -0.000680

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

