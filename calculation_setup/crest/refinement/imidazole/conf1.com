%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
N -1.076928 -0.322684 0.000114
N 1.125950 -0.394042 0.000105
C -0.643258 0.971415 -0.000048
C 0.026925 -1.100897 -0.000156
C 0.719971 0.906006 -0.000082
H -2.036219 -0.629440 0.000244
H -1.313298 1.804758 -0.000217
H -0.002750 -2.171320 -0.000119
H 1.436071 1.702735 0.000459

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

