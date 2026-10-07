%nprocshared=8
%mem=6GB
%chk=calc_6e9313e7d80730f3.chk
# opt wb97xd/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

calc_6e9313e7d80730f3; regenerated study input

0 1
N -0.748422 0.798838 -0.000065
N 0.128267 -1.220102 0.000009
C 0.608803 0.975082 0.000103
C -0.987724 -0.530026 0.000068
C 1.134559 -0.282502 -0.000086
H -1.443186 1.525149 -0.000323
H 1.061938 1.949805 0.000217
H -1.985516 -0.936500 0.000098
H 2.174025 -0.564938 -0.000112

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_6e9313e7d80730f3.chk
# freq wb97xd/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read

