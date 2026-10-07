%nprocshared=8
%mem=6GB
%chk=calc_1aa4f60ffba16a6a.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_1aa4f60ffba16a6a; regenerated study input

0 1
O 2.375601 -1.233849 -0.001461
O -2.375607 -1.233840 -0.001195
O 0.000004 2.619549 0.000234
N 1.170640 0.673142 0.000250
N -1.170636 0.673145 -0.000242
C -0.000002 -1.477560 0.002014
C 1.282022 -0.696062 -0.000035
C -1.282025 -0.696059 -0.000045
C 0.000002 1.402313 0.000094
H -0.000012 -2.135669 -0.868707
H 0.000006 -2.130000 0.877093
H 2.027326 1.215494 -0.000382
H -2.027319 1.215502 -0.000842

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_1aa4f60ffba16a6a.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

