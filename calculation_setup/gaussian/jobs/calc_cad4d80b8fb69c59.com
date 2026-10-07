%nprocshared=8
%mem=6GB
%chk=calc_cad4d80b8fb69c59.chk
# opt wb97xd/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

calc_cad4d80b8fb69c59; regenerated study input

0 1
O 2.369693 -1.229719 -0.001048
O -2.369694 -1.229718 -0.000925
O 0.000000 2.612701 -0.000059
N 1.167657 0.672197 0.000254
N -1.167657 0.672197 0.000039
C -0.000000 -1.478038 0.001339
C 1.286107 -0.697046 0.000011
C -1.286107 -0.697047 -0.000006
C -0.000000 1.406131 0.000229
H 0.000001 -2.134028 -0.869278
H -0.000001 -2.130255 0.874859
H 2.021364 1.211707 -0.000293
H -2.021364 1.211707 -0.000519

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_cad4d80b8fb69c59.chk
# freq wb97xd/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read

