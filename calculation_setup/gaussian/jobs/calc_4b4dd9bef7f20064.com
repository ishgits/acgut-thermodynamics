%nprocshared=8
%mem=6GB
%chk=calc_4b4dd9bef7f20064.chk
# opt wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_4b4dd9bef7f20064; regenerated study input

0 1
O 2.229778 -1.132639 0.000001
N -0.011544 -0.915728 -0.000003
N -2.316389 -0.846243 0.000003
N 1.375906 0.986446 0.000001
C -1.156380 -0.197889 -0.000001
C -1.030763 1.189959 -0.000001
C 1.263647 -0.357108 -0.000000
C 0.247619 1.697767 0.000000
H -0.056156 -1.926236 -0.000003
H -1.903701 1.822323 -0.000001
H 0.375164 2.774334 -0.000001
H -3.176163 -0.325528 -0.000004
H -2.357923 -1.851481 -0.000001

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_4b4dd9bef7f20064.chk
# freq wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

