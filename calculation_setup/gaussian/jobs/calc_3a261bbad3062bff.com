%nprocshared=8
%mem=6GB
%chk=calc_3a261bbad3062bff.chk
# opt wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_3a261bbad3062bff; regenerated study input

0 1
N -0.759298 0.788617 -0.000070
N 0.144192 -1.220208 0.000012
C 0.596108 0.981484 0.000118
C -0.980616 -0.540500 0.000075
C 1.138729 -0.266687 -0.000096
H -1.465627 1.505915 -0.000374
H 1.037047 1.962194 0.000238
H -1.972921 -0.958882 0.000100
H 2.181911 -0.533873 -0.000139

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_3a261bbad3062bff.chk
# freq wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

