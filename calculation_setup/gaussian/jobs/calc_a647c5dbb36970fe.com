%nprocshared=8
%mem=6GB
%chk=calc_a647c5dbb36970fe.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_a647c5dbb36970fe; regenerated study input

0 1
P -0.036688 -0.043655 0.106731
O -1.064133 -0.009784 -1.106851
O 0.496800 1.447595 0.136138
O 1.240177 -0.859384 -0.350013
O -0.577336 -0.586973 1.368442
H -1.652791 -0.777680 -1.146955
H -0.166051 2.111497 0.376146
H 1.605095 -0.610623 -1.211898

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_a647c5dbb36970fe.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

