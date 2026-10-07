%nprocshared=8
%mem=6GB
%chk=calc_365dba8199d0b0ba.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_365dba8199d0b0ba; regenerated study input

0 1
O -0.000000 0.000000 0.117945
H 0.000000 -0.762275 -0.471782
H -0.000000 0.762275 -0.471782

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_365dba8199d0b0ba.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

