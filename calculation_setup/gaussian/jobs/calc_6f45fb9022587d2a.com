%nprocshared=8
%mem=6GB
%chk=calc_6f45fb9022587d2a.chk
# opt wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_6f45fb9022587d2a; regenerated study input

0 1
O -0.000000 0.117347 0.000000
H 0.759484 -0.469388 -0.000000
H -0.759484 -0.469388 0.000000

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_6f45fb9022587d2a.chk
# freq wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

