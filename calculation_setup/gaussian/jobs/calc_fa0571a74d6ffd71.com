%nprocshared=8
%mem=6GB
%chk=calc_fa0571a74d6ffd71.chk
# opt wb97xd/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

calc_fa0571a74d6ffd71; regenerated study input

0 1
O 0.000000 0.116984 -0.000000
H 0.759192 -0.467938 0.000000
H -0.759192 -0.467938 0.000000

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_fa0571a74d6ffd71.chk
# freq wb97xd/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read

