%nprocshared=8
%mem=6GB
%chk=calc_c2c9d57073565a33.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

calc_c2c9d57073565a33; regenerated study input

0 1
O 0.000000 0.117577 0.000000
H 0.762056 -0.470309 -0.000000
H -0.762056 -0.470309 -0.000000

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_c2c9d57073565a33.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read

