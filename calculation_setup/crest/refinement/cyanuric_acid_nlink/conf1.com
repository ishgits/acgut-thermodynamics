%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O 1.211484 2.324524 -0.000046
O -2.618943 -0.113285 0.000003
O 1.407529 -2.211381 -0.000064
N -0.716562 1.125489 0.000013
N 1.332878 0.057817 0.000035
N -0.616357 -1.183184 0.000026
C 0.653017 1.253612 0.000050
C -1.412278 -0.061242 -0.000068
C 0.759208 -1.192352 0.000079
H -1.262448 1.982304 0.000002
H 2.347770 0.102458 -0.000054
H -1.085203 -2.084429 0.000004

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

