%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

dft_seed; regenerated study input

0 1
O 2.371937 1.144109 -0.000004
O -2.176942 1.481919 -0.000001
O -0.195079 -2.626196 -0.000001
N 0.099075 1.332806 0.000002
N 1.104762 -0.752065 0.000001
N -1.203644 -0.580675 0.000000
C 1.281109 0.617903 0.000005
C -1.175689 0.800648 0.000000
C -0.105510 -1.418445 0.000001
H 0.173710 2.339712 0.000000
H 1.939325 -1.320369 -0.000001
H -2.113176 -1.019094 -0.000002

--Link1--
%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

