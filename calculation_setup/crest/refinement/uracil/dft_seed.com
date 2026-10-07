%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

dft_seed; regenerated study input

0 1
O -2.249949 -1.035424 -0.000000
O 2.309317 -1.001811 0.000000
N 0.036577 -0.983160 -0.000000
N -1.173095 0.980360 0.000000
C -1.210636 -0.397940 0.000000
C 1.279552 -0.338956 -0.000000
C -0.011183 1.701274 0.000000
C 1.196157 1.104024 -0.000000
H 0.048355 -1.994094 -0.000000
H -2.065428 1.449943 0.000000
H 2.107455 1.677575 -0.000000
H -0.133034 2.773640 0.000000

--Link1--
%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

