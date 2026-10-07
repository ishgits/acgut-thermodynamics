%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

dft_seed; regenerated study input

0 1
O -0.195410 2.645199 0.000000
O -2.974731 -0.987932 -0.000000
N -0.729232 -1.401586 0.000001
N 2.146573 0.586237 -0.000000
N -1.533900 0.785253 -0.000000
N 1.700728 -1.597551 0.000000
C 0.773302 0.465058 0.000000
C 0.543994 -0.893455 0.000000
C -0.283547 1.424563 0.000000
C -1.826427 -0.572935 -0.000000
C 2.647386 -0.664046 -0.000000
H -0.894828 -2.396712 0.000001
H 2.674736 1.444396 -0.000001
H -2.342141 1.391969 -0.000001
H 3.705931 -0.859378 -0.000001

--Link1--
%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

