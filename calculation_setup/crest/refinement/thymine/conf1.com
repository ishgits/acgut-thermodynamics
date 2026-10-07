%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O -1.494335 1.793508 0.000129
O 2.785602 0.299049 0.000220
N 0.625033 1.013091 -0.000203
N 1.113285 -1.235518 -0.000072
C -1.169875 -0.567474 -0.000053
C -0.748478 0.824028 -0.000127
C -0.220882 -1.526282 -0.000088
C -2.635281 -0.854292 0.000126
C 1.592675 0.041158 -0.000098
H -0.473977 -2.574793 0.000635
H 0.952012 1.974243 -0.000013
H 1.800485 -1.979098 0.000310
H -3.098251 -0.398102 0.873016
H -3.097360 -0.404857 -0.876769
H -2.820765 -1.922746 0.003957

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

