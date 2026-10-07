%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
P 0.075018 0.030487 -0.059344
O 0.309055 -1.238098 0.872803
O -0.602071 1.108374 0.890577
O -1.098592 -0.399069 -1.003148
O 1.285955 0.493523 -0.754311
H 1.238723 -1.434148 1.013535
H 0.016433 1.753050 1.244066
H -1.883363 -0.695931 -0.527972

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

