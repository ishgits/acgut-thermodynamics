%nprocshared=8
%mem=6GB
%chk=calc_67ddb66290e3abfb.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

calc_67ddb66290e3abfb; regenerated study input

0 1
O -2.052096 -1.076665 0.000002
O -1.956931 1.173847 -0.000002
N 0.864207 1.041702 -0.000002
N 0.766230 -1.170513 -0.000001
C 0.055224 -0.053847 -0.000000
C 2.149616 0.602821 0.000001
C 2.069815 -0.773406 -0.000001
C -1.400734 0.094397 0.000002
H 0.549700 1.999764 0.000014
H 2.987587 1.275822 0.000003
H 2.881470 -1.480566 0.000001
H -3.003130 -0.890580 -0.000003

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_67ddb66290e3abfb.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read

