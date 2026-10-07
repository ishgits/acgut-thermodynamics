%nprocshared=8
%mem=6GB
%chk=calc_f25e3854ff8505f5.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_f25e3854ff8505f5; regenerated study input

0 1
N -0.758596 0.794598 -0.000008
N 0.140102 -1.223962 0.000017
C 0.602525 0.982401 0.000055
C -0.987057 -0.539450 0.000027
C 1.141543 -0.272013 -0.000063
H -1.462016 1.517070 -0.000290
H 1.046842 1.961204 0.000147
H -1.980287 -0.954128 0.000039
H 2.182849 -0.544234 -0.000077

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_f25e3854ff8505f5.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

