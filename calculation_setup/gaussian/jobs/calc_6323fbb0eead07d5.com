%nprocshared=8
%mem=6GB
%chk=calc_6323fbb0eead07d5.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

calc_6323fbb0eead07d5; regenerated study input

0 1
O 1.926025 1.097286 0.000002
O 2.082682 -1.151024 -0.000001
N -0.788360 1.033826 -0.000001
N -0.847693 -1.179609 0.000001
C -0.057655 -0.118151 -0.000001
C -2.101833 0.686671 -0.000000
C -2.119128 -0.691509 0.000001
C 1.407573 -0.148120 -0.000000
H -0.414569 1.969944 -0.000004
H -2.890031 1.417452 -0.000001
H -2.978679 -1.339516 0.000002
H 2.892241 1.029153 0.000001

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_6323fbb0eead07d5.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read

