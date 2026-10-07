%nprocshared=8
%mem=6GB
%chk=calc_0f12ea6b595a93ec.chk
# opt wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_0f12ea6b595a93ec; regenerated study input

0 1
O -2.236238 -1.036072 -0.000006
O 2.299670 -0.997192 0.000010
N 0.038354 -0.981047 -0.000003
N -1.170753 0.974239 0.000009
C -1.198075 -0.391351 -0.000000
C 1.265915 -0.335272 -0.000003
C -0.013583 1.694910 0.000001
C 1.190886 1.101429 -0.000012
H 0.048027 -1.992562 -0.000002
H -2.063057 1.445608 0.000025
H 2.102624 1.675250 -0.000023
H -0.139111 2.767172 0.000003

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_0f12ea6b595a93ec.chk
# freq wb97xd/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

