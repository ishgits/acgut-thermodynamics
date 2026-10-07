%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O 2.191409 -1.076296 -0.000010
N -0.059600 -0.920873 0.000010
N -2.377021 -0.882888 -0.000014
N 1.283158 1.013629 -0.000001
C -1.217434 -0.227509 0.000015
C -1.113406 1.167610 -0.000002
C 1.219914 -0.323480 0.000006
C 0.158082 1.706620 0.000002
H -0.063646 -1.934753 -0.000003
H -1.989464 1.788954 -0.000084
H 0.296849 2.779742 0.000022
H -3.238981 -0.365255 0.000020
H -2.406658 -1.888297 0.000028

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

