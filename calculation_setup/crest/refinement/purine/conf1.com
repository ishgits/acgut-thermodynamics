%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
N -1.559462 1.089680 -0.000022
N -1.524969 -1.147513 -0.000024
N 0.886812 -1.401462 0.000019
N 2.091003 0.652887 -0.000022
C -0.240475 0.735278 0.000021
C -0.238500 -0.694733 0.000024
C 0.982709 1.378092 0.000018
C -2.260278 -0.071515 0.000009
C 1.985457 -0.671893 -0.000020
H -1.936960 2.025028 -0.000072
H 1.104603 2.452381 0.000094
H -3.332098 -0.099121 0.000054
H 2.918270 -1.218282 -0.000002

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

