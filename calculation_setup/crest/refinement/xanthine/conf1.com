%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O -0.195890 2.637233 -0.000021
O -2.933049 -0.987699 0.000097
N -0.705814 -1.412607 -0.000076
N 2.162251 0.591906 0.000102
N -1.506672 0.775171 0.000031
N 1.722086 -1.587258 -0.000028
C 0.795626 0.468193 -0.000048
C 0.558332 -0.905295 -0.000149
C -0.273697 1.416684 -0.000117
C -1.785451 -0.577227 0.000031
C 2.655736 -0.657451 0.000081
H -0.861415 -2.413321 -0.000032
H 2.673155 1.462067 0.000261
H -2.323306 1.378880 0.000278
H 3.702784 -0.881716 0.000286

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

