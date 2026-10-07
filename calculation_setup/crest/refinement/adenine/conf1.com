%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
N -1.240865 -1.391523 0.000548
N -2.098894 0.668971 -0.000048
N 0.004728 1.883223 -0.000212
N 1.938206 0.490749 0.001689
N 1.898975 -1.807309 -0.010699
C -0.171499 -0.530475 0.001007
C -0.736842 0.765352 0.000570
C 1.227745 -0.645261 0.003252
C -2.346870 -0.612003 -0.000926
C 1.299050 1.652311 -0.001048
H -1.221581 -2.398793 -0.006686
H -3.339816 -1.015300 -0.002129
H 1.934511 2.528635 -0.002105
H 2.902224 -1.766105 0.029537
H 1.426500 -2.688397 0.068551

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

