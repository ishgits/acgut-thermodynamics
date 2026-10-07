%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O -2.002088 -1.086535 0.000103
O -1.934299 1.149104 -0.000103
N 0.878095 1.067480 0.000153
N 0.815876 -1.145257 -0.000051
C 0.087454 -0.042830 0.000015
C 2.158555 0.641442 -0.000022
C 2.097408 -0.734703 -0.000083
C -1.365352 0.067307 -0.000066
H 0.537769 2.017467 0.000287
H 2.994118 1.308926 -0.000110
H 2.898660 -1.444715 0.000227
H -2.973847 -0.974404 0.000034

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

