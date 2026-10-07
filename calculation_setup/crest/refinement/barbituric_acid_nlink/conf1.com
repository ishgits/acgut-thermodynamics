%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O 2.370045 -1.238179 -0.070131
O -2.370203 -1.237925 -0.070036
O 0.000159 2.597418 0.001267
N 1.163111 0.654529 0.006124
N -1.163037 0.654660 0.005944
C -0.000076 -1.489216 0.120601
C 1.281461 -0.713088 0.003716
C -1.281541 -0.712954 0.003708
C 0.000076 1.387196 0.008900
H -0.000113 -2.300456 -0.613320
H -0.000090 -1.993519 1.094288
H 2.025508 1.191852 -0.037601
H -2.025388 1.192060 -0.037848

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

