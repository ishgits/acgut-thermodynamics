%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

dft_seed; regenerated study input

0 1
O 2.378550 -1.233814 -0.000014
O -2.378548 -1.233817 -0.000013
O -0.000001 2.622644 0.000001
N 1.170350 0.673571 0.000003
N -1.170351 0.673570 -0.000000
C 0.000001 -1.481905 0.000020
C 1.290211 -0.699608 -0.000000
C -1.290211 -0.699609 -0.000001
C -0.000001 1.410973 0.000002
H 0.000001 -2.138359 -0.871061
H 0.000000 -2.138303 0.871142
H 2.025639 1.213730 -0.000005
H -2.025640 1.213728 -0.000008

--Link1--
%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

