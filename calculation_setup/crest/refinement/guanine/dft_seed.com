%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

dft_seed; regenerated study input

0 1
O -0.079332 2.654923 0.000792
N 2.195924 0.511102 -0.000163
N -1.457550 0.820640 0.000579
N -0.741120 -1.431253 0.005050
N 1.673034 -1.658145 0.001646
N -2.996502 -0.922951 -0.057966
C 0.820612 0.436981 0.000295
C 0.517509 -0.920680 0.000971
C -0.191078 1.430186 0.001025
C -1.694655 -0.534401 -0.000835
C 2.646495 -0.762873 0.000573
H 2.758126 1.346744 -0.001623
H -2.239259 1.460852 -0.016195
H 3.698936 -0.992248 0.000153
H -3.170797 -1.897174 0.124424
H -3.722142 -0.288571 0.230714

--Link1--
%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

