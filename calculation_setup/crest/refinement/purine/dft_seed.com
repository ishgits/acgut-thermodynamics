%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

dft_seed; regenerated study input

0 1
N -1.573815 1.024681 -0.000002
N -1.493691 -1.216601 0.000001
N 0.927455 -1.408306 0.000001
N 2.097824 0.685615 -0.000001
C -0.239508 0.694914 0.000002
C -0.212714 -0.715944 0.000001
C 0.959865 1.382222 0.000001
C -2.261575 -0.149263 -0.000001
C 2.023868 -0.655773 -0.000001
H -1.978470 1.947428 0.000007
H 1.022836 2.463252 0.000001
H -3.339242 -0.164073 -0.000004
H 2.970858 -1.181272 -0.000005

--Link1--
%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

