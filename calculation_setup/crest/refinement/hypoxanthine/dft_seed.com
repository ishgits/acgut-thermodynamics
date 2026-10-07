%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

dft_seed; regenerated study input

0 1
O 1.583598 2.040230 0.000003
N -1.421096 1.215089 -0.000002
N 1.932755 -0.234339 -0.000002
N -1.996803 -0.940684 0.000002
N 0.246393 -1.882811 0.000001
C -0.251038 0.500166 -0.000001
C -0.631982 -0.837877 -0.000000
C 1.114622 0.905509 -0.000001
C -2.427312 0.306996 0.000000
C 1.492587 -1.526179 0.000000
H -1.518386 2.218172 -0.000004
H 2.929422 -0.064075 -0.000002
H -3.461419 0.608013 -0.000000
H 2.271589 -2.276424 -0.000000

--Link1--
%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

