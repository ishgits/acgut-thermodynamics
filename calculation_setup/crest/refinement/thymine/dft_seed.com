%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

dft_seed; regenerated study input

0 1
O -1.369137 1.894020 -0.000001
O 2.869659 0.222491 -0.000000
N 0.726844 1.024956 0.000001
N 1.111448 -1.239189 -0.000001
C -1.153430 -0.478367 0.000000
C -0.663765 0.892138 -0.000000
C -0.241145 -1.471936 -0.000000
C -2.631724 -0.717522 0.000001
C 1.665764 0.017108 0.000000
H -0.521915 -2.514898 -0.000000
H 1.091040 1.967944 0.000002
H 1.758107 -2.012360 0.000001
H -3.101314 -0.269178 0.876630
H -3.101314 -0.269181 -0.876631
H -2.851019 -1.783311 0.000002

--Link1--
%nprocshared=8
%mem=6GB
%chk=dft_seed.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

