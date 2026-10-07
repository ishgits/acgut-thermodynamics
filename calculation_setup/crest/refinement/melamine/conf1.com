%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
N -0.757454 -1.142201 0.005342
N -0.610641 1.227146 -0.016418
N 1.367810 -0.084774 0.006344
N -2.638099 0.163748 0.019772
N 1.177466 -2.366693 -0.023070
N 1.461022 2.202627 0.019041
C -1.298825 0.080547 -0.004371
C 0.579584 -1.164925 0.006065
C 0.718962 1.084625 -0.004175
H -3.063374 1.067243 -0.077056
H -3.172593 -0.680907 -0.066579
H 0.607460 -3.185961 0.080227
H 2.174410 -2.406798 0.083296
H 0.998870 3.087794 -0.077746
H 2.457110 2.117735 -0.065598

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

