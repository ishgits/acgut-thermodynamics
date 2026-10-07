%nprocshared=8
%mem=6GB
%chk=conf1.chk
# opt=(tight,calcfc,maxcycles=200) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

conf1; regenerated study input

0 1
O 2.368913 -1.238458 -0.078734
O -2.368348 -1.239381 -0.078805
O -0.000586 2.597478 0.001867
N 1.162999 0.654857 0.006253
N -1.163252 0.654368 0.006285
C 0.000295 -1.487693 0.136323
C 1.281229 -0.712875 0.003955
C -1.280918 -0.713386 0.004093
C -0.000295 1.387257 0.010059
H 0.000313 -2.315070 -0.579153
H 0.000278 -1.968799 1.121801
H 2.025071 1.192238 -0.043004
H -2.025545 1.191383 -0.043047

--Link1--
%nprocshared=8
%mem=6GB
%chk=conf1.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read

