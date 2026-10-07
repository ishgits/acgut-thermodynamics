%nprocshared=8
%mem=6GB
%chk=calc_9245a11b85b2f147.chk
# opt wb97xd/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

calc_9245a11b85b2f147; regenerated study input

0 1
P -0.056140 0.033509 -0.109302
O -0.931823 -0.166946 1.198342
O 0.777039 -1.307962 -0.152669
O 1.055232 1.100817 0.227554
O -0.807009 0.401557 -1.310365
H -1.696078 0.414708 1.241540
H 0.260526 -2.092915 -0.355612
H 1.530144 0.955842 1.050708

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_9245a11b85b2f147.chk
# freq wb97xd/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read

