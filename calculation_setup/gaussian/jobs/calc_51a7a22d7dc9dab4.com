%nprocshared=8
%mem=6GB
%chk=calc_51a7a22d7dc9dab4.chk
# opt b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water)

calc_51a7a22d7dc9dab4; regenerated study input

0 1
O -2.246464 -1.038952 -0.000000
O 2.309857 -0.999872 0.000001
N 0.037944 -0.983031 -0.000000
N -1.173329 0.977638 0.000001
C -1.202703 -0.393697 0.000000
C 1.271461 -0.333895 -0.000000
C -0.014139 1.699748 0.000000
C 1.194084 1.102295 -0.000001
H 0.049114 -1.995963 -0.000000
H -2.067081 1.449692 0.000001
H 2.105626 1.676137 -0.000001
H -0.139322 2.771768 0.000000

--Link1--
%nprocshared=8
%mem=6GB
%chk=calc_51a7a22d7dc9dab4.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(smd,solvent=water) temperature=298 geom=allchk guess=read

