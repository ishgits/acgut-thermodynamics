from pathlib import Path
import pytest


@pytest.fixture
def root():
    return Path(__file__).resolve().parents[1]


@pytest.fixture
def water_log(tmp_path):
    # Entirely synthetic compact output, not original Gaussian output bytes.
    orient = """ Standard orientation:
 --------------------
 Center Atomic Atomic Coordinates
 Number Number Type X Y Z
 --------------------
 1 8 0 0.000000 0.000000 0.000000
 2 1 0 0.957200 0.000000 0.000000
 3 1 0 -0.239987 0.926627 0.000000
 --------------------
"""
    opt = "# opt b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)\n --------------------\n"
    freq = "# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 geom=allchk guess=read\n --------------------\n"
    text = opt + orient + " Optimization completed.\n Normal termination of Gaussian\n" + freq + orient + """ Stoichiometry H2O
 Charge = 0 Multiplicity = 1
 SCF Done: E(RB3LYP) = -76.400000000 A.U.
 Frequencies -- 1600.0 3500.0 3600.0
 Temperature 298.000 Kelvin. Pressure 1.00000 Atm.
 Thermal correction to Gibbs Free Energy= 0.005000
 Sum of electronic and thermal Free Energies= -76.395000
 Normal termination of Gaussian
"""
    path = tmp_path / "synthetic_water.txt"; path.write_text(text)
    return path
