"""Independent 60-decimal evaluation of the documented equations.

Standard-library Decimal only: no production imports. Run explicitly to inspect
fixture derivation; never invoked by tests or reference-regeneration workflows.
Constants have SI units except light speed (cm/s) and frequencies (cm^-1).
"""
from decimal import Decimal as D, localcontext
import json
from pathlib import Path


def derive():
    with localcontext() as ctx:
        ctx.prec = 60
        h, k, c, R, T = map(D, ['6.62607015e-34', '1.380649e-23', '2.99792458e10', '8.31446261815324', '298'])
        pi = D('3.14159265358979323846264338327950288419716939937510582097494459')
        bav, JH = D('1e-44'), D('4184') * D('627.5094740631')
        def ho(v):
            x = h * c * v / (k * T)
            return R * (x / (x.exp() - 1) - (1 - (-x).exp()).ln())
        def mode(v):
            mu = h / (8 * pi**2 * v * c)
            effective = 1 / (1 / mu + 1 / bav)
            rotor = R / 2 * (1 + (8 * pi**3 * effective * k * T / h**2).ln())
            w = v**4 / (v**4 + D(100)**4)
            return {'harmonic_entropy_j_mol_k': ho(v), 'damping': w, 'rotor_entropy_j_mol_k': rotor,
                    'floor_50': T * (ho(v) - ho(max(v, D(50)))) / JH,
                    'floor_100': T * (ho(v) - ho(max(v, D(100)))) / JH,
                    'grimme_qrrho_100': T * (1-w) * (ho(v)-rotor) / JH}
        frequencies = ['15', '49.999', '50', '50.001', '99.999', '100', '100.001', '3000']
        rows = [{'frequency_cm1': f, **{key: str(value) for key,value in mode(D(f)).items()}} for f in frequencies]
        spectrum = ['15', '100', '3000']
        expected = {key: str(sum(mode(D(f))[key] for f in spectrum)) for key in ['floor_50','floor_100','grimme_qrrho_100']}
        return {'provenance': 'Independent Decimal precision=60 evaluation; equations in docs/EQUATIONS.md. Algebra uses reciprocal reduced inertia and v^4/(v^4+cutoff^4). No production functions.',
                'constants': dict(h=str(h), k=str(k), c_cm_s=str(c), R=str(R), temperature_k=str(T), bav_kg_m2=str(bav), joule_mol_per_hartree=str(JH)),
                'absolute_tolerance_hartree': 1e-15, 'absolute_tolerance_entropy_j_mol_k': 1e-12,
                'modes': rows, 'spectrum_cm1': spectrum, 'spectrum_shifts_hartree': expected}

if __name__ == '__main__':
    Path(__file__).with_name('entropy.json').write_text(json.dumps(derive(), indent=2) + '\n')
