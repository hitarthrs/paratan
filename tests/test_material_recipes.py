"""Material recipes must import without accidental void-fraction warnings."""
import os
from pathlib import Path
import subprocess
import sys
import unittest


class MaterialRecipeTests(unittest.TestCase):
    def test_import_without_fraction_warning_and_rebco_density(self):
        code = '''
import warnings
warnings.filterwarnings('error', message='Warning: sum of fractions do not add to 1.*')
from src.paratan.materials import material as m
fractions = [0.1312508, 0.05250033, 0.02625016, 0.000525,
             0.00131251, 0.00013125, 0.000525, 0.7875049]
components = [m.copper, m.silver, m.rebco, m.lamno3, m.MgO,
              m.yttrium_oxide, m.alumina, m.hastelloy]
expected = sum(f * c.get_mass_density() for f, c in zip(fractions, components))/sum(fractions)
assert abs(m.rebco_tape.get_mass_density()-expected) < 1e-10
'''
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        result = subprocess.run([sys.executable, '-c', code],
                                cwd=Path(__file__).resolve().parents[1],
                                env=environment, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)


if __name__ == '__main__':
    unittest.main()
