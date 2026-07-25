#!/usr/bin/env python3
"""Regenerate the test fixtures. They are generated, not shipped.

    python3 mkfixtures.py && python3 vh.py && python3 VALIDATE_ALL.py
"""
import numpy as np
import harvest

R = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
H = np.zeros((3, 2, 2), complex)
H[1, 0, 1] = 1.0
H[2, 1, 0] = 1.0
harvest.write_hr('dimer_hr.dat', R, H, np.ones(3), 'dimerized chain')
print("wrote dimer_hr.dat")

import mktest  # noqa: F401  (writes cb_M0.0_hr.dat, cb_M0.9_hr.dat on import)
