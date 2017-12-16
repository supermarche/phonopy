[![Version Badge](https://anaconda.org/atztogo/phonopy/badges/version.svg)](https://anaconda.org/atztogo/phonopy)
[![Downloads Badge](https://anaconda.org/atztogo/phonopy/badges/downloads.svg)](https://anaconda.org/atztogo/phonopy)
[![PyPI](https://img.shields.io/pypi/dm/phonopy.svg?maxAge=2592000)](https://pypi.python.org/pypi/phonopy)
[![Build Status](https://travis-ci.org/atztogo/phonopy.svg?branch=master)](https://travis-ci.org/atztogo/phonopy)

phonopy | FPLO interface
========================

Phonon code. For general details, see http://atztogo.github.io/phonopy/

For the time being, the fplo/fedit commands have to be hard coded in interfaces/fplo.py.
Only spacegroup 1 is supported, i.e. no symmetry information is used.


FPLO computational setup
------------------------

The input has to be an =.in file with the structure represented in spacegroup 1.

Be aware that the default basis set of FPLO has to be increased (see
[here](https://molmod.ugent.be/sites/default/files/deltadftcodes/supplmat/SupplMat-FPLO_TplusF.pdf))
in order to achieve plane-wave / (l)apw accuracy. For metals, please choose the
methfessel-paxton smearing for the BZ integration.
