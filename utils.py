import re
import numpy as np


def toVacuum(lam_air):
    #
    # n-1 = 0.057921/(238.0185 - 1/lam^2) + 0.001679/(57.362 - 1/lam^2) * (P/101325)*(288.15)/T
    # But setting the pressure and temp terms to unity.
    #
    one_over_lambda_sq = 1/(lam_air*lam_air)
    n_minus_one = (0.057921/(238.0185 - one_over_lambda_sq)
                   + 0.001679/(57.362 - one_over_lambda_sq))

    return lam_air * (1. + n_minus_one)


def readLineList(filename, erin=True):

    sep = ',' if erin else '\s+'

    wavelengthArr = []
    intensityArr = []
    with open(filename) as fd:
        for ii, line in enumerate(fd):
            if erin and ii == 0:
                print(f"An Erin linelist file. Skipping line {ii}")
                continue
            line = re.sub(r"\s*#.*$", "", line).rstrip()  # strip comments
            if not line:
                continue
            fields = re.split(sep, line)
            try:
                if erin:
                    wavelength = fields[1]
                    intensity = fields[2]
                else:
                    wavelength = fields[0]
                    intensity = fields[1]

            except Exception as ex:
                raise RuntimeError(f"Unable to parse line {ii} of {filename}: {ex}")

            try:
                intensity = float(intensity)
            except ValueError:
                intensity = np.nan

            wavelengthArr.append(wavelength)
            intensityArr.append(intensity)
    return wavelengthArr, intensityArr

def createLineList():

    
