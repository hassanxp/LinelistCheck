import re
import numpy as np
import pandas as pd
from pfs.drp.stella.referenceLine import ReferenceLineSet, ReferenceLine


def referenceLineSetToDataFrame(refLineSet):
    d = []
    for rl in refLineSet:
        d.append(
            {
                'description': rl.description,
                'wavelength': rl.wavelength,
                'intensity': rl.intensity,
                'status': rl.status,
                # 'transition': rl.transition,
                # 'source': rl.source,
            })
    return pd.DataFrame(d)


def toVacuum(lam_air_um):
    """
    Converts a wavelength measured in air to vacuum.

    Following is taken from
    https://pfspipe.ipmu.jp/jira/browse/PIPE2D-935
    Originally from Jim Gunn, but quoted by N Caplar:

    The issue of the vacuum vs air wavelengths is complicated, not
    by the physics particularly, but by convention. A simple Sellmeier
    equation for the index is
    n-1 = 0.057921/(238.0185 - 1/lam^2) + 0.001679/(57.362 - 1/lam^2) * (P/101325)*(288.15)/T
    with the wavelength lam in microns, pressure in pascals, temeperature in
    Kelvin; if you set the pressure and temperature terms to unity,
    you get the index for 15C and the standard sea level pressure, at which I think
    the difference in vacuum and `air' wavelength is DEFINED, and this
    is what you want to use, I think, not the ACTUAL air wavelength, on MK,
    which is at 0.6 bar pressure and 0C. A very useful reference is the SDSS
    APOGEE one, which uses a different Sellmeier form with 5 constants
    instead of 4 (there is a constant term of about 7e-5),
    and with a nice discussion, is to be found at
    https://www.as.utexas.edu/~hebe/apogee/docs/air_vacuum.pdf,
    also attached.

    Using

    n-1 = 0.057921/(238.0185 - 1/lam^2)
    + 0.001679/(57.362 - 1/lam^2) * (P/101325)*(288.15)/T
    But setting the pressure and temp terms to unity.

    Parameters:
        lam_air_um : `float`
            The input wavelength, measured in air [micron]

    Return:
        lam_vac_um : `float`
            The converted wavelength in vacuum [micron]
    """

    one_over_lambda_sq = 1/(lam_air_um*lam_air_um)
    n_minus_one = (0.057921/(238.0185 - one_over_lambda_sq)
                   + 0.001679/(57.362 - one_over_lambda_sq))

    return lam_air_um * (1. + n_minus_one)


def readLineList(filename, erin=True):

    sep = r',' if erin else r'\s+'

    wavelengthArr = []
    intensityArr = []
    speciesArr = []
    statusArr = []
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
                    species = 'UNKNOWN'
                    status = 0
                else:
                    wavelength = float(fields[0])
                    intensity = float(fields[1])
                    species = fields[2]
                    status = int(fields[3])

            except Exception as ex:
                raise RuntimeError(f"Unable to parse line {ii} of {filename}: {ex}")

            try:
                intensity = float(intensity)
            except ValueError:
                intensity = np.nan

            wavelengthArr.append(wavelength)
            intensityArr.append(intensity)
            speciesArr.append(species)
            statusArr.append(status)
    return wavelengthArr, intensityArr, speciesArr, statusArr


def generateAirLineList(filename, outFile,
                        infoFileName, lambdaNm=True,
                        erin=True):
    wavelength_air, intensity, species, status = readLineList(filename,
                                                              erin=erin)
    linelist = []
    with open(infoFileName, 'w') as infoFile:
        if lambdaNm is False:
            infoFile.write("lambda_air[angstrom], lambda_vac[nm]\n")
        else:
            infoFile.write("lambda_air[nm], lambda_vac[nm]\n")

        infoFile.write("\n")
        for w_air, ii, sp, st in zip(wavelength_air, intensity, species, status):

            if lambdaNm is False:
                # convert wavelengths from angstrom to microns
                w_air_micron = float(w_air)/10000.0
            else:
                # convert wavelengths from nm to microns
                w_air_micron = float(w_air)/1000.00

            w_vac_nm = toVacuum(w_air_micron) * 1000.0
            infoFile.write(f'{w_air}, {w_vac_nm}\n')
            refLine = ReferenceLine(sp, w_vac_nm, ii, st)
            linelist.append(refLine)

    lineListSorted = sorted(linelist, key=lambda refLine: refLine.wavelength)

    df = referenceLineSetToDataFrame(lineListSorted)
    rls = ReferenceLineSet(df)
    rls.writeLineList(outFile)
