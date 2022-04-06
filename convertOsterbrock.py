import os
from utils import toVacuum, referenceLineSetToDataFrame
from pfs.drp.stella.referenceLine import ReferenceLineSet, ReferenceLine
from pfs.drp.stella.referenceLine import ReferenceLineStatus, ReferenceLineSource


def convert(filename: str, outfile: str, infofile: str):
    with open(filename, 'r') as f, open(infofile, 'w') as iff:
        iff.write("lambda_air[angstrom], lambda_vac[nm], transition\n")
        iff.write("\n")

        linelist = []
        for line in f:
            fields = line.strip().split()
            waveAirAngs, trans1, trans2 = fields
            trans = f'{trans1}_{trans2}'
            waveVacNm = toVacuum(float(waveAirAngs)/10000.)*1000.
            iff.write(f'{waveAirAngs}, {waveVacNm}, {trans}\n')
            refLine = ReferenceLine(description='OH',
                                    wavelength=waveVacNm,
                                    intensity=1.0,
                                    status=ReferenceLineStatus.GOOD,
                                    transition=trans,
                                    source=ReferenceLineSource.OSTERBROCK97)
            linelist.append(refLine)

        lineListSorted = sorted(linelist,
                                key=lambda refLine: refLine.wavelength)

        df = referenceLineSetToDataFrame(lineListSorted)
        rls = ReferenceLineSet(df)
        rls.writeLineList(outfile)


def main():
    filename = 'oh.dat'
    fullfilename = os.path.join('data', 'Osterbrock+1997', filename)
    print(f'Reading input file {fullfilename}...')
    outfile = 'osterbrock-linelist.csv'
    infofile = 'osterbrock-info.csv'
    convert(fullfilename, outfile, infofile)


if __name__ == "__main__":
    main()
