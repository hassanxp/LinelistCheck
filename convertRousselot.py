import os
from utils import referenceLineSetToDataFrame
from pfs.drp.stella.referenceLine import ReferenceLineSet, ReferenceLine
from pfs.drp.stella.referenceLine import ReferenceLineStatus, ReferenceLineSource


def convert(filename: str, outfile: str, infofile: str):
    with open(filename, 'r') as f:
        linelist = []
        for line in f:
            lineStrp = line.strip()
            if lineStrp.startswith('#'):
                continue
            fields = lineStrp.split()
            waveVacAngs, intensity = fields
            trans = "UNKNOWN"
            waveVacNm = float(waveVacAngs)/10.
            refLine = ReferenceLine(description='OH',
                                    wavelength=waveVacNm,
                                    intensity=float(intensity),
                                    status=ReferenceLineStatus.GOOD,
                                    transition=trans,
                                    source=ReferenceLineSource.ROUSSELOT2000)
            linelist.append(refLine)

        lineListSorted = sorted(linelist,
                                key=lambda refLine: refLine.wavelength)

        df = referenceLineSetToDataFrame(lineListSorted)
        rls = ReferenceLineSet(df)
        rls.writeLineList(outfile)


def main():
    filename = 'Rousellot_list_v2.0.dat'
    fullfilename = os.path.join('data', 'Rousselot+2000', filename)
    print(f'Reading input file {fullfilename}...')
    outfile = 'rousselot-linelist.txt'
    infofile = 'rousselot-info.txt'
    convert(fullfilename, outfile, infofile)


if __name__ == "__main__":
    main()
