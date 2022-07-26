import os

from pfs.drp.stella.referenceLine import ReferenceLineSet, ReferenceLine
from pfs.drp.stella.referenceLine import ReferenceLineStatus, ReferenceLineSource

from utils import referenceLineSetToDataFrame


def combineDoublesFromRouOst(rouOstFile, outfile, tolerance=0.5):
    """Combine close doublets from rousselot-osterbrock linelist
    """
    rouOstDict = {}
    rouOstLineList = ReferenceLineSet.fromLineList(rouOstFile)
    rouOstLineListSorted = sorted(rouOstLineList,
                                  key=lambda refLine: refLine.wavelength)

    prev = None
    for line in rouOstLineListSorted:
        if prev is not None:
            if line.description == prev.description and abs(line.wavelength - prev.wavelength) < tolerance:
                # combine lines
                transitions = []
                for ll in [prev, line]:
                    ll.status |= ReferenceLineStatus.MERGED
                    transitions.append(ll.transition)
                wavelength = (prev.wavelength + line.wavelength)/2.0
                intensity = (prev.intensity + line.intensity)/2.0
                combinedLine = ReferenceLine(line.description, wavelength,
                                             intensity,
                                             ReferenceLineStatus.GOOD)
                combinedLine.transition = f'{transitions[0]}|{transitions[1]}'
                combinedLine.source = ReferenceLineSource.ROUSSELOT2000 | ReferenceLineSource.OSTERBROCK97
                combinedLine.status = ReferenceLineStatus.COMBINED
                if wavelength in rouOstDict:
                    raise ValueError(f'Wavelength {wavelength} is already added.')
                rouOstDict[wavelength] = combinedLine
            else:
                rouOstDict[prev.wavelength] = prev
        prev = line

    linelist = sorted(rouOstDict.values(),
                      key=lambda refLine: refLine.wavelength)

    df = referenceLineSetToDataFrame(linelist)
    rls = ReferenceLineSet(df)
    print(f'Writing output to file {outfile}.')
    rls.writeLineList(outfile)


def main():
    rouOstFileName = os.path.join('derived-data', 'rousselot-osterbrock-merged-linelist.txt')
    outfile = 'rousselot-osterbrock-combined.txt'

    print(f'Reading input file {inputfilename}...')
    print(f'Reading file {rouOstFileName} for orig Rousselot Osterbrock info...')

    outfile = 'skyLine_mergedDoublets.txt'
    combineDoublesFromRouOst(rouOstFileName, outfile)


if __name__ == "__main__":
    main()
