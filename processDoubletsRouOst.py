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
    for roLine in rouOstLineListSorted:
        if prev is not None:
            if roLine.description == prev.description and abs(roLine.wavelength - prev.wavelength) < tolerance:
                # combine lines
                transitions = []
                for ll in [prev, roLine]:
                    ll.status |= ReferenceLineStatus.MERGED
                    transitions.append(ll.transition)
                wavelength = (prev.wavelength + roLine.wavelength)/2.0
                intensity = (prev.intensity + roLine.intensity)/2.0
                transition = f'{transitions[0]}|{transitions[1]}'
                source = prev.source | roLine.source
                status = ReferenceLineStatus.COMBINED
                combinedLine = ReferenceLine(roLine.description, wavelength,
                                             intensity,
                                             status,
                                             transition,
                                             source)
                if wavelength in rouOstDict:
                    raise ValueError(f'Wavelength {wavelength} is already added.')

                # Add new combined line AND original lines (which have been marked as MERGED)
                # These will be removed if necessary later.
                rouOstDict[prev.wavelength] = prev
                rouOstDict[roLine.wavelength] = roLine
                rouOstDict[wavelength] = combinedLine

                prev = None
                continue
            else:
                rouOstDict[prev.wavelength] = prev
        prev = roLine

    linelist = sorted(rouOstDict.values(),
                      key=lambda refLine: refLine.wavelength)

    df = referenceLineSetToDataFrame(linelist)
    rls = ReferenceLineSet(df)
    print(f'Writing output to file {outfile}.')
    rls.writeLineList(outfile)


def main():
    rouOstFileName = os.path.join('derived-data', 'rousselot-osterbrock-merged-linelist.txt')
    # rouOstFileName = 'rousselot-osterbrock-test.txt'
    outfile = 'rousselot-osterbrock-combined.txt'

    print(f'Reading file {rouOstFileName} for orig Rousselot Osterbrock info...')
    combineDoublesFromRouOst(rouOstFileName, outfile)


if __name__ == "__main__":
    main()
