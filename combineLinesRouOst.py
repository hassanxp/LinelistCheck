import os
from statistics import mean
from pfs.drp.stella.referenceLine import ReferenceLineSet, ReferenceLine
from pfs.drp.stella.referenceLine import ReferenceLineStatus
from pfs.drp.stella.referenceLine import ReferenceLineSource

from utils import referenceLineSetToDataFrame


def combineDoublesFromRouOst(rouOstFile, outfile, separation=0.35):
    """Combine close doublets from rousselot-osterbrock linelist
    """
    rouOstDict = {}
    rouOstLineList = ReferenceLineSet.fromLineList(rouOstFile)
    rouOstLineListSorted = sorted(rouOstLineList,
                                  key=lambda refLine: refLine.wavelength)

    groupedLines = []
    # Marker for first line in wavelength for group
    start = None
    for roLine in rouOstLineListSorted:
        if start is None:
            start = roLine
            continue
        if (roLine.description == start.description and
           abs(roLine.wavelength - start.wavelength) < separation):
            groupedLines.append(roLine)
            continue
        if len(groupedLines) == 0:
            # Just write start line to output as uncombined line
            rouOstDict[start.wavelength] = start
        else:
            groupedLines.append(start)
            combineLines(groupedLines, rouOstDict)
            groupedLines = []
        start = None

    linelist = sorted(rouOstDict.values(),
                      key=lambda refLine: refLine.wavelength)

    df = referenceLineSetToDataFrame(linelist)
    rls = ReferenceLineSet(df)
    print(f'Writing output to file {outfile}.')
    rls.writeLineList(outfile)


def combineLines(groupedLines, rouOstDict, transition_str_limit=128):
    """Combine grouped lines into a single line.
    The combined line is flagged as COMBINED.
    The input lines are retained, with flag MERGED.

    Note applying limit to transition strings as
    some strings are too long to fit
    into an afw.table when processing (limit 128 char).
    """
    if len(groupedLines) == 0:
        raise ValueError('Line group is unexpectedly empty.')

    # Limit to the allowed length in butler registry of transition string
    transStr = ""
    truncated = False
    wavelengths = []
    intensities = []
    source = ReferenceLineSource.NONE
    for ll in groupedLines:
        ll.status |= ReferenceLineStatus.MERGED
        # Add original lines (which have been marked as MERGED)
        # These will be removed if necessary later.
        rouOstDict[ll.wavelength] = ll
        wavelengths.append(ll.wavelength)
        intensities.append(ll.intensity)
        source |= ll.source
        if truncated is True:
            continue
        if (len(transStr) + len(ll.transition)) > (transition_str_limit-4):
            transStr += '|...'
            truncated = True
            continue
        if len(transStr) != 0:
            transStr += '|'
        transStr += ll.transition

    wavelength = mean(wavelengths)
    intensity = mean(intensities)

    status = ReferenceLineStatus.COMBINED
    combinedLine = ReferenceLine(groupedLines[0].description,
                                 wavelength,
                                 intensity,
                                 status,
                                 transStr,
                                 source)
    if wavelength in rouOstDict:
        # Add a tiny, otherwise insignificant, value
        # to avoid clashes with other lines
        wavelength += 1e-13
        # Check again
        if wavelength in rouOstDict:
            raise ValueError(f'Line with (adjusted) wavelength {wavelength} '
                             'already in dict')
        combinedLine.wavelength = wavelength
    rouOstDict[wavelength] = combinedLine


def main():
    rouOstFileName = os.path.join('derived-data',
                                  'rousselot-osterbrock-merged-linelist.txt')
    # rouOstFileName = 'rousselot-osterbrock-test.txt'
    outfile = 'rousselot-osterbrock-combined.txt'

    print(f'Reading file {rouOstFileName} '
          'for orig Rousselot Osterbrock info...')
    combineDoublesFromRouOst(rouOstFileName, outfile)


if __name__ == "__main__":
    main()
