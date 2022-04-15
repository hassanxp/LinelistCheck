import os
from utils import toVacuum, referenceLineSetToDataFrame
from pfs.drp.stella.referenceLine import ReferenceLineSet, ReferenceLine
from pfs.drp.stella.referenceLine import ReferenceLineStatus, ReferenceLineSource


def merge(rouFilename: str, ostFilename: str, outfile: str):

    # If lines from the two catalogues are within this limit,
    # they will be merged.
    # separation = 0.0001
    separation = 0.05

    rouRLS = ReferenceLineSet.fromLineList(rouFilename)
    ostRLS = ReferenceLineSet.fromLineList(ostFilename)

    mergedDict = {}

    rouRLS.sort()
    ostRLS.sort()

    rouIdx = 0
    prevRouIdx = 0
    ostIdx = 0
    prevOstIdx = 0
    mergedLines = 0
    while rouIdx < len(rouRLS)-1 and ostIdx < len(ostRLS)-1:

        # Get PAIRS of lines from each catalog.
        # From the Roussolet catalog, get two consecutive lines of the same intensity.
        # From the Osterbrock catalog, get two consecutive lines of the same transition information (e or f).
        rouRL = rouRLS[rouIdx]
        ostRL = ostRLS[ostIdx]
        prevRouRL = rouRLS[prevRouIdx]
        prevOstRL = ostRLS[prevOstIdx]

        if (rouRL.intensity == prevRouRL.intensity and
            ostRL.transition != 'UNKNOWN' and
            prevOstRL.transition != 'UNKNOWN' and
            ostRL.transition[:6] == prevOstRL.transition[:6] and
            ostRL.transition[7:] == prevOstRL.transition[7:]):

            print(f'Looking at ROU lines w={prevRouRL.wavelength}, {rouRL.wavelength} and intensities {prevRouRL.intensity}, {rouRL.intensity}')
            print(f'           OST lines w={prevOstRL.wavelength}, {ostRL.wavelength} and transitions {prevOstRL.transition}, {ostRL.transition}')
            sep = abs(ostRL.wavelength - rouRL.wavelength)
            print(f'sep={sep}')
            if sep < separation:
                print(f'       merging lines')
                prevRouRL.transition = prevOstRL.transition
                rouRL.transition = ostRL.transition
                prevRouRL.source = ReferenceLineSource.ROUSSELOT2000 | ReferenceLineSource.OSTERBROCK97
                rouRL.source = ReferenceLineSource.ROUSSELOT2000 | ReferenceLineSource.OSTERBROCK97
                if prevOstRL.wavelength in mergedDict.keys():
                    mergedDict.pop(prevOstRL.wavelength)
                mergedDict[prevRouRL.wavelength] = prevRouRL
                mergedDict[rouRL.wavelength] = rouRL
                prevOstIdx = ostIdx
                ostIdx += 1
                prevRouIdx = rouIdx
                rouIdx += 1
                mergedLines += 1
                continue

        if ostRL.wavelength < rouRL.wavelength:
            mergedDict[ostRL.wavelength] = ostRL
            prevOstIdx = ostIdx
            ostIdx += 1
            continue

        mergedDict[rouRL.wavelength] = rouRL
        prevRouIdx = rouIdx
        rouIdx += 1

    # Complete list with remaining lines
    # from either outstanding list
    if rouIdx < len(rouRLS):
        while rouIdx < len(rouRLS):
            rouRL = rouRLS[rouIdx]
            mergedDict[rouRL.wavelength] = rouRL
            rouIdx += 1
    else:
        if ostIdx < len(ostRLS):
            while ostIdx < len(ostRLS):
                ostRL = ostRLS[ostIdx]
                mergedDict[ostRL.wavelength] = ostRL
                ostIdx += 1

    mergedLineSet = ReferenceLineSet.fromRows(mergedDict.values())
    mergedLineSet.writeLineList(outfile)
    print(f'Number of merged lines: {mergedLines}.')
    print(f'Output written to {outfile}.')


def main():
    rouFilename = os.path.join('derived-data', 'rousselot-linelist.txt')
    ostFilename = os.path.join('derived-data', 'osterbrock-linelist.txt')
    outfile = 'rousselot-osterbrock-merged-linelist.txt'
    merge(rouFilename, ostFilename, outfile)


if __name__ == "__main__":
    main()
