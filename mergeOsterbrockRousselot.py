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

    mergedList = []

    rouRLS.sort()
    ostRLS.sort()

    rouIdx = 0
    ostIdx = 0
    mergedLines = 0
    while rouIdx < len(rouRLS) and ostIdx < len(ostRLS):
        rouRL = rouRLS[rouIdx]
        ostRL = ostRLS[ostIdx]
        print(f'Looking at Ost wavelength {ostRL.wavelength} and Rou wavelength {rouRL.wavelength}..')

        if abs(ostRL.wavelength - rouRL.wavelength) < separation:
            print(f'    merging line w. Ost wavelength {ostRL.wavelength} and Rou wavelength {rouRL.wavelength})')
            rouRL.transition = ostRL.transition
            rouRL.source = ReferenceLineSource.ROUSSELOT2000 | ReferenceLineSource.OSTERBROCK97
            mergedList.append(rouRL)
            ostIdx += 1
            rouIdx += 1
            mergedLines += 1
            continue

        if ostRL.wavelength < rouRL.wavelength:
            mergedList.append(ostRL)
            ostIdx += 1
            continue

        mergedList.append(rouRL)
        rouIdx += 1

    # Complete list with remaining lines
    # from either outstanding list
    if rouIdx < len(rouRLS):
        while rouIdx < len(rouRLS):
            rouRL = rouRLS[rouIdx]
            mergedList.append(rouRL)
            rouIdx += 1
    else:
        if ostIdx < len(ostRLS):
            while ostIdx < len(ostRLS):
                ostRL = ostRLS[ostIdx]
                mergedList.append(ostRL)
                ostIdx += 1

    mergedLineSet = ReferenceLineSet.fromRows(mergedList)
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
