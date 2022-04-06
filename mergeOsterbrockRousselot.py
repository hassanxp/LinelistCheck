import os
from utils import toVacuum, referenceLineSetToDataFrame
from pfs.drp.stella.referenceLine import ReferenceLineSet, ReferenceLine
from pfs.drp.stella.referenceLine import ReferenceLineStatus, ReferenceLineSource


def merge(rouFilename: str, ostFilename: str, outfile: str):
    rouRLS = ReferenceLineSet.fromLineList(rouFilename)
    ostRLS = ReferenceLineSet.fromLineList(ostFilename)

    mergedList = []

    rouRLS.sort()
    ostRLS.sort()

    rouIdx = 0
    ostIdx = 0
    while rouIdx < len(rouRLS) and ostIdx < len(ostRLS):
        rouRL = rouRLS[rouIdx]
        ostRL = ostRLS[ostIdx]
        print(f'Looking at Ost wavelength {ostRL.wavelength} and Rou wavelength {rouRL.wavelength}..')

        if abs(ostRL.wavelength - rouRL.wavelength) < 0.0001:
            print(f'    merging line w. Ost wavelength {ostRL.wavelength} and Rou wavelength {rouRL.wavelength})')
            rouRL.transition = ostRL.transition
            rouRL.source = ReferenceLineSource.ROUSSELOT2000EXT
            mergedList.append(rouRL)
            ostIdx += 1
            rouIdx += 1
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


def main():
    rouFilename = 'rousselot-linelist.csv'
    ostFilename = 'osterbrock-linelist.csv'
    outfile = 'rousselot-osterbrock-merged-linelist.csv'
    merge(rouFilename, ostFilename, outfile)


if __name__ == "__main__":
    main()
