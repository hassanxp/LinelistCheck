from distutils.command.build_scripts import first_line_re
import os
import wave
from utils import toVacuum, referenceLineSetToDataFrame
from pfs.drp.stella.referenceLine import ReferenceLineSet, ReferenceLine
from pfs.drp.stella.referenceLine import ReferenceLineStatus, ReferenceLineSource


def parse(line: str) -> ReferenceLine:
    fields = line.split()
    wavelength, intensity, species, status = fields[:4]
    return ReferenceLine(description=species,
                         wavelength=float(wavelength),
                         intensity=float(intensity),
                         status=int(status),
                         transition=None,
                         source=None)


def process(prevSkyLineName: str, rouOstFile: str, nonOHFile: str, outfile: str):

    rouOstLineList = ReferenceLineSet.fromLineList(rouOstFile)
    rouOstDict = {}
    for rl in rouOstLineList:
        rouOstDict[rl.wavelength] = rl

    # lineDict = processPrevSkyLine(prevSkyLineName, rouOstDict)
    lineDict = {}

    # Add MERGED Rousselot-Osterbrock lines to linelist
    for refLine in rouOstDict.values():
        if refLine.wavelength in lineDict:
            print(f'WARNING: Wavelength {refLine.wavelength} is already added. Skipping MERGED line')
        lineDict[refLine.wavelength] = refLine

    linelist = sorted(lineDict.values(),
                      key=lambda refLine: refLine.wavelength)

    # linelist = combineExplicitGroups(linelist)

    # Add non-OH lines.
    nonOHLineList = ReferenceLineSet.fromLineList(nonOHFile)
    for rl in nonOHLineList:
        linelist.append(rl)

    # Re-sort lines
    linelist = sorted(linelist,
                      key=lambda refLine: refLine.wavelength)

    # Mark lines that are not visible (faint, out of PFS range) to be used
    print('Flagging non-visible lines in full list...')
    nFaintLines = 0
    for line in linelist:
        # if line.intensity < 0.1 or line.wavelength < 630.00 or line.source == ReferenceLineSource.OSTERBROCK97:
        if line.intensity < 0.1:
            line.status |= ReferenceLineStatus.NOT_VISIBLE
            nFaintLines += 1

    print(f'Flagged {nFaintLines} not visible lines.')

    # Check for duplicates
    checkDuplicates(linelist)

    # Remove all lines marked as MERGED temporarily until DRP updated to support these
    filteredLineList = []
    for line in linelist:
        if line.status & ReferenceLineStatus.MERGED == 0:
            filteredLineList.append(line)

    df = referenceLineSetToDataFrame(filteredLineList)
    rls = ReferenceLineSet(df)
    print(f'Writing output to file {outfile}.')
    rls.writeLineList(outfile)


def checkDuplicates(linelist):
    print('Checking for any duplicates in wavelength..')
    wavelengthDict = {}
    for line in linelist:
        wavelength = line.wavelength
        if wavelength in wavelengthDict:
            print(f'Wavelength {wavelength} is already added.')
        wavelengthDict[wavelength] = line


def processPrevSkyLine(filename, rouOstDict):
    with open(filename, 'r') as f:
        lineDict = {}
        commentedLines = []

        print('Combining doublets..')
        for line in f:
            # Look for sequence consisting of:
            # 1. Two commented lines
            # 2. A uncommented single entry that ends with [5]*, eg
            #
            #  # 613.82390      0.126  OH         0 # [5]
            #  # 613.82650      0.126  OH         0 # [5]
            #    613.82520      0.126  OH         0 # [5]*
            lineStripped = line.strip()
            if lineStripped.startswith('#'):
                commentedLines.append(lineStripped[1:])
                continue

            if len(commentedLines) == 2:
                if lineStripped.endswith('[5]*'):
                    transitions = []
                    for ii in [0, 1]:
                        commentedRefLine = parse(commentedLines[ii])
                        rouOstLine = rouOstDict[commentedRefLine.wavelength]
                        rouOstLine.status = ReferenceLineStatus.MERGED
                        transition = rouOstLine.transition
                        transitions.append(transition)
                    refLine = parse(lineStripped)
                    refLine.transition = f'{transitions[0]}|{transitions[1]}'
                    refLine.source = ReferenceLineSource.ROUSSELOT2000 | ReferenceLineSource.OSTERBROCK97
                    refLine.status = ReferenceLineStatus.COMBINED
                    wavelength = refLine.wavelength
                    if wavelength in lineDict:
                        raise ValueError(f'Wavelength {wavelength} is already added.')
                    lineDict[refLine.wavelength] = refLine
            commentedLines = []
    return lineDict


def combineExplicitGroups(linelist):
    """Combine additional lines that cause fitting issues
    Note some strings are too long to fit into an afw.table when processing (limit 128 char).
    Using ellipses to truncate.
    """
    wavelengthData = [
                    # (682.93300, 683.18290, '7-2_R1f(3.5)|7-2_R1e(3.5)|7-2_R1f(4.5)|7-2_R1e(4.5)|7-2_R1f(2.5)|7-2_R1e(2.5)|7-2_R2e(3.5)|7-2_R2f(3.5)|7-2_R2e(2.5)|7-2_R2f(2.5)'),
                    (682.93300, 683.18290, '7-2_R1f(3.5)|7-2_R1e(3.5)|...|7-2_R2f(3.5)|7-2_R2e(2.5)|7-2_R2f(2.5)'),
                    (724.21610, 724.39640, '8-3_R1f(2.5)|8-3_R1e(2.5)|8-3_R1e(4.5)|8-3_R2e(2.5)|8-3_R2e(3.5)|8-3_R2f(3.5)'),
                    (724.68910, 724.86120, '8-3_R1f(1.5)|8-3_R1e(1.5)|8-3_R2e(1.5)|8-3_R2f(1.5)|8-3_R1e(5.5)|8-3_R1f(5.5)|8-3_R2e(4.5)|8-3_R2f(4.5)'),
                    (752.47230, 752.69161, '4-0_Q2f(0.5)|4-0_Q2e(0.5)|4-0_Q1f(1.5)|4-0_Q1e(1.5)|UNKNOWN'),
                    (771.43290, 771.47490, '4-0_P1e(6.5)|4-0_P1f(6.5)|9-4_R1f(3.5)|9-4_R1e(3.5)|9-4_R1f(2.5)|9-4_R1e(2.5)'),
                    (771.71340, 771.90879, '9-4_R1f(1.5)|9-4_R1e(1.5)|9-4_R2e(2.5)|9-4_R2f(2.5)|9-4_R2e(3.5)|9-4_R2f(3.5)|9-4_R1e(4.5)|9-4_R1f(4.5)'),
                    (772.56500, 772.58400, '9-4_R2e(1.5)|9-4_R2f(1.5)|9-4_R2e(4.5)|9-4_R2f(4.5)|9-4_R1e(5.5)|9-4_R1f(5.5)'),
                    (785.12800, 785.27320, '9-4_P2e(3.5)|9-4_P2f(3.5)|5-1_R1f(5.5)|5-1_R1e(5.5)|'),
                    (785.41180, 785.58040, '5-1_R2e(3.5)|5-1_R2f(3.5)|9-4_P1e(4.5)|9-4_P1f(4.5)|5-1_R1f(3.5)|5-1_R1e(3.5)'),
                    (792.31610, 792.33981, '5-1_Q1e(2.5)|5-1_Q1f(2.5)|9-4_P2e(5.5)|9-4_P2f(5.5)|5-1_R1f(4.5)|5-1_R1e(4.5)|5-1_R2e(4.5)|5-1_R2f(4.5)'),
                    (802.78750, 803.12090, '5-1_P1e(4.5)|5-1_P1f(4.5)|9-4_P1e(8.5)|9-4_P1f(8.5)'),
                    (840.14290, 840.21970, '6-2_P1e(2.5)|6-2_P1f(2.5)|UNKNOWN'),
                    (828.05420, 828.08830, '6-2_R1f(4.5)|6-2_R1e(4.5)|6-2_R1f(5.5)|6-2_R1e(5.5)'),
                    (828.22160, 828.54940, '6-2_R2e(4.5)|6-2_R2e(3.5)|6-2_R2f(3.5)|6-2_R2f(4.5)|6-2_R1e(6.5)|6-2_R1f(6.5)|6-2_R2e(5.5)|6-2_R2f(5.5)'),
                    # (876.10960, 876.89330, '7-3_R1f(4.5)|7-3_R1e(4.5)|7-3_R1f(5.5)|7-3_R1e(5.5)|7-3_R1f(3.5)|7-3_R1e(3.5)|7-3_R1f(3.5)|7-3_R1e(3.5)|7-3_R2e(3.5)|7-3_R2f(3.5)|7-3_R2e(2.5)|7-3_R2f(2.5)'),
                    (876.10960, 876.61320, '7-3_R1f(4.5)|7-3_R1e(4.5)|7-3_R1f(5.5)|7-3_R1e(5.5)|...|7-3_R2e(3.5)|7-3_R2f(3.5)|6-2_P2f(9.5)|6-2_P2e(9.5)'),
                    (876.84700, 877.27710, '7-3_R1e(6.5)|7-3_R1f(6.5)|7-3_R1f(2.5)|7-3_R1e(2.5)'),
                    (886.97680, 887.02760, '7-3_Q1f(4.5)|7-3_P2e(1.5)|7-3_P2f(1.5)|7-3_Q1e(4.5)'),
                    (930.85350, 931.19550, '8-4_R1f(4.5)|8-4_R1e(4.5)|8-4_R1f(3.5)|8-4_R1e(3.5)|8-4_R2e(3.5)|8-4_R2f(3.5)'),
                    (931.32560, 931.59310, '8-4_R2e(4.5)|8-4_R2f(4.5)|8-4_R2e(2.5)|8-4_R2f(2.5)|8-4_R1f(2.5)|8-4_R1e(2.5)'),
                    (932.02200, 932.05450, '8-4_R1e(6.5)|8-4_R1f(6.5)|8-4_R2e(5.5)|8-4_R2f(5.5)'),
                    (944.21860, 944.43420, '8-4_P1e(2.5)|8-4_P1f(2.5)|UNKNOWN'),
                    (947.93479, 948.20250, '8-4_P1e(3.5)|8-4_P1f(3.5)|UNKNOWN'),
                    (967.09880, 967.13220, '8-4_P2e(6.5)|8-4_P2f(6.5)')]

    for idx, wD in enumerate(wavelengthData):
        transition = wD[2]
        if len(transition) >= 128:
            raise ValueError(f'entry {idx} has a transition of length > 128 ({len(transition)})')
        if ' ' in transition:
            raise ValueError(f'transition [{transition}] contains spaces.')

    for line in linelist:
        wavelength = line.wavelength
        for wmin, wmax, _, in wavelengthData:
            if wavelength >= wmin and wavelength <= wmax:
                print(f'Assigning wavelength {wavelength} as MERGED')
                line.status |= ReferenceLineStatus.MERGED

    # Add replacement COMBINED lines for the above
    for wmin, wmax, transition in wavelengthData:
        linelist.append(ReferenceLine(description='OH',
                                      wavelength=0.5 * (wmin + wmax),
                                      intensity=1.0,
                                      status=ReferenceLineStatus.COMBINED,
                                      transition=transition,
                                      source=ReferenceLineSource.ROUSSELOT2000|ReferenceLineSource.OSTERBROCK97))

    return linelist


def main():
    inputfilename = os.path.join('derived-data', 'skyLines-1fc7b67.txt')
    rouOstFileName = os.path.join('derived-data', 'rousselot-osterbrock-combined.txt')
    nonOHfile = os.path.join('derived-data', 'nonOHlines.txt')

    print(f'Reading input file {inputfilename}...')
    print(f'Reading file {rouOstFileName} for orig Rousselot Osterbrock info...')

    outfile = 'skyLine_mergedDoublets.txt'
    process(inputfilename, rouOstFileName, nonOHfile, outfile)


if __name__ == "__main__":
    main()
