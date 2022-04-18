from distutils.command.build_scripts import first_line_re
import os
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


def process(filename: str, rouOstFile: str, outfile: str):

    rouOstLineList = ReferenceLineSet.fromLineList(rouOstFile)
    rouOstDict = {}
    for rl in rouOstLineList:
        rouOstDict[rl.wavelength] = rl

    with open(filename, 'r') as f:
        linelist = []
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
                    linelist.append(refLine)
            commentedLines = []

        # Added modified Rousselot-Osterbrock lines to linelist
        for key, value in rouOstDict.items():
            linelist.append(value)

        # Combine additional lines that cause fitting issues
        wavelengthData = [(771.71340, 771.90879, '9-4_R1f(1.5)|9-4_R1e(1.5)|9-4_R2e(2.5)|9-4_R2f(2.5)|9-4_R2e(3.5)|9-4_R2f(3.5)|9-4_R1e(4.5)|9-4_R1f(4.5)'),
                          (792.3161, 792.33981, '5-1_Q1e(2.5)|5-1_Q1f(2.5)|9-4_P2e(5.5)|9-4_P2f(5.5)'),
                          (876.36880, 876.41160, '7-3_R1f(3.5)|7-3_R1e(3.5)|7-3_R2e(3.5)|7-3_R2f(3.5)'),
                          (967.09880, 967.1322, '8-4_P2e(6.5)|8-4_P2f(6.5)')]

        linelistSorted = sorted(linelist,
                                key=lambda refLine: refLine.wavelength)

        for line in linelistSorted:
            for wmin, wmax, _, in wavelengthData:
                if line.wavelength >= wmin and line.wavelength <= wmax:
                    line.status |= ReferenceLineStatus.MERGED

        # Add replacement COMBINED lines for the above
        for wmin, wmax, transition in wavelengthData:
            linelistSorted.append(ReferenceLine(description='OH',
                                  wavelength=0.5 * (wmin + wmax),
                                  intensity=1.0,
                                  status=ReferenceLineStatus.COMBINED,
                                  transition=transition,
                                  source=ReferenceLineSource.ROUSSELOT2000|ReferenceLineSource.OSTERBROCK97))

        # Re-sort lines
        lineListSorted = sorted(linelistSorted,
                                key=lambda refLine: refLine.wavelength)

        # Mark lines that are not visible (faint, out of PFS range) to be used
        print('Flagging non-visible lines in full list...')
        nFaintLines = 0
        for line in lineListSorted:
            if line.intensity < 0.1 or line.wavelength < 630.00 or line.source == ReferenceLineSource.OSTERBROCK97:
                line.status = ReferenceLineStatus.NOT_VISIBLE
                nFaintLines += 1

        print(f'Flagged {nFaintLines} not visible lines.')

        df = referenceLineSetToDataFrame(lineListSorted)
        rls = ReferenceLineSet(df)
        print(f'Writing output to file {outfile}.')
        rls.writeLineList(outfile)


def main():
    inputfilename = os.path.join('derived-data', 'skyLines-1fc7b67.txt')
    rouOstFileName = os.path.join('derived-data', 'rousselot-osterbrock-merged-linelist.txt')

    print(f'Reading input file {inputfilename}...')
    print(f'Reading file {rouOstFileName} for orig Rousselot Osterbrock info...')

    outfile = 'skyLine_mergedDoublets.txt'
    process(inputfilename, rouOstFileName, outfile)


if __name__ == "__main__":
    main()
