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
                        transition = rouOstDict[commentedRefLine.wavelength]
                        transitions.append(transition)
                    refLine = parse(lineStripped)
                    refLine.transition = f'{transitions[0]}|{transitions[1]}'
                    refLine.source = ReferenceLineSource.ROUSSELOT2000 | ReferenceLineSource.OSTERBROCK97
                    refLine.status = ReferenceLineStatus.MERGED
                    linelist.append(refLine)
            commentedLines = []

        df = referenceLineSetToDataFrame(linelist)
        rls = ReferenceLineSet(df)
        print(f'Writing output to file {outfile}.')
        rls.writeLineList(outfile)


def main():
    inputfilename = os.path.join('derived-data', 'skyLines-1fc7b67.txt')
    rouOstFileName = os.path.join('derived-data', 'rousselot-osterbrock-merged-linelist.txt')

    print(f'Reading input file {inputfilename}...')
    print(f'Reading file {rouOstFileName} for orig Rousselot Osterbrock info...')

    outfile = 'out.txt'
    process(inputfilename, rouOstFileName, outfile)


if __name__ == "__main__":
    main()