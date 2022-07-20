import os

from pfs.drp.stella.referenceLine import ReferenceLineSet


def main():
    inputfile = 'skyLines-master.txt'
    inputpath = os.path.join('derived-data', inputfile)

    print(f'Reading input file {inputpath}...')

    rls = ReferenceLineSet.fromLineList(inputpath)
    outfile = f'{inputfile}-conv'
    rls.writeLineList(outfile)

    print(f'Wrote to file {outfile}.')


if __name__ == "__main__":
    main()
