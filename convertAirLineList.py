from utils import generateAirLineList


def main():
    filename = 'osterbrock1996_order_44_45.csv'
    print(f'Reading input file {filename}...')
    outfile = f'{filename}-vac'
    infoFile = f'{filename}-info.csv'
    generateAirLineList(filename, outfile, infoFile,
                        lambdaNm=False, erin=False)
    print(f'Output written to {outfile}. Conversion table in {infoFile}')


if __name__ == "__main__":
    main()
