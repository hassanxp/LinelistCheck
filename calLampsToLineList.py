import os
from utils import generateAirLineList


def main():
    pfilamp_linelist_dir = r'./data'

    filename = os.path.join(pfilamp_linelist_dir, f'new_{lamp}_linemeas.csv')
    generateAirLineList(filename)


if __name__ == "__main__":
    main()
