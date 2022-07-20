import os
from utils import generateAirLineList


def main():
    pfilamp_linelist_dir = r'./data'
    lamp = 'HgCd'
    filename = os.path.join(pfilamp_linelist_dir, f'new_{lamp}_linemeas.csv')
    generateAirLineList(filename, 'out.csv', 'info.txt')


if __name__ == "__main__":
    main()
