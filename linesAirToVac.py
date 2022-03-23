from utils import toVacuum

def main():
    # for airWl in [557.734, 630.003, 636.378, 704.12510]:
    for airWl in [
                    704.12510,
                    704.14730,
                    704.78460,
                    704.83760,
                    708.23700,
                    708.26830,
                    708.89920,
                    712.72480,
                    712.76610,
                    713.25530,
                    713.33460,
                    717.60050,
                    717.65280,
                    718.08080,
                    718.17500,
                    ]:
        vacWl = toVacuum(airWl)
        # print(f'air wavelength: {airWl}, vacuum wavelength: {vacWl:0.5f}')
        print(f'{vacWl:0.5f}')


if __name__ == "__main__":
    main()
