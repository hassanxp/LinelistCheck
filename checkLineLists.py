import os
import .utils


def orderMatches(value, data):
    #
    # Ranks values in order of their distance to the input value
    #
    return sorted([float(d) for d in data], key=lambda x: abs(x-value))


def main():
    ref_linelist_dir = r'/Users/hassans/PFS/visualStudioCode/obs_pfs/pfs/lineLists'
    pfilamp_linelist_dir = r'/Users/hassans/Downloads'

    for lamp in ['Ar', 'Kr', 'Ne', 'Xe']:
        print(f'Lamp: {lamp}-------------')
        wavelength_ref, _ = readLineList(os.path.join(ref_linelist_dir,
                                                      f'{lamp}.txt'),
                                         erin=False)
        wavelength_air, _ = readLineList(os.path.join(pfilamp_linelist_dir,
                                                      f'new_{lamp}_linemeas.csv'),
                                         erin=True)
        for w in wavelength_air:
            w_air = float(w)
            w_vac = toVacuum(w_air)
            closest_match = orderMatches(w_vac, wavelength_ref)[0]
            diff = abs(closest_match-w_vac)
            print(f'wavelength_air={w_air}, wavelength_vac={w_vac}, '
                  f'matches={closest_match}, diff={diff}')


if __name__ == "__main__":
    main()
