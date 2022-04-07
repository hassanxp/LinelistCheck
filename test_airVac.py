import lsst.utils.tests
import utils


class AirVacTestCase(lsst.utils.tests.TestCase):
    """Tests the reading of lamp metadata"""

    def setUp(self):
        pass

    def check(self, wavelength_air_ang, wavelength_vac_ang):
        actual_wavelength_vac_ang = utils.toVacuum(wavelength_air_ang/10000.)*10000.
        self.assertAlmostEqual(actual_wavelength_vac_ang, wavelength_vac_ang,
                               delta=1e-3)

    def testConversion(self):
        """Tests air-vac conversion is accurate.
            Using sample data from table 3 of Rousselot+2000
            https://articles.adsabs.harvard.edu/pdf/2000A%26A...354.1134R
        """

        # Note: Rousselot provides values in angstroms.
        # These are converted to microns prior to the actual checks.
        self.check(9972.349, 9975.083)
        self.check(10728.775, 10731.714)
        self.check(12127.748, 12131.067)


if __name__ == "__main__":

    a = AirVacTestCase()
    a.testConversion()
