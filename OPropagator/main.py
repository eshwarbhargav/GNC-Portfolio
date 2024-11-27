import orekit

from orekit.pyhelpers import setup_orekit_curdir

from tests import test_propagator


def main():
    """Main function"""
    orekit.initVM()
    setup_orekit_curdir()

    # SOYUZ-MS 26
    tle_line_1 = "1 61043U 24162A   24323.88118015  .00018295  00000+0  32819-3 0  9995"
    tle_line_2 = "2 61043  51.6392 268.9875 0007344 233.7782 126.2529 15.49944020482499"

    # SOYUZ-MS 26 - Mass
    sc_mass = 6900.0
    test_propagator.run(tle_line_1, tle_line_2, mass=sc_mass)


if __name__ == "__main__":
    main()
