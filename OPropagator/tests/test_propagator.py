from src.orbit.initial_conditions import InitOrbit
from src.orbit import utils
from src.orbit.models import InitPropagator
from org.orekit.propagation.analytical.tle import TLE

# Utils
from org.orekit.utils import Constants


def tle_orbit_setup(tle_line_1, tle_line_2, orbit_type="TLE", timezone="UTC"):
    # Initialize orbit
    initOrbit = InitOrbit(
        R_i=Constants.WGS84_EARTH_EQUATORIAL_RADIUS, mu_i=Constants.WGS84_EARTH_MU
    )
    initOrbit.set_reference_frame()  # Reference frame
    initOrbit.set_timescale(timezone)

    tle_string = TLE(tle_line_1, tle_line_2)
    tle_epoch = tle_string.getDate()

    components = tle_epoch.getComponents(initOrbit.timescale)
    initOrbit.set_epoch(
        year=components.getDate().getYear(),
        month=components.getDate().getMonth(),
        day=components.getDate().getDay(),
        hour=components.getTime().getHour(),
        minute=components.getTime().getMinute(),
        second=int(components.getTime().getSecond()),
        microsecond=(
            components.getTime().getSecond() - int(components.getTime().getSecond())
        )
        * 1e6,
    )
    initOrbit.set_orbit(
        orbit_type, anomaly_type="MEAN", tle_str=tle_string
    )  # Orbital parameters

    return initOrbit, tle_string


def propagator_setup(orbit_obj, pName, mass=6900.0, **kwargs):
    initPropagator = InitPropagator(
        orbit_obj.orbit, mass, step_type="FIXED"
    )  # Model setup
    if pName in [
        "Euler",
        "Midpoint",
        "Gill",
        "ThreeEighthes",
        "Luther",
        "HighamHall54",
        "DormandPrince54",
        "DormandPrince853",
        "GraggBulirschStoer",
        "AdamsBashforth",
        "AdamsMoulton",
    ]:
        # Numerical model - DormandPrince853
        initPropagator.set_propagation(
            "NUMERICAL",
            pName,
            iStep=60.0,
            iType="ADAPTIVE",
            iSolver="NONSTIFF",
        )
        propagator, states = initPropagator.set_model(
            model="N_BODY", bodies=[], forces=["GRAVITY"]
        )
    elif pName in ["Kepler", "Eckstein-Heschler", "Brouwer-Lyddane", "TLE"]:
        # High fidelity model - SGP4
        initPropagator.set_propagation(
            "ANALYTICAL",
            pName,
            tle_str=kwargs.get("tle_str"),
            iStep=60.0,
        )
        propagator, states = initPropagator.set_model(model="TWO_BODY")

    return propagator, states


def run(line_1, line_2, mass):
    """Run the test propagator script"""
    test_orbit, test_orbit_tle = tle_orbit_setup(
        line_1, line_2, orbit_type="TLE"
    )  # Orbit setup

    # NUMERICAL METHOD (INTEGRATOR)
    numerical_propagator, numerical_states = propagator_setup(
        test_orbit,
        "DormandPrince853",
        mass,
    )  # Propagator setup
    numerical_propagator.propagate(
        test_orbit.orbit.getDate(),
        test_orbit.orbit.getDate().shiftedBy(60.0 * 60.0 * 24.0 * 1.0),
    )  # Propagate

    # ANALYTICAL METHOD (TLE)
    analytical_propagator, analytical_states = propagator_setup(
        test_orbit,
        "TLE",
        mass,
        tle_str=test_orbit_tle,
    )  # Propagator setup
    analytical_propagator.propagate(
        test_orbit.orbit.getDate(),
        test_orbit.orbit.getDate().shiftedBy(60.0 * 60.0 * 24.0 * 1.0),
    )  # Propagate

    # Plot
    utils.plot_error(estimated=numerical_states, actual=analytical_states)
