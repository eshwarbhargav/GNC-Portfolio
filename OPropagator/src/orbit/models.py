# Standards
from org.orekit.frames import FramesFactory
from org.orekit.utils import Constants, IERSConventions
from org.orekit.time import TimeScalesFactory

# Hipparchus library
from org.hipparchus.ode.nonstiff import (
    HighamHall54Integrator,
    DormandPrince54Integrator,
    DormandPrince853Integrator,
    GraggBulirschStoerIntegrator,
    AdamsBashforthIntegrator,
    AdamsMoultonIntegrator,
)
from org.hipparchus.ode.nonstiff import (
    EulerIntegrator,
    MidpointIntegrator,
    ClassicalRungeKuttaIntegrator,
    GillIntegrator,
    ThreeEighthesIntegrator,
    LutherIntegrator,
)

# File-specific
from orekit import JArray_double
from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
from org.orekit.forces.drag import DragForce, IsotropicDrag
from org.orekit.forces.gravity import (
    ThirdBodyAttraction,
    HolmesFeatherstoneAttractionModel,
    Relativity,
    OceanTides,
    SolidTides,
)
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.models.earth.atmosphere import NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.analytical import (
    KeplerianPropagator,
    EcksteinHechlerPropagator,
    BrouwerLyddanePropagator,
)
from org.orekit.propagation.analytical.tle import TLEPropagator
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.propagation.sampling import (
    PythonOrekitFixedStepHandler,
    PythonOrekitStepHandler,
)

from src.orbit import utils
from src.orbit.initial_conditions import orbit_wrapper


class FixedStepHandler(PythonOrekitFixedStepHandler):
    """
    FixedStepHandler The purpose of PythonOrekitFixedStepHandler is to perform step handling during propagation of NUMERICAL and ANALYTICAL models at desired fixed step.

    Args:
        PythonOrekitFixedStepHandler (_type_): _description_
    """

    def __init__(self, orbit_type):
        self.states = {
            "TIME": [],
            "CARTESIAN": {
                "position": [],
                "velocity": [],
                "acceleration": [],
            },
            "KEPLERIAN": {
                "a": [],
                "e": [],
                "i": [],
                "raan": [],
                "pa": [],
                "true_anomaly": [],
            },
        }
        self.orbit_type = orbit_type
        self.counter = 0
        super(FixedStepHandler, self).__init__()

    def init(self, s0, t, step):
        print(f"{self.orbit_type} STATES PROPAGATION")
        utils.print_state_header(self.orbit_type)

    def handleStep(self, currentState):
        wrappedOrbit, keplerian_orbit = orbit_wrapper(currentState)
        utils.update_states(currentState, keplerian_orbit, self.states)

        if self.counter > -1:
            utils.print_state(wrappedOrbit)
            self.counter += 1

    def finish(self, finalState):
        print("this was the last step")
        # utils.plot_states(self.states, save_fig=True, pType=self.orbit_type)


class VariableStepHandler(PythonOrekitStepHandler):
    """
    Note: The purpose of PythonOrekitStepHandler is to perform step handling during propagation of NUMERICAL models not ANALYTICAL ones.
    """

    def __init__(self, orbit_type):
        self.states = {
            "TIME": [],
            "CARTESIAN": {
                "position": [],
                "velocity": [],
                "acceleration": [],
            },
            "KEPLERIAN": {
                "a": [],
                "e": [],
                "i": [],
                "raan": [],
                "pa": [],
                "true_anomaly": [],
            },
        }
        self.orbit_type = orbit_type
        self.counter = 0
        super(VariableStepHandler, self).__init__()

    def init(self, s0, t):
        print(f"{self.orbit_type} STATES PROPAGATION")
        utils.print_state_header(self.orbit_type)

    def handleStep(self, interpolator):
        currentState = interpolator.getCurrentState()
        wrappedOrbit, keplerian_orbit = orbit_wrapper(currentState)
        utils.update_states(currentState, keplerian_orbit, self.states)

        if self.counter > -1:
            utils.print_state(wrappedOrbit)
            self.counter += 1

    def finish(self, finalState):
        print("this was the last step")
        # utils.plot_states(self.states, save_fig=True, pType=self.orbit_type)


class InitPropagator:
    """
    Initialize Orbit Propagator
    """

    def __init__(self, iOrbit, mass, step_type: str):
        """
        __init__ Initialize InitPropagator class

        Args:
            iOrbit (xxxOrbit): Orbit defined from initial conditions for propagation
            step_type (str): Propagation step type. Must be selected between 'FIXED' or 'VARIABLE'
        """
        self.initialOrbit = iOrbit  # Set the initial orbit
        self.initialState = SpacecraftState(iOrbit, mass)  # Set the initial state
        self.propagation_type = step_type
        print("\nPROPAGATOR STEP TYPE")
        print(f"    {self.propagation_type}")

    def set_integrator(self, iSetup: dict):
        """
        set_integrator Set up a Numerical Integrator

        Args:
            iSetup (dict): Contains keys such as
                iStep (float): A float value for the Numerical integrator's step size
                iType (str): A string representing the type of the step size:
                    FIXED
                    ADAPTIVE
                iSolver (str): A string representing the solver type:
                    STIFF
                    NONSTIFF

        Raises:
            NameError: if the desired 'iSolver' is not found
            NameError: if the desired 'iType' is not found
            NameError: if the desired 'pName' is not found in the 'iType' -> 'FIXED'
            NameError: if the desired 'pName' is not found in the 'iType' -> 'ADAPTIVE'

        Returns:
            ODEIntegrator: Numerical integrator to solve ODE
        """
        # Set up integrator variables
        self.integratorSolver = iSetup.get("iSolver", "NONSTIFF")
        self.integratorType = iSetup.get("iType", "ADAPTIVE")
        self.integratorminStep = iSetup.get("iminStep", 0.001)
        self.integratormaxStep = iSetup.get("imaxStep", 1000.0)
        self.integratordP = iSetup.get("idP", 1.0)
        self.integratornSteps = iSetup.get("inSteps", 1000)

        print(f"INTEGRATOR SETUP")
        print(f"    Solver: {self.integratorSolver}")
        print(f"    Step control: {self.integratorType}")
        print(f"    Step size: {self.propagatorStep} s")
        # Set up the integrator
        if self.integratorSolver == "NONSTIFF":
            if self.integratorType == "FIXED":
                if self.propagatorName == "Euler":
                    self.integrator = EulerIntegrator(self.propagatorStep)
                elif self.propagatorName == "Midpoint":
                    self.integrator = MidpointIntegrator(self.propagatorStep)
                elif self.propagatorName == "ClassicalRungeKutta":
                    self.integrator = ClassicalRungeKuttaIntegrator(self.propagatorStep)
                elif self.propagatorName == "Gill":
                    self.integrator = GillIntegrator(self.propagatorStep)
                elif self.propagatorName == "ThreeEighthes":
                    self.integrator = ThreeEighthesIntegrator(self.propagatorStep)
                elif self.propagatorName == "Luther":
                    self.integrator = LutherIntegrator(self.propagatorStep)
                else:
                    raise NameError("The requested integrator is not defined")
            elif self.integratorType == "ADAPTIVE":
                self.integratorTolerances = NumericalPropagator.tolerances(
                    self.integratordP, self.initialOrbit, self.initialOrbit.getType()
                )
                if self.propagatorName == "HighamHall54":
                    self.integrator = HighamHall54Integrator(
                        self.integratorminStep,
                        self.integratormaxStep,
                        JArray_double.cast_(self.integratorTolerances[0]),
                        JArray_double.cast_(self.integratorTolerances[1]),
                    )
                elif self.propagatorName == "DormandPrince54":
                    self.integrator = DormandPrince54Integrator(
                        self.integratorminStep,
                        self.integratormaxStep,
                        JArray_double.cast_(self.integratorTolerances[0]),
                        JArray_double.cast_(self.integratorTolerances[1]),
                    )
                elif self.propagatorName == "DormandPrince853":
                    self.integrator = DormandPrince853Integrator(
                        self.integratorminStep,
                        self.integratormaxStep,
                        JArray_double.cast_(self.integratorTolerances[0]),
                        JArray_double.cast_(self.integratorTolerances[1]),
                    )
                elif self.propagatorName == "AdamsBashforth":
                    self.integrator = AdamsBashforthIntegrator(
                        self.integratornSteps,
                        self.integratorminStep,
                        self.integratormaxStep,
                        JArray_double.cast_(self.integratorTolerances[0]),
                        JArray_double.cast_(self.integratorTolerances[1]),
                    )
                else:
                    raise NameError("The requested integrator is not defined")
                print(f"    Min. step: {self.integrator.getMinStep()} s")
                print(f"    Max. step: {self.integrator.getMaxStep()} s")
                print(f"    Position tolerance: {self.integratordP} m")
            else:
                raise NameError("The requested integrator type is not defined")
        elif self.integratorSolver == "STIFF":
            if self.integratorType == "FIXED":
                raise NameError(
                    "The requested integrator type 'FIXED' does not exist in 'STIFF' solvers"
                )
            elif self.integratorType == "ADAPTIVE":
                if self.propagatorName == "AdamsMoulton":
                    self.integrator = AdamsMoultonIntegrator(
                        self.integratornSteps,
                        self.integratorminStep,
                        self.integratormaxStep,
                        JArray_double.cast_(self.integratorTolerances[0]),
                        JArray_double.cast_(self.integratorTolerances[1]),
                    )
                elif self.propagatorName == "GraggBulirschStoer":
                    self.integrator = GraggBulirschStoerIntegrator(
                        self.integratorminStep,
                        self.integratormaxStep,
                        JArray_double.cast_(self.integratorTolerances[0]),
                        JArray_double.cast_(self.integratorTolerances[1]),
                    )
                else:
                    raise NameError("The requested integrator is not defined")
            else:
                raise NameError("The requested integrator type is not defined")
        else:
            raise NameError("The requested solver is not defined")
        print(f"    Integrator: {self.integrator.getName()}")

        return self.integrator

    def set_propagation(self, pModel: str, pName: str, **kwargs):
        """
        set_propagation Defines propagator model

        Args:
            pModel (str): A string representing the propagator model. Must be selected between 'NUMERICAL' and 'ANALYTICAL'. Defaults to "NUMERICAL":
                NUMERICAL: Integrators use numerical models to solve the ODE
                ANALYTICAL: Integrators use analytical models to solve the ODE

            pName (str): A string representing the name of pModel to be used. Defaults to "DormandPrince853":
                NUMERICAL
                    FIXED <-- iType
                        Euler
                        Midpoint
                        ClassicalRungeKutta
                        Gill
                        ThreeEighthes
                        Luther
                    ADAPTIVE <-- iType
                        HighamHall54
                        DormandPrince54
                        DormandPrince853
                        GraggBulirschStoer
                        AdamsBashforth
                        AdamsMoulton
                ANALYTICAL
                    Kepler
                    Eckstein-Heschler
                    Brouwer-Lyddane
                    TLE
            **kwargs: Arbitrary keyword arguments:
                iStep (float, optional): A float value for the Numerical integrator step size
                iType (str, optional): A string representing the type of the step size:
                    FIXED
                    ADAPTIVE
                iSolver (str, optional): A string representing the solver type:
                    STIFF
                    NONSTIFF
                inSteps (int, optional): A int indicating number of computing steps (excluding the final one)
        Raises:
            NameError: if the desired 'pModel' is not found
            NameError: if the desired 'pName' is not found
            NameError: if the desired 'step_type' is not found

        Returns:
            xxxPropagator: Propagator model with the desired integrator or analytical model:
            Updates integrator (ODEIntegrator) and propagator (xxxPropagator)
        """
        self.propagatorStep = kwargs.get("iStep", 60.0)
        self.propagatorModel = pModel
        self.propagatorName = pName
        print(f"PROPAGATOR MODEL")
        print(f"    {self.propagatorModel}")
        if (
            self.propagation_type == "FIXED" or self.propagation_type == "VARIABLE"
        ) and self.propagatorModel == "NUMERICAL":
            self.set_integrator(kwargs)
            self.propagator = NumericalPropagator(self.integrator)
            self.propagator.setInitialState(self.initialState)
            self.propagator.setOrbitType(self.initialOrbit.getType())
            self.propagator_orbit_type = self.initialOrbit.getType().toString()
        elif (self.propagation_type == "FIXED") and (
            self.propagatorModel == "ANALYTICAL"
        ):
            if self.propagatorName == "Kepler":
                self.propagator = KeplerianPropagator(self.initialOrbit)
                self.propagator_orbit_type = self.initialOrbit.getType().toString()
            elif self.propagatorName == "Eckstein-Heschler":
                self.propagator = EcksteinHechlerPropagator(
                    self.initialOrbit,
                    Constants.EIGEN5C_EARTH_EQUATORIAL_RADIUS,
                    Constants.EIGEN5C_EARTH_MU,
                    Constants.EIGEN5C_EARTH_C20,
                    Constants.EIGEN5C_EARTH_C30,
                    Constants.EIGEN5C_EARTH_C40,
                    Constants.EIGEN5C_EARTH_C50,
                    Constants.EIGEN5C_EARTH_C60,
                )
                self.propagator_orbit_type = "CARTESIAN"
            elif self.propagatorName == "Brouwer-Lyddane":
                self.propagator = BrouwerLyddanePropagator(
                    self.initialOrbit,
                    Constants.EIGEN5C_EARTH_EQUATORIAL_RADIUS,
                    Constants.EIGEN5C_EARTH_MU,
                    Constants.EIGEN5C_EARTH_C20,
                    Constants.EIGEN5C_EARTH_C30,
                    Constants.EIGEN5C_EARTH_C40,
                    Constants.EIGEN5C_EARTH_C50,
                    Constants.EIGEN5C_EARTH_C60,
                )
                self.propagator_orbit_type = "KEPLERIAN"
            elif self.propagatorName == "TLE":
                self.propagator = TLEPropagator.selectExtrapolator(
                    kwargs.get("tle_str")
                )
                self.propagator_orbit_type = "CARTESIAN"
            else:
                raise NameError("The requested propagator name is not defined")
        else:
            raise NameError("The requested propagator model is not defined")

        # Setup Handler
        self.propagator.getMultiplexer().clear()
        if self.propagation_type == "FIXED":
            self.handler = FixedStepHandler(self.propagator_orbit_type)
            self.propagator.getMultiplexer().add(self.propagatorStep, self.handler)
        elif self.propagation_type == "VARIABLE":
            self.handler = VariableStepHandler(self.propagator_orbit_type)
            self.propagator.getMultiplexer().add(self.handler)
        else:
            raise NameError("The requested propagation step type is not defined")

        self.propagator_states = self.handler.states

        return self.propagator

    def set_model(self, model: str, **kwargs):
        """
        set_model Return a Orekit's PROPAGATOR for the chosen model with the ThirdBody perturbations and external forces

        Args:
            model (str): A string representing the Equation of motion:
                TWO_BODY: Contains only Spacecraft and Central attraction body in the EOM
                N_BODY: Contains TWO_BODY model and other N bodies specified in the EOM
            bodies (list): A list of celestial bodies perturbations:
                EARTH: Planet Earth
                EARTH_MOON: Earth-Moon barycenter
                JUPITER: Planet Jupiter
                MARS: Planet Mars
                MERCURY: Planet
                MOON: Planet Earth's Moon
                NEPTUNE: Planet Neptune
                PLUTO: Planet Pluto
                SATURN: Planet Saturn
                SOLAR_SYSTEM_BARYCENTER: Solar System Barycenter
                SUN: Sun
                URANUS: Planet Uranus
                VENUS: Planet Venus
            forces (list): A list of active and passive forces acting on the spacecraft. Must be selected from below:
                GRAVITY: Gravity Field

        Raises:
            NameError: if the desired 'bodies' is not found
            NameError: if the desired 'model' is not defined
            NameError: if the desired 'force' is not defined

        Returns:
            Propagator: Propagator loaded with the necessary celestial bodies and external forces contributions
        """
        bodies = kwargs.get("bodies", [])
        forces = kwargs.get("forces", [])
        print("DYNAMIC MODEL")
        print(f"    {model}")
        n_bodies = [
            "EARTH",
            "JUPITER",
            "MARS",
            "MERCURY",
            "MOON",
            "NEPTUNE",
            "PLUTO",
            "SATURN",
            "SUN",
            "URANUS",
            "VENUS",
        ]
        non_planets = ["EARTH_MOON", "SOLAR_SYSTEM_BARYCENTER"]
        if model == "TWO_BODY":
            pass
        elif model == "N_BODY":
            if bodies:
                for b in bodies:
                    if b == "ALL":
                        for bb in n_bodies:
                            if bb in non_planets or bb == "EARTH":
                                pass
                            else:
                                body = CelestialBodyFactory.getBody(bb)
                                print(f"        {bb}")
                        self.propagator.addForceModel(ThirdBodyAttraction(body))
                    elif b in (n_bodies + non_planets):
                        body = CelestialBodyFactory.getBody(b)
                        self.propagator.addForceModel(ThirdBodyAttraction(body))
                        print(f"        {b}")
                    else:
                        raise NameError("The requested body is not defined")
        else:
            raise NameError("The requested model is not defined")

        if forces and (self.propagatorModel == "NUMERICAL"):
            print(f"    ACTIVE FORCES")
            for f in forces:
                # Body-centric Frame
                ITRF_frame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
                # Normalized Spherical Harmonics Provider
                degree = 100
                order = 100
                gravityProvider = GravityFieldFactory.getNormalizedProvider(
                    degree, order
                )
                if f == "GRAVITY":
                    self.propagator.addForceModel(
                        HolmesFeatherstoneAttractionModel(
                            ITRF_frame,
                            gravityProvider,
                        )
                    )
                    print(f"        {f}")
                elif f == "SOLID":
                    # Add solid tides
                    solid_tides = SolidTides(
                        ITRF_frame,
                        gravityProvider.getAe(),
                        gravityProvider.getMu(),
                        gravityProvider.getTideSystem(),
                        IERSConventions.IERS_2010,
                        TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True),
                        [CelestialBodyFactory.getSun(), CelestialBodyFactory.getMoon()],
                    )
                    self.propagator.addForceModel(solid_tides)
                    print(f"        {f}")
                elif f == "OCEAN":
                    # Add solid tides
                    ocean_tides = OceanTides(
                        ITRF_frame,
                        gravityProvider.getAe(),
                        gravityProvider.getMu(),
                        degree,
                        order,
                        IERSConventions.IERS_2010,
                        TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True),
                    )
                    self.propagator.addForceModel(ocean_tides)
                    print(f"        {f}")
                elif f == "DRAG":
                    # Define earth body
                    earthBody = OneAxisEllipsoid(
                        Constants.WGS84_EARTH_EQUATORIAL_RADIUS,  # Equatorial radius (meters)
                        Constants.WGS84_EARTH_FLATTENING,  # Flattening
                        ITRF_frame,  # Reference frame
                    )
                    # Atmosphere model
                    atmosphere = NRLMSISE00(
                        CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt"),
                        CelestialBodyFactory.getSun(),
                        earthBody,
                    )
                    # Define spacecraft shape and drag coefficient
                    drag_shape = IsotropicDrag(
                        2.2, 11.0
                    )  # CD = 2.2, Surface Area = 11.0 m²
                    self.propagator.addForceModel(DragForce(atmosphere, drag_shape))
                    print(f"        {f}")
                elif f == "RELATIVITY":
                    relativity = Relativity(Constants.WGS84_EARTH_MU)
                    self.propagator.addForceModel(relativity)
                    print(f"        {f}")
                else:
                    raise NameError(f"The requested force model {f} is not defined")

        return self.propagator, self.propagator_states
