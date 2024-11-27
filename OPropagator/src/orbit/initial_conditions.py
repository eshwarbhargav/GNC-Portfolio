# Standard library
from datetime import datetime, timezone
from math import radians, sin, cos, tan, pi

# Standards
from org.orekit.frames import FramesFactory, Predefined
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import PVCoordinates

# Hipparchus library
from org.hipparchus.geometry.euclidean.threed import Vector3D

# File-specific
from org.orekit.orbits import (
    KeplerianOrbit,
    CartesianOrbit,
    EquinoctialOrbit,
    CircularOrbit,
    PositionAngleType,
    OrbitType,
)
from src.orbit import utils


def transform_orbit(currentOrbit):
    """
    Orbit transformation: Keplerian-to-Cartesian, Cartesian-to-Keplerian
    """
    if currentOrbit.getType().toString() == "KEPLERIAN":
        transformed_orbit = CartesianOrbit(
            OrbitType.CARTESIAN.convertType(currentOrbit)
        )
    elif currentOrbit.getType().toString() == "CARTESIAN":
        transformed_orbit = KeplerianOrbit(
            OrbitType.KEPLERIAN.convertType(currentOrbit)
        )
    else:
        transformed_orbit = None

    return transformed_orbit


def orbit_wrapper(currentState):
    currentOrbit = currentState.getOrbit()
    keplerian_orbit = KeplerianOrbit(currentOrbit)
    orbit_type = currentOrbit.getType().toString()
    if orbit_type == "KEPLERIAN":
        wrappedOrbit = KeplerianOrbit(currentOrbit)
    elif orbit_type == "CARTESIAN":
        wrappedOrbit = CartesianOrbit(currentOrbit)
    elif orbit_type == "EQUINOCTIAL":
        wrappedOrbit = EquinoctialOrbit(currentOrbit)
    elif orbit_type == "CIRCULAR":
        wrappedOrbit = CircularOrbit(currentOrbit)
    else:
        raise NameError("The requested orbit type is not defined")

    return wrappedOrbit, keplerian_orbit


class InitOrbit:
    """
    Initialize Orbit Setup
    """

    def __init__(self, R_i: float, mu_i: float):
        """
        __init__ Initialize InitOrbit class

        Args:
            R_i (float): Radius if the body from it's center
            mu_i (float): Attraction coefficient of the primary body
        """
        # Propagator setup
        print("\nPROPAGATOR SETUP")
        self.R_i = R_i
        self.mu_i = mu_i

    def set_reference_frame(self, frame_type: str = "EME2000"):
        """
        Return a propagator's reference frame from OREKIT library. Default: EME2000

                Parameters:
                        frame (str): A string representing the reference frame
                Options:
                        ECLIPTIC_CONVENTIONS_2010: Ecliptic framen IERS 2010 conventions
                        ICRF: International Celestial Reference Frame
                        GCRF: Geocentric Celestial Reference Frame
                        ITRF_EQUINOX_CONV_2010_SIMPLE_EOP: Equinox-based ITRF, IERS 2010 conventions with simple EOP interpolation.
                        EME2000: Earth Mean Equator and Equinox at J2000 Frame

                Returns:
                        ref_frame (Frame): Extracted reference frame from OREKIT's FramesFactory
        """
        # Define frame
        if frame_type in [
            "ECLIPTIC_CONVENTIONS_2010",
            "ICRF",
            "GCRF",
            "ITRF_EQUINOX_CONV_2010_SIMPLE_EOP",
            "EME2000",
        ]:
            self.frame = FramesFactory.getFrame(Predefined.valueOf(frame_type))
        else:
            raise NameError("The requested frame_type is not defined")

        print("\nREFERENCE FRAME")
        print(f"    {self.frame}")

        return self.frame

    def set_timescale(self, timezone: str = "UTC"):
        """
        set_timescale Returns the timescale.

        Args:
            timezone (str, optional): A string representing timescale. Defaults to "UTC":
                UTC: Coordinated Universal Time
                TAI: International Atomic Time
                TDT: Terrestrial Dynamic Time
                TT: Terrestrial Time
                TDB: Barycentric Dynamic Time
                GPS: Global Positioning System

        Raises:
            NameError: if the argument string doesn't match the available ones.

        Returns:
            timescale (TimeScale): A string with a timescale abbrevation
        """

        # Define timescale
        if timezone == "UTC":
            self.timescale = TimeScalesFactory.getUTC()
        elif timezone == "TAI":
            self.timescale = TimeScalesFactory.getTAI()
        elif timezone == "TDT" or timezone == "TT":
            self.timescale = TimeScalesFactory.getTT()
        elif timezone == "TDB":
            self.timescale = TimeScalesFactory.getTDB()
        elif timezone == "GPS":
            self.timescale = TimeScalesFactory.getGPS()
        else:
            raise NameError("The requested time coordinate is not defined")

        return self.timescale

    def set_epoch(
        self,
        year: int = datetime.now(timezone.utc).year,
        month: int = datetime.now(timezone.utc).month,
        day: int = datetime.now(timezone.utc).day,
        hour: int = datetime.now(timezone.utc).hour,
        minute: int = datetime.now(timezone.utc).minute,
        second: int = datetime.now(timezone.utc).second,
        microsecond: float = datetime.now(timezone.utc).microsecond,
    ) -> AbsoluteDate:
        """
        Return an AbsoluteDate of Epoch and timescale. Defaults: Current DDMMYY:hhmmss in UTC.

                Parameters:
                        year: Epoch year
                        month: Epoch month
                        day: Epoch day
                        hour: Epoch hour
                        minute: Epoch minute
                        second: Epoch second
                        microsecond: Epoch microsecond

                Returns:
                        epoch_date (AbsoluteDate): Transformed epoch in OREKIT's AbsoluteDate format
        """

        # Define epoch/initialDate
        self.epoch = AbsoluteDate(
            year, month, day, hour, minute, second + microsecond / 1e6, self.timescale
        )
        print("EPOCH")
        print(f"    {self.epoch.toString()} {self.timescale.toString()}")

        return self.epoch, self.timescale

    def set_orbit(self, orbit_type: str = "KEPLERIAN", **kwargs):
        """
        set_orbit Returns an Orbit of desired 'orbit_type'

        Args:
            orbit_type (str, optional): Must be selected between. Defaults to 'KEPLERIAN':
                    KEPLERIAN: Keplerian orbital parameters (i.e., a, e, i, ω, Ω, ν),
                    CARTESIAN: Cartesian orbital parameters (i.e., x, y, z, xDot, yDot, zDot),
                    EQUINOCTIAL: Equinoctial orbital parameters that support both circular and equatorial orbits (i.e., a, ex = e cos(ω + Ω), ey = e sin(ω + Ω), hx = tan(i/2) cos(Ω), hy = tan(i/2) sin(Ω), lv = v + ω + Ω),
                    CIRCULAR: Circular orbital parameters (i.e., a, ex = e cos(ω), ey = e sin(ω), i, Ω, α_{v} = v + ω).
            anomaly_type (str, optional): Sets which type of anomaly we use
            tle_string (str, optional): TLE string class

        Raises:
            NameError: if the argument string doesn't match the available ones.

        Returns:
            xxxOrbit: Orbit made up of defined orbital elements either in 'KEPLERIAN' or 'CARTESIAN'
        """

        # Get input parameters
        print("ORBITAL PARAMETERS (user-input) <- Press enter for default values")

        if orbit_type == "TLE":
            tle = kwargs.get("tle_str")
            mean_motion = tle.getMeanMotion()
            e = tle.getE()
            i = tle.getI()
            pa = tle.getPerigeeArgument()
            raan = tle.getRaan()
            anomaly = tle.getMeanAnomaly()
            a = (self.mu_i / (mean_motion**2)) ** (1 / 3)
            anomaly_type = kwargs.get("anomaly_type", "MEAN")
            orbit_type = "KEPLERIAN"
        else:
            # rp = self.R_i + float(input("    Perigee radius (in m): ").strip() or 400000)
            # ra = self.R_i + float(input("    Apogee radius (in m): ").strip() or 3789068.5)
            # i = radians(float(input("    Inclination (in deg): ").strip() or 97.6))
            # pa = radians(float(input("    Argument of perigee (in deg): ").strip() or 90.0))
            # raan = radians(float(input("    Right Ascension of Ascending Node (in deg): ").strip() or 180.0))
            # anomaly = radians(float(input("    Anomaly (in deg): ").strip() or 10.0))

            rp = self.R_i + float(400000)
            ra = self.R_i + float(3789068.5)
            i = radians(float(97.6))
            pa = radians(float(90.0))
            raan = radians(float(180.0))
            anomaly = radians(float(10.0))

            a = (rp + ra) / 2.0  # Semimajor Axis (m)
            e = 1.0 - (rp / a)  # Eccentricity

            anomaly_type = kwargs.get("anomaly_type", "TRUE")

        if e < 0 + 1e-3:
            orbit_type = "CIRCULAR"
            if i < 0 + 1e-3:
                orbit_type = "EQUINOCTIAL"

        if anomaly_type == "MEAN":
            a_type = PositionAngleType.MEAN
        elif anomaly_type == "ECCENTRIC":
            a_type = PositionAngleType.ECCENTRIC
        else:
            a_type = PositionAngleType.TRUE

        self.orbit = KeplerianOrbit(
            a,
            e,
            i,
            pa,
            raan,
            anomaly,
            a_type,  # Sets which type of anomaly we use
            self.frame,  # The frame in which the parameters are defined (must be a pseudo-inertial frame)
            self.epoch,  # Sets the epoch of the orbital parameters
            self.mu_i,  # Sets the central attraction coefficient (m³/s²)
        )
        # Set 'orbit' from input parameters
        if orbit_type == "KEPLERIAN":
            self.orbit
        elif orbit_type == "CARTESIAN":
            self.orbit = CartesianOrbit(self.orbit)
        elif orbit_type == "EQUINOCTIAL":
            self.orbit = EquinoctialOrbit(self.orbit)
        elif orbit_type == "CIRCULAR":
            self.orbit = CircularOrbit(self.orbit)
        else:
            raise NameError("The requested orbit type is not defined")

        print(f"{orbit_type} COORDINATES (initial-orbit)")
        utils.print_state_header(orbit_type)
        utils.print_state(self.orbit)

        return self.orbit
