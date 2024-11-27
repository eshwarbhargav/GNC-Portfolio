import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec


def update_states(tmp_state, keplerian_state, states):
    # Initialize
    pv_coords = tmp_state.getPVCoordinates()
    states["TIME"].append(pv_coords.getDate().toString())

    # CARTESIAN
    position = pv_coords.getPosition()
    states["CARTESIAN"]["position"].append(
        [position.getX(), position.getY(), position.getZ()]
    )
    velocity = pv_coords.getVelocity()
    states["CARTESIAN"]["velocity"].append(
        [velocity.getX(), velocity.getY(), velocity.getZ()]
    )
    accleration = pv_coords.getAcceleration()
    states["CARTESIAN"]["acceleration"].append(
        [accleration.getX(), accleration.getY(), accleration.getZ()]
    )

    # KEPLERIAN
    states["KEPLERIAN"]["a"].append(keplerian_state.getA())
    states["KEPLERIAN"]["e"].append(keplerian_state.getE())
    states["KEPLERIAN"]["i"].append(np.rad2deg(keplerian_state.getI()) % 360)
    states["KEPLERIAN"]["raan"].append(
        np.rad2deg(keplerian_state.getRightAscensionOfAscendingNode()) % 360
    )
    states["KEPLERIAN"]["pa"].append(
        np.rad2deg(keplerian_state.getPerigeeArgument()) % 360
    )
    states["KEPLERIAN"]["true_anomaly"].append(
        np.rad2deg(keplerian_state.getTrueAnomaly()) % 360
    )


def print_state_header(orbit_type):
    if orbit_type == "KEPLERIAN":
        header = (
            "          date                        a          e"
            "         i          ω          Ω          ν     "
        )
    elif orbit_type == "CARTESIAN":
        header = (
            "          date                        x          y"
            "          z            xDot            yDot            zDot"
        )
    elif orbit_type == "EQUINOCTIAL":
        header = (
            "          date                        a          ex"
            "          ey            hx            hy            lv"
        )
    elif orbit_type == "CIRCULAR":
        header = (
            "          date                        a          ex"
            "          ey            i            Ω            αv"
        )
    else:
        raise NameError("The desired orbit type header is not defined")
    return print(header)


def print_state(orbit):
    if orbit.getType().toString() == "KEPLERIAN":
        body = "{:<30} {:12.3f} {:10.8f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
            orbit.getDate().toString(),
            orbit.getA(),
            orbit.getE(),
            np.rad2deg(orbit.getI()) % 360,
            np.rad2deg(orbit.getPerigeeArgument()) % 360,
            np.rad2deg(orbit.getRightAscensionOfAscendingNode()) % 360,
            np.rad2deg(orbit.getTrueAnomaly()) % 360,
        )
    elif orbit.getType().toString() == "CARTESIAN":
        body = "{:<30} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f}".format(
            orbit.getDate().toString(),
            orbit.getPVCoordinates().getPosition().getX(),
            orbit.getPVCoordinates().getPosition().getY(),
            orbit.getPVCoordinates().getPosition().getZ(),
            orbit.getPVCoordinates().getVelocity().getX(),
            orbit.getPVCoordinates().getVelocity().getY(),
            orbit.getPVCoordinates().getVelocity().getZ(),
        )
    elif orbit.getType().toString() == "EQUINOCTIAL":
        body = "{:<30} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f}".format(
            orbit.getDate().toString(),
            orbit.getA(),
            np.rad2deg(orbit.getEquinoctialEx()) % 360,
            np.rad2deg(orbit.getEquinoctialEy()) % 360,
            np.rad2deg(orbit.getHx()) % 360,
            np.rad2deg(orbit.getHy()) % 360,
            np.rad2deg(orbit.getLv()) % 360,
        )
    elif orbit.getType().toString() == "CIRCULAR":
        body = "{:<30} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f}".format(
            orbit.getDate().toString(),
            orbit.getA(),
            np.rad2deg(orbit.getCircularEx()) % 360,
            np.rad2deg(orbit.getCircularEy()) % 360,
            np.rad2deg(orbit.getI()) % 360,
            np.rad2deg(orbit.getRightAscensionOfAscendingNode()) % 360,
            np.rad2deg(orbit.getAlphaV()) % 360,
        )
    else:
        raise NameError("The desired orbit's state is not defined")
    return print(body)


def plot_orbit_3D(Pstart, Pend, Pstates):
    # Create a figure for the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the orbit path in 3D
    ax.plot(
        np.divide(Pstates[0], 1000),
        np.divide(Pstates[1], 1000),
        np.divide(Pstates[2], 1000),
        color="b",
        label="Orbit Path",
    )

    # Mark the starting point
    ax.scatter(
        Pstart.getX() / 1000,
        Pstart.getY() / 1000,
        Pstart.getZ() / 1000,
        color="r",
        marker="o",
        label="Start Position",
    )

    # Mark the final position
    ax.scatter(
        Pend.getX() / 1000,
        Pend.getY() / 1000,
        Pend.getZ() / 1000,
        color="g",
        marker="x",
        label="End Position",
    )

    # Labels for axes
    ax.set_xlabel("X Position (km)")
    ax.set_ylabel("Y Position (km)")
    ax.set_zlabel("Z Position (km)")

    # Title and legend
    ax.set_title("3D Orbit Plot")
    ax.legend()
    plt.show()


def states_to_cartesian(states):
    # CARTESIAN
    # Position
    pos = np.array(states["CARTESIAN"]["position"])
    # Velocity
    vel = np.array(states["CARTESIAN"]["velocity"])
    # Acceleration
    acc = np.array(states["CARTESIAN"]["acceleration"])

    return pos, vel, acc


def states_to_keplerian(states):
    # KEPLERIAN
    a_axis = np.array(states["KEPLERIAN"]["a"])
    ecc = np.array(states["KEPLERIAN"]["e"])
    inc = np.array(states["KEPLERIAN"]["i"])
    Omega = np.array(states["KEPLERIAN"]["raan"])
    omega = np.array(states["KEPLERIAN"]["pa"])
    thetha = np.array(states["KEPLERIAN"]["true_anomaly"])

    return a_axis, ecc, inc, Omega, omega, thetha


# def plot_states(states, save_fig: bool = False, pType: str = "NUMERICAL"):
#     # Time
#     time = np.arange(len(states["TIME"]))
#     # Cartesian
#     position, velocity, acceleration = states_to_cartesian(states)

#     # Plotting
#     fig, ax = plt.subplots(3, sharex=True)
#     # fig.suptitle("Propagated states")

#     ax[0].plot(time, position[:, 0], label=r"$r_{x}$")
#     ax[0].plot(time, position[:, 1], label=r"$r_{y}$")
#     ax[0].plot(time, position[:, 2], label=r"$r_{z}$")
#     ax[0].set_ylabel(r"$\vec{r}\ [m]$")
#     ax[0].legend()

#     ax[1].plot(time, velocity[:, 0], label=r"$v_{x}$")
#     ax[1].plot(time, velocity[:, 1], label=r"$v_{y}$")
#     ax[1].plot(time, velocity[:, 2], label=r"$v_{z}$")
#     ax[1].set_ylabel(r"$\vec{v}\ [m/s]$")
#     ax[1].legend()

#     ax[2].plot(time, acceleration[:, 0], label=r"$a_{x}$")
#     ax[2].plot(time, acceleration[:, 1], label=r"$a_{y}$")
#     ax[2].plot(time, acceleration[:, 2], label=r"$a_{z}$")
#     ax[2].set_xlabel("Integration time [s]")
#     ax[2].set_ylabel(r"$\vec{a}\ [m/s^{2}]$")
#     ax[2].legend()

#     plt.tight_layout()

#     if save_fig is True:
#         plt.savefig(
#             f"outputs/cartesian_coords_{pType}.png",
#             dpi=300,
#             format="png",
#             bbox_inches="tight",
#             pad_inches=0.1,
#         )

# semimajor_axis, eccentricity, inclination, raan, perigee_argument, true_anomaly = (
#     states_to_keplerian(states)
# )

# # Plotting
# fig = plt.figure()

# gs = gridspec.GridSpec(4, 2)
# # Add subplots
# ax1 = fig.add_subplot(gs[0, 0])  # First row, first column
# ax2 = fig.add_subplot(gs[0, 1])  # First row, second column
# ax3 = fig.add_subplot(gs[1, 0])  # Second row, first column
# ax4 = fig.add_subplot(gs[1, 1])  # Second row, second column
# ax5 = fig.add_subplot(gs[2, 0])  # Third row, first column
# ax6 = fig.add_subplot(gs[2, 1])  # Third row, second column
# ax7 = fig.add_subplot(gs[3, :])  # Last row, spans both columns

# ax1.plot(time, semimajor_axis, linewidth=1, label="Semimajor axis")
# # ax1.set_xlabel("Integration time [s]")
# ax1.set_ylabel(r"$a$")

# ax2.plot(time, eccentricity, linewidth=1, label="Eccentricity")
# # ax[0, 1].set_xlabel("Integration time [s]")
# ax2.set_ylabel(r"$e$")

# ax3.plot(time, inclination, linewidth=1, label="Inclination")
# # ax[1, 0].set_xlabel("Integration time [s]")
# ax3.set_ylabel(r"$i$")

# ax4.plot(time, raan, linewidth=1, label="RAAN")
# # ax[1, 1].set_xlabel("Integration time [s]")
# ax4.set_ylabel(r"$\Omega$")

# ax5.plot(time, perigee_argument, linewidth=1, label="Perigee Argument")
# # ax[2, 0].set_xlabel("Integration time [s]")
# ax5.set_ylabel(r"$\omega$")

# ax6.plot(time, true_anomaly, linewidth=1, label="True Anomaly")
# # ax[2, 1].set_xlabel("Integration time [s]")
# ax6.set_ylabel(r"$\theta$")

# ax7.plot(
#     time,
#     (perigee_argument + true_anomaly) % 360,
#     linewidth=1,
#     label="Perigee + Anomaly",
# )
# # ax[2, 1].set_xlabel("Integration time [s]")
# ax7.set_ylabel(r"$\alpha_{v} = \omega + \theta$")

# plt.tight_layout()

# if save_fig is True:
#     plt.savefig(
#         f"outputs/keplerian_coords_{pType}.png",
#         dpi=300,
#         format="png",
#         bbox_inches="tight",
#         pad_inches=0.1,
#     )


def plot_error(estimated, actual, save_fig: bool = True):
    # Time
    time = np.arange(len(actual["TIME"]))

    # Cartesian
    position_est, velocity_est, acceleration_est = states_to_cartesian(estimated)
    position_actual, velocity_actual, acceleration_actual = states_to_cartesian(actual)

    pos_norm_est = np.linalg.norm(position_est, axis=1)
    pos_norm_actual = np.linalg.norm(position_actual, axis=1)

    vel_norm_est = np.linalg.norm(velocity_est, axis=1)
    vel_norm_actual = np.linalg.norm(velocity_actual, axis=1)

    acc_norm_est = np.linalg.norm(acceleration_est, axis=1)
    acc_norm_actual = np.linalg.norm(acceleration_actual, axis=1)

    # Plotting
    fig, ax = plt.subplots(3, sharex=True)
    fig.suptitle("Cartesian propagation")

    ax[0].plot(time, pos_norm_est, label="Numerical", linestyle="--", linewidth=1)
    ax[0].plot(time, pos_norm_actual, label="Analytical", linestyle="-", linewidth=1)
    ax[0].set_ylabel(r"$|\vec{r}|\ [m]$")
    # ax[0].legend()

    ax[1].plot(time, vel_norm_est, label="Numerical", linestyle="--", linewidth=1)
    ax[1].plot(time, vel_norm_actual, label="Analytical", linestyle="-", linewidth=1)
    ax[1].set_ylabel(r"$|\vec{v}|\ [m/s]$")
    # ax[1].legend()

    ax[2].plot(time, acc_norm_est, label="Numerical", linestyle="--", linewidth=1)
    ax[2].plot(time, acc_norm_actual, label="Analytical", linestyle="-", linewidth=1)
    ax[2].set_xlabel("Integration time [s]")
    ax[2].set_ylabel(r"$|\vec{a}|\ [m/s^{2}]$")
    # ax[2].legend()

    fig.legend(
        labels=["Numerical", "Analytical"],
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.05),
    )
    plt.tight_layout()

    if save_fig is True:
        plt.savefig(
            f"outputs/cartesian_coords.png",
            dpi=300,
            format="png",
            bbox_inches="tight",
            pad_inches=0.1,
        )

    # Compute error
    delta_pos = position_est - position_actual
    delta_vel = velocity_est - velocity_actual
    delta_acc = acceleration_est - acceleration_actual

    pos_err = np.linalg.norm(delta_pos, axis=1)
    vel_err = np.linalg.norm(delta_vel, axis=1)
    acc_err = np.linalg.norm(delta_acc, axis=1)

    # Plotting
    fig, ax = plt.subplots(3, sharex=True)
    fig.suptitle("Cartesian error propagation")

    ax[0].plot(time, pos_err, label="Position")
    ax[0].set_ylabel(r"$|\Delta\vec{r}|\ [m]$")
    ax[0].legend()

    ax[1].plot(time, vel_err, label="Velocity")
    ax[1].set_ylabel(r"$|\Delta\vec{v}|\ [m/s]$")
    ax[1].legend()

    ax[2].plot(time, acc_err, label="Accleration")
    ax[2].set_xlabel(r"Integration time $[s]$")
    ax[2].set_ylabel(r"$|\Delta\vec{a}|\ [m/s^{2}]$")
    ax[2].legend()

    plt.tight_layout()

    if save_fig is True:
        plt.savefig(
            f"outputs/cartesian_error_coords.png",
            dpi=300,
            format="png",
            bbox_inches="tight",
            pad_inches=0.1,
        )

    (
        semimajor_axis_est,
        eccentricity_est,
        inclination_est,
        raan_est,
        perigee_argument_est,
        true_anomaly_est,
    ) = states_to_keplerian(estimated)

    (
        semimajor_axis_actual,
        eccentricity_actual,
        inclination_actual,
        raan_actual,
        perigee_argument_actual,
        true_anomaly_actual,
    ) = states_to_keplerian(actual)

    # Plotting
    fig = plt.figure()
    fig.suptitle("Keplerian propagation")

    gs = gridspec.GridSpec(4, 2)
    # Add subplots
    ax1 = fig.add_subplot(gs[0, 0])  # First row, first column
    ax2 = fig.add_subplot(gs[0, 1])  # First row, second column
    ax3 = fig.add_subplot(gs[1, 0])  # Second row, first column
    ax4 = fig.add_subplot(gs[1, 1])  # Second row, second column
    ax5 = fig.add_subplot(gs[2, 0])  # Third row, first column
    ax6 = fig.add_subplot(gs[2, 1])  # Third row, second column
    ax7 = fig.add_subplot(gs[3, :])  # Last row, spans both columns

    ax1.plot(time, semimajor_axis_est, label="Numerical", linestyle="--", linewidth=1)
    ax1.plot(
        time, semimajor_axis_actual, label="Analytical", linestyle="-", linewidth=1
    )
    # ax1.set_xlabel("Integration time [s]")
    ax1.set_ylabel(r"$a$")

    ax2.plot(time, eccentricity_est, label="Numerical", linestyle="--", linewidth=1)
    ax2.plot(time, eccentricity_actual, label="Analytical", linestyle="-", linewidth=1)
    # ax[0, 1].set_xlabel("Integration time [s]")
    ax2.set_ylabel(r"$e$")

    ax3.plot(time, inclination_est, label="Numerical", linestyle="--", linewidth=1)
    ax3.plot(time, inclination_actual, label="Analytical", linestyle="-", linewidth=1)
    # ax[1, 0].set_xlabel("Integration time [s]")
    ax3.set_ylabel(r"$i$")

    ax4.plot(time, raan_est, label="Numerical", linestyle="--", linewidth=1)
    ax4.plot(time, raan_actual, label="Analytical", linestyle="-", linewidth=1)
    # ax[1, 1].set_xlabel("Integration time [s]")
    ax4.set_ylabel(r"$\Omega$")

    ax5.plot(time, perigee_argument_est, label="Numerical", linestyle="--", linewidth=1)
    ax5.plot(
        time, perigee_argument_actual, label="Analytical", linestyle="-", linewidth=1
    )
    # ax[2, 0].set_xlabel("Integration time [s]")
    ax5.set_ylabel(r"$\omega$")

    ax6.plot(time, true_anomaly_est, label="Numerical", linestyle="--", linewidth=1)
    ax6.plot(time, true_anomaly_actual, label="Analytical", linestyle="-", linewidth=1)
    # ax[2, 1].set_xlabel("Integration time [s]")
    ax6.set_ylabel(r"$\theta$")

    ax7.plot(
        time,
        (perigee_argument_est + true_anomaly_est) % 360,
        linewidth=1,
        label="Numerical",
        linestyle="--",
    )
    ax7.plot(
        time,
        (perigee_argument_actual + true_anomaly_actual) % 360,
        linewidth=1,
        label="Analytical",
        linestyle="-",
    )
    ax7.set_xlabel(r"Integration time $[s]$")
    ax7.set_ylabel(r"$\alpha_{v} = \omega + \theta$")

    fig.legend(
        labels=["Numerical", "Analytical"],
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.05),
    )
    plt.tight_layout()

    if save_fig is True:
        plt.savefig(
            f"outputs/keplerian_coords.png",
            dpi=300,
            format="png",
            bbox_inches="tight",
            pad_inches=0.1,
        )

    # Compute error
    semimajor_axis_err = np.subtract(semimajor_axis_est, semimajor_axis_actual)
    eccentricity_err = np.subtract(eccentricity_est, eccentricity_actual)
    inclination_err = np.subtract(inclination_est, inclination_actual)
    raan_err = np.subtract(raan_est, raan_actual)
    perigee_argument_err = np.subtract(perigee_argument_est, perigee_argument_actual)
    true_anomaly_err = np.subtract(true_anomaly_est, true_anomaly_actual)
    pa_ta_err = np.subtract(
        (perigee_argument_est + true_anomaly_est) % 360,
        (perigee_argument_actual + true_anomaly_actual) % 360,
    )

    # Plotting
    fig = plt.figure()
    fig.suptitle("Keplerian error propagation")

    gs = gridspec.GridSpec(4, 2)
    # Add subplots
    ax1 = fig.add_subplot(gs[0, 0])  # First row, first column
    ax2 = fig.add_subplot(gs[0, 1])  # First row, second column
    ax3 = fig.add_subplot(gs[1, 0])  # Second row, first column
    ax4 = fig.add_subplot(gs[1, 1])  # Second row, second column
    ax5 = fig.add_subplot(gs[2, 0])  # Third row, first column
    ax6 = fig.add_subplot(gs[2, 1])  # Third row, second column
    ax7 = fig.add_subplot(gs[3, :])  # Last row, spans both columns

    ax1.plot(time, semimajor_axis_err, linewidth=1, label="Semimajor axis")
    # ax[0, 0].set_xlabel("Integration time [s]")
    ax1.set_ylabel(r"$\Delta{a}$")

    ax2.plot(time, eccentricity_err, linewidth=1, label="Eccentricity")
    # ax[0, 1].set_xlabel("Integration time [s]")
    ax2.set_ylabel(r"$\Delta{e}$")

    ax3.plot(time, inclination_err, linewidth=1, label="Inclination")
    # ax[1, 0].set_xlabel("Integration time [s]")
    ax3.set_ylabel(r"$\Delta{i}$")

    ax4.plot(time, raan_err, linewidth=1, label="RAAN")
    # ax[1, 1].set_xlabel("Integration time [s]")
    ax4.set_ylabel(r"$\Delta\Omega$")

    ax5.plot(time, perigee_argument_err, linewidth=1, label="Perigee Argument")
    # ax[2, 0].set_xlabel("Integration time [s]")
    ax5.set_ylabel(r"$\Delta\omega$")

    ax6.plot(time, true_anomaly_err, linewidth=1, label="True Anomaly")
    # ax[2, 1].set_xlabel("Integration time [s]")
    ax6.set_ylabel(r"$\Delta\theta$")

    ax7.plot(
        time, (pa_ta_err + 180) % 360 - 180, linewidth=1, label="Perigee + Anomaly"
    )
    ax7.set_xlabel(r"Integration time $[s]$")
    ax7.set_ylabel(r"$\Delta\alpha_{v} = \omega + \theta$")

    plt.tight_layout()

    if save_fig is True:
        plt.savefig(
            f"outputs/keplerian_error_coords.png",
            dpi=300,
            format="png",
            bbox_inches="tight",
            pad_inches=0.1,
        )
