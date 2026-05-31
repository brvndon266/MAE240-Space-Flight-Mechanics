import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# MAE 240 - Homework 4 - Problem 3
mu = 1.0

# Initial Conditions
r0 = np.array([1.0, 0.0, 0.0])
v0 = np.array([0.0, 1.0, 0.0])

# Disturbing acceleration magnitudes
a_values = [0.001, 0.01, 0.1, 1.0]

# Circular orbit period
T = 2 * np.pi

# Simulate for 10 revolutions
t_final = 10 * T

# Time vector
t_eval = np.linspace(0, t_final, 3000)

# Equations of Motion
def eom(t, state, mu, a_d):

    r = state[0:3]
    v = state[3:6]

    r_norm = np.linalg.norm(r)

    # Central gravity acceleration
    a_gravity = -mu * r / r_norm**3

    # Total acceleration
    a_total = a_gravity + a_d

    return np.hstack((v, a_total))

# Cartesian State -> Orbital Elements
def cartesian_to_elements(r, v, mu):

    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)

    # Angular Momentum
    h_vec = np.cross(r, v)
    h = np.linalg.norm(h_vec)

    # Inclination
    i = np.arccos(h_vec[2] / h)

    # Node Vector
    k_hat = np.array([0.0, 0.0, 1.0])

    n_vec = np.cross(k_hat, h_vec)
    n = np.linalg.norm(n_vec)

    # Eccentricity Vector
    e_vec = np.cross(v, h_vec) / mu - r / r_norm
    e = np.linalg.norm(e_vec)

    # Specific Orbital Energy
    energy = v_norm**2 / 2 - mu / r_norm

    # Semi-major Axis
    if abs(energy) > 1e-12:
        a = -mu / (2 * energy)
    else:
        a = np.nan

    # RAAN
    if n > 1e-12:

        Omega = np.arccos(n_vec[0] / n)

        if n_vec[1] < 0:
            Omega = 2 * np.pi - Omega

    else:
        Omega = 0.0

    # Argument of Periapsis
    if n > 1e-12 and e > 1e-12:

        omega = np.arccos(
            np.dot(n_vec, e_vec) / (n * e)
        )

        if e_vec[2] < 0:
            omega = 2 * np.pi - omega

    else:
        omega = 0.0

    # True Anomaly
    if e > 1e-12:

        f = np.arccos(
            np.dot(e_vec, r) / (e * r_norm)
        )

        if np.dot(r, v) < 0:
            f = 2 * np.pi - f

    else:
        f = 0.0

    # Mean Anomaly
    if e < 1 and e > 1e-12:

        E = 2 * np.arctan2(
            np.sqrt(1 - e) * np.sin(f / 2),
            np.sqrt(1 + e) * np.cos(f / 2)
        )

        M = E - e * np.sin(E)

        M = M % (2 * np.pi)

    else:
        M = 0.0

    return a, e, i, Omega, omega, M

# Run Simulation Cases
for accel_mag in a_values:

    print("\n======================================")
    print(f"Running case: a_d = [{accel_mag}, 0, 0]")
    print("======================================")

    # Disturbing acceleration
    a_d = np.array([accel_mag, 0.0, 0.0])

    # Initial state vector
    state0 = np.hstack((r0, v0))

    # Numerical integration
    sol = solve_ivp(
        eom,
        [0, t_final],
        state0,
        t_eval=t_eval,
        args=(mu, a_d),
        rtol=1e-10,
        atol=1e-12
    )

    t = sol.t
    X = sol.y.T

    r_hist = X[:, 0:3]
    v_hist = X[:, 3:6]

    # Compute Orbital Elements
    elements = np.zeros((len(t), 6))

    for k in range(len(t)):

        elements[k, :] = cartesian_to_elements(
            r_hist[k, :],
            v_hist[k, :],
            mu
        )

    a_hist = elements[:, 0]
    e_hist = elements[:, 1]
    i_hist = np.rad2deg(elements[:, 2])
    Omega_hist = np.rad2deg(elements[:, 3])
    omega_hist = np.rad2deg(elements[:, 4])
    M_hist = np.rad2deg(elements[:, 5])

    r_mag = np.linalg.norm(r_hist, axis=1)
    v_mag = np.linalg.norm(v_hist, axis=1)

    # Figure 1: Trajectory
    plt.figure(figsize=(7, 7))

    plt.plot(
        r_hist[:, 0],
        r_hist[:, 1],
        linewidth=1.5,
        label='Trajectory'
    )

    plt.scatter(
        r_hist[0, 0],
        r_hist[0, 1],
        marker='o',
        s=50,
        label='Start'
    )

    plt.scatter(
        0,
        0,
        marker='x',
        s=120,
        label='Central Body'
    )

    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)

    plt.title(
        f'Trajectory for $a_d = [{accel_mag}, 0, 0]$',
        fontsize=16
    )

    plt.axis('equal')

    plt.grid(True)
    plt.legend(fontsize=10)

    plt.tight_layout()

    plt.show()

    # Figure 2: Orbital Elements vs Time
    fig, axs = plt.subplots(
        4,
        2,
        figsize=(16, 14),
        sharex=True
    )

    # Main figure title
    fig.suptitle(
        f'Orbital Elements vs Time for '
        f'$a_d = [{accel_mag}, 0, 0]$',
        fontsize=18,
        y=0.98
    )

    # LEFT COLUMN
    # Semi-major Axis
    
    axs[0, 0].plot(
        t,
        a_hist,
        linewidth=1.5
    )

    axs[0, 0].set_title('Semi-major Axis')
    axs[0, 0].set_ylabel('a')

    # Eccentricity
    axs[1, 0].plot(
        t,
        e_hist,
        linewidth=1.5
    )

    axs[1, 0].set_title('Eccentricity')
    axs[1, 0].set_ylabel('e')

    # Inclination
    axs[2, 0].plot(
        t,
        i_hist,
        linewidth=1.5
    )

    axs[2, 0].set_title('Inclination')
    axs[2, 0].set_ylabel('i [deg]')

    # RAAN
    axs[3, 0].plot(
        t,
        Omega_hist,
        linewidth=1.5
    )

    axs[3, 0].set_title('RAAN')
    axs[3, 0].set_ylabel('Ω [deg]')
    axs[3, 0].set_xlabel('Time')

    # RIGHT COLUMN
    # Argument of Periapsis
    axs[0, 1].plot(
        t,
        omega_hist,
        linewidth=1.5
    )

    axs[0, 1].set_title('Argument of Periapsis')
    axs[0, 1].set_ylabel('ω [deg]')

    # Mean Anomaly
    axs[1, 1].plot(
        t,
        M_hist,
        linewidth=1.5
    )

    axs[1, 1].set_title('Mean Anomaly')
    axs[1, 1].set_ylabel('M [deg]')

    # Radius Magnitude
    axs[2, 1].plot(
        t,
        r_mag,
        linewidth=1.5
    )

    axs[2, 1].set_title('Radius Magnitude')
    axs[2, 1].set_ylabel('|r|')

    # Velocity Magnitude
    axs[3, 1].plot(
        t,
        v_mag,
        linewidth=1.5
    )

    axs[3, 1].set_title('Velocity Magnitude')
    axs[3, 1].set_ylabel('|v|')
    axs[3, 1].set_xlabel('Time')

    # Formatting
    for row in axs:

        for ax in row:

            ax.grid(True)

            ax.tick_params(
                axis='both',
                labelsize=10
            )

    # Adjust subplot spacing
    fig.subplots_adjust(
        hspace=0.55,
        wspace=0.30,
        top=0.88
    )

    plt.show()