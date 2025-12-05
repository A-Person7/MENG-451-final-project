#!/usr/bin/env python3

"""
    simulation.py

    A simulation of inverted pendulum dynamics to determine what system
    parameters are required for a physically plausible system.

    [units] refers to the units of the RHS. Commas may be used to
    describe a vector or matrix, but may not. These data structures
    are expected to contain the same units throughout. This is primarily
    done to catch and avoid physically nonsensical operations.

    d# and #_dot will be used interchangeably to denote the time
    derivative of some quantity #.

    Running this script:
        You have three options. You can 1. run this script by executing
        it directly, 2. calling python on the script, or 3. using your
        preferred GUI IDE.

        1. Set the script to executable. You only need to run this once.
            $ chmod +x work.py 
          Then, simply execute the script 
            $ ./work.py 
        2. Ensure python3 is in your PATH 
            $ python3 work.py 
        3. Refer to the guides and manuals provided by your prefered GUI
           IDE.

    Using this script as a module:
        - See __all__ for a list of functions that can be used outside
          this module.

    Dependencies:
        See requirements.txt
"""

__all__ = ["get_eom", "simulate"]

from typing import Tuple
from collections.abc import Callable

import time

import numpy as np
from scipy.integrate import RK45
from matplotlib import pyplot as plt
import sympy
from sympy import Symbol, Function, Matrix, solve, diff, pprint, latex
from sympy import simplify
from sympy.utilities.lambdify import lambdify
# Add leading s to differentiate it as symbolic and avoid namespace
# conflict
from sympy import sin as ssin
from sympy import cos as scos

# PEP 008 says you can ignore typical naming conventions to match
# notation used in math formulas.
# pylint: disable=invalid-name


# May require LaTeX installed on your machine, just comment this out
# if it gives you grief.
plt.rcParams['text.usetex'] = True
# Interactive mode
plt.ion()

# Setup for pretty printing
sympy.init_printing(use_unicode=True)

# --- System parameters ---

# The length of the pendulum
length = 50e-3 # [m]
# The length of the pendulum to the center of gravity
length_CG = length/2 # [m]
# The mass of the pendulum
mass = 25e-3 # [kg]

# The radius of the motor
radius = 12.5e-3 # [m]

# Mass moment of inertia of thin rod about its end. Assume our pendulum
# can be approximated by a rod.
inertia = 1/3 * mass * length**2 # [kg m^2]
grav = 9.80665 # [m/s^2]

def get_eom() -> Callable[[float,float,float,float,float,float], float]:
    """
    Produces the equation of motion analytically, and applies defined
    numerical values for system parameters to produce an evaluatable
    expression.

    Returns:
        (Callable[[float, float, float, float, float, float], float]): a
        function that takes in
        (phi, phi_dot, theta, theta_dot, theta_ddot) respectively and
        returns phi_ddot.

        Units on input are [rad, rad/s, rad, rad/s, rad/s^2], and
        units of output are [rad/s^2].
    """
    m = Symbol("m", real=True, positive=True) # [kg]
    I = Symbol("I", real=True, positive=True) # [kg m^2]
    l = Symbol("l", real=True, positive=True) # [m]
    l_CG = Symbol("l_CG", real=True, positive=True) # [m]
    r = Symbol("r", real=True, positive=True) # [m]
    g = Symbol("g", real=True, positive=True) # [m/s^2]

    t = Symbol("t", real=True, positive=True) # [s]

    # Chosen (can be picked)
    theta = Function("theta", real=True)(t) # [rad]

    theta_dot = diff(theta, t) # [rad/s]
    theta_ddot = diff(theta_dot, t) # [rad/s^2]

    # Free (determined by system dynamics)
    phi = Function("phi ", real=True)(t) # [rad]

    phi_dot = diff(phi, t) # [rad/s]

    # Directions
    i_hat = Matrix([1, 0])
    j_hat = Matrix([0, 1])

    # position of contact with respect to inertial frame
    pos_C = r * Matrix([ssin(theta), -scos(theta)]) # [m]
    # position of CG of pendulum with respect to inertial frame
    # Assume pendulum is symmetrical
    pos_G = pos_C + l_CG * Matrix([ssin(phi), scos(phi)]) # [m]

    vel_G = diff(pos_G, t) # [m/s]

    # Potential Energy
    T = m*g * pos_G.dot(j_hat) # [J]
    # Kinetic Energy
    V = 1/2 * m * vel_G.dot(vel_G) + 1/2 * I * phi_dot**2 # [J]

    # Lagrangian
    L = T - V # [J]

    # Right hand side is 0
    eom_lhs = diff(L,phi_dot, t) - diff(L, phi) # [J/rad]
    eom_lhs = simplify(eom_lhs) # [J/rad]

    pprint(eom_lhs)
    print(latex(eom_lhs))

    # Provide numerical values for constants
    eom = eom_lhs.subs([
        (l, length),
        (l_CG, length_CG),
        (m, mass),
        (r, radius),
        (I, inertia),
        (g, grav),
    ])

    # Use the first (and only) solution
    ddphi_eq = solve(eom, diff(phi_dot))[0] # [rad/s^2]

    # Evaluates ddphi_eq to a float systematically.
    func = lambdify(
        [phi, phi_dot, theta, theta_dot, theta_ddot],
        ddphi_eq,
        "numpy",
    )

    return func


def get_thetas(
    t: float,
    y: np.ndarray[float]
) -> Tuple[float, float, float]:
    """
    Provides (dtheta, ddtheta) value for a given (observed)
    system state.
    """
    phi = y[0]
    dphi = y[1]

    k_p = 0.1
    k_e = 0.5

    return (dphi * k_e, 0)
    # f_0 = 1/(2*np.pi) * np.sqrt(grav/length)
    # omega = 40 * f_0
    # # omega = 10000
    # amp = np.pi / 32
    # amp = ((0.450 + 1.799/omega**2) + (np.sqrt(2)/omega))/2
    # # omega = 10
    # theta = amp * np.sin(omega * t)
    # dtheta = omega * amp * np.cos(100 * t)
    # ddtheta = omega**2 * amp * -np.sin(100 * t)
    # return (theta, dtheta, ddtheta)

def simulate() -> None:
    """Simulates the EOM"""
    ddphi: Callable = get_eom()

    # Let y = [phi, dphi, theta]

    # Initial conditions
    y_0 = np.array([0.1, 0, 0]) # [rad, rad/s, rad]

    t_0 = 0 # [s]
    # No rest for the wicked.
    t_f = float("inf") # [s]

    # Bad practice
    t_old = t_0 # [s]

    def get_pos(
        y_state: np.ndarray[float]
    ) -> Tuple[np.ndarray[float], np.ndarray[float]]:
        """Returns (x, y) as two arrays of coordinates to graph."""
        phi = y_state[0]
        theta = y_state[2]

        orig_x = 0 # [m]
        orig_y = 0 # [m]
        C_x = radius * np.sin(theta) # [m]
        C_y = radius * -np.cos(theta) # [m]

        # position of CG of pendulum with respect to inertial frame
        # Assume pendulum is symmetrical
        G_x = C_x + length_CG * np.sin(phi) # [m]
        G_y = C_y + length_CG * np.cos(phi) # [m]

        x_coord = np.array([orig_x, C_x, G_x]) # [m]
        y_coord = np.array([orig_y, C_y, G_y]) # [m]

        return (x_coord, y_coord)

    fig, ax = plt.subplots(figsize=(8,8))

    phi_0 = y_0[0]
    line = ax.plot(*get_pos(y_0), "-")[0]
    plt.xlim(-radius-length, radius+length)
    plt.ylim(-radius-length, radius+length)
    ax.set_aspect("equal", adjustable="box")

    plt.title("Plot")


    def redraw(t: float, y_state: np.ndarray[float]) -> None:
        """Redraws the canvas to plot the state of the pendulum."""
        x, y = get_pos(y_state)
        line.set_xdata(x)
        line.set_ydata(y)
        fig.canvas.draw()
        fig.canvas.flush_events()
        # time.sleep(t - t_old)

    def dy(t: float, y: np.ndarray[float]) -> np.ndarray[float]:
        """Does what you should expect (returns dy/dt given t, y)"""
        phi = y[0] # [rad]
        dphi = y[1] # [rad/s]

        theta = y[2] # [rad]
        dtheta, ddtheta = get_thetas(t, y) # [rad, rad/s, rad/s^2]

        print(t)

        return np.array([
            dphi,
            ddphi(phi, dphi, theta, dtheta, ddtheta),
            dtheta,
        ])


    ode = RK45(dy, t_0, y_0, t_f, max_step=1e-3)

    while ode.t <= ode.t_bound:
        status = ode.step()
        if status is not None:
            raise RuntimeError(f"Integration failed. {status}")
        redraw(ode.t, ode.y)
    # plt.plot(sln.t, sln.y[0,:])

    # plt.show()

def main() -> None:
    """Simulates the dynamic system"""
    simulate()

if __name__ == "__main__":
    main()
