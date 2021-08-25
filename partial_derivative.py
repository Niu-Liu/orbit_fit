#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: partial_derivative.py
"""
Created on Fri Jul 30 16:17:18 2021

@author: Neo(niu.liu@nju.edu.cn)

Calculate the partial derivative of tangental coordinates (x, y) wrt.
orbital elements (a, e, i, Omega, omega, T0, P), where

a     : apparent semimajor axis of the photocenter 
e     : eccentricity
i     : orbit inclination in radian (zero for face-on)
omega : pariastron longitude in radian
Omega : position angle of node in radian
T0    : periastron time
P     : orbital period in year

The main references are

Goldin &Makarov, APJSS, 2006 (GM2006)

I use tool from PyAstronomy to solve the Kepler's equation.
(https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/keplerOrbit.html)

"""

import numpy as np
from numpy import sin, cos, pi

from PyAstronomy import pyasl
# Instantiate the solver
ks = pyasl.MarkleyKESolver()


def calc_ThieleInnes_const(a, i, omega, Omega):
    """Calculate Thiele-Innes constants

    The definition could be found at Eq.(A5) in GM2006.

    Parameters
    ----------
    TBD

    Returns
    -------
    TBD
    """

    A = a * (cos(omega) * cos(Omega) - sin(omega) * sin(Omega) * cos(i))
    B = a * (cos(omega) * sin(Omega) + sin(omega) * cos(Omega) * cos(i))
    F = a * (-sin(omega) * cos(Omega) - cos(omega) * sin(Omega) * cos(i))
    G = a * (-sin(omega) * sin(Omega) + cos(omega) * cos(Omega) * cos(i))

    return A, B, F, G


def par_E_P(T, e, T0, E, P):

    return 2 * pi * (T-T0) / P*P / (e*cos(E) - 1)


def par_E_e(e, E):

    return sin(E) / (1 - e * cos(E))


def par_E_T0(e, E, P):

    return -2 * pi / P / (1 - e * cos(E))


def par_ABFG_i(a, i, omega, Omega):
    """

    """

    par_A_i = a * sin(omega) * sin(Omega) * sin(i)
    par_B_i = -a * sin(omega) * cos(Omega) * sin(i)
    par_F_i = a * cos(omega) * sin(Omega) * sin(i)
    par_G_i = -a * cos(omega) * cos(Omega) * sin(i)

    return par_A_i, par_B_i, par_F_i, par_G_i


def par_xy_wrt_a(x, y, a, e, i, omega, Omega, P, T0, T):
    """

    """

    # Compute the mean anomaly
    M = 2 * pi * (T-T0) / P

    # Solve the Kepler's equation to get the eccentric anomaly
    E = ks.getE(M, e)

    # Calculate the Thiele-Innes constants
    A, B, F, G = calc_ThieleInnes_const(a, i, omega, Omega)

    fac1 = sin(E)
    fac2 = sqrt(1-e*e) * cos(E)
    fac3 = cos(E) - e
    fac4 = 1 + fac1 * par_E_e
    fac5 = fac1 * fac3 / sqrt(1-e*e) / (1-e*cos(E))
    fac6 = sqrt(1-e*e) * sin(E)

    par_x_a = x / a
    par_y_a = y / a

    par_x_P = (-A * fac1 + F * fac2) * par_E_P(T, e, T0, E, P)
    par_y_P = (-B * fac1 + G * fac2) * par_E_P(T, e, T0, E, P)

    par_x_e = -A * fac4 + F * fac5
    par_y_e = -B * fac4 + G * fac5

    par_x_T0 = (-A * fac1 + F * fac2) * par_E_T0(e, E, P)
    par_y_T0 = (-B * fac1 + G * fac2) * par_E_T0(e, E, P)

    par_x_omega = fac3 * F + fac6 * (-A)
    par_y_omega = fac3 * G + fac6 * (-B)

    par_x_Omega = fac3 * (-B) + fac6 * (-G)
    par_y_Omega = fac3 * A + fac6 * F

    par_A_i, par_B_i, par_F_i, par_G_i = par_ABFG_i(a, i, omega, Omega)
    par_x_i = fac3 * par_A_i + fac6 * par_F_i
    par_y_i = fac3 * par_B_i + fac6 * par_G_i
