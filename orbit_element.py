#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: orbit_element.py
"""
Created on Tue Jul 20 15:01:48 2021

@author: Neo(niu.liu@nju.edu.cn)
"""


import numpy as np


def apriori_init(**kwargs):
    """Initialize a priori value of orbit elements

    """

    # List of parameters to be estimated
    est_pmt_list = ["a", "e", "i", "omega", "Omega", "P", "T0"]

    # To see if the a priori value of a is assigned
    if "a0" in kwargs.keys():
        a0 = kwargs["a0"]
        est_pmt_list.remove("a")
    else:
        a0 = 1

    # To see if the a priori value of e is assigned
    if "e0" in kwargs.keys():
        e0 = kwargs["e0"]
        est_pmt_list.remove("e")
    else:
        e0 = 0

    # To see if the a priori value of i is assigned
    if "i0" in kwargs.keys():
        i0 = kwargs["i0"]
        est_pmt_list.remove("i")
    else:
        i0 = 0

    # To see if the a priori value of omega is assigned
    if "omega0" in kwargs.keys():
        omega0 = kwargs["omega0"]
        est_pmt_list.remove("omega")
    else:
        omega0 = 0

    # To see if the a priori value of Omega is assigned
    if "Omega0" in kwargs.keys():
        Omega0 = kwargs["Omega0"]
        est_pmt_list.remove("Omega")
    else:
        Omega0 = 0

    # To see if the a priori value of P is assigned
    if "P0" in kwargs.keys():
        P0 = kwargs["P0"]
        est_pmt_list.remove("P")
    else:
        P0 = 1

    # To see if the a priori value of T0 is assigned
    if "T00" in kwargs.keys():
        T00 = kwargs["T00"]
        est_pmt_list.remove("T0")
    else:
        T00 = 0

    return a0, e0, i0, omega0, Omega0, P0, T00, est_pmt_list


vim orbele_2_tancoo(orb_ele):
    """Calculate tangential coordinates from orbital elements.

    Parameters
    ----------
    Seven orbital elements orb_ele=(a, e, i, omega, Omega, P, T0)^T

    ########## TBD  ##########
    ## 1) Add a input parameters of time interval
    ## 2) orb_ele could be a list or ndarray, dict
    ########## TBD  ##########

    Returns
    -------
    tangential coordinates (x, y)
    """

    # just assume that orb_ele is a list
    a, e, i, omega, Omega, P, T0 = orb_ele

    # I need a function to calculate the Thiele-Innes constants
    # Anyway, assume now I had one
    # Calculate the Thiele-Innes constants
    A, B, F, G = calc_thiele_innes_con(a, i, omega, Omega)

    # Calcukate the mean anomaly and then eccentric anomaly
    T = np.linspace(0, P, 100)
    M = 2 * np.pi * (T-T0) / P
    E = some_func(M)

    fac1 = np.cos(E) - e
    fac2 = np.sqrt(1-e*e) * np.sin(E)
    x = A * fac1 + F * fac2
    y = B * fac1 + G * fac2

    return x, y
