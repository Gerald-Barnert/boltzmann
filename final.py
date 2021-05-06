# -*- coding: utf-8 -*-
import camb
import classy
from classy import Class
from camb import model
import numpy as np
import matplotlib.pyplot as plt



def P_k(code, redshifts, H0, omega_b, omega_cdm, kmax):
    if code == "camb":
        