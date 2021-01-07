#!/usr/bin/env python
# coding: utf-8
# Author: Bora Isildak
# Date  : 16/12/2020

#To estimate the cross section for the 2 → 3 bremsstrahlung process pp → pA0X, we use
#the Fermi-Weizsacker-Williams method of virtual quanta 
#https://arxiv.org/pdf/1311.3870.pdf

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import ROOT
import pandas as pd

def h(z, pt, m_a):
    m_p = 0.938 # GeV
    return pt**2 + (1-z) * m_a**2 + z**2 + m_p**2

def w_z(theta, pt, m_a, epsilon):
    m_p      = 0.938 # GeV
    alpha    = 1./1.
    p_proton = 13000 # GeV momentum of the incoming proton
    
    p = pt / np.sin(theta)
    z = p * np.cos(theta) / p_proton
    
    a = epsilon**2 * alpha / (2 * np.pi * h(z, pt, m_a))
    b = (1 + (1 - z)**2) / z
    c = 2 * z * (1 - z) 
    d = ((2 * m_p**2 + m_a**2) / h(z, pt, m_a)) - (z**2 * 2 * m_p**4 / h(z, pt, m_a)**2)
    e = 2 * z * (1 - z) * (z + (1 - z)**2) * m_p**2 * m_a**2 / h(z, pt, m_a)**2
    f = 2 * z * (1 - z)**2 * m_a**4 / h(z, pt, m_a)**2

    return a * (b - c * d + e + f)

def q_min(theta, pt, m_a):
    m_p      = 0.938 # GeV
    E_p      = 13000 # GeV is the incident proton energy in the rest frame of the other proton.
    p_proton = 13000 # GeV momentum of the incoming proton
    
    p = pt / np.sin(theta)
    z = ((pt**2 + p**2 * np.cos(theta)**2)/ p_proton**2)**0.5
    return (1 / (4 * z * E_p**2 * (1 - z)**2)) * (pt**2 + (1 - z) * m_a**2 + z**2 * m_p**2)

def sqrt_s_prime(theta, pt, m_a):
    m_p = 0.938 # GeV
    E_p = 6500 # GeV is the incident proton energy.
    
    p = pt / np.sin(theta)
    E_a = (m_a**2 + p**2)**0.5
    
    return (2 * m_p * (E_p - E_a))**0.5

def get_decay_prob(theta, pt, m_a, epsilon, br=1.0):
    p = pt / np.sin(theta)
    E_a = (m_a**2 + p**2)**0.5
    decay_length = 80 * br * (1e-5/epsilon)**2 * (E_a/1e+3) * (100/(m_a*1e+3))**2 
    # https://arxiv.org/pdf/1708.09389.pdf - page 6
    L1 = 400
    L2 = 410
    prob = ROOT.TMath.Exp(-L1/decay_length)-ROOT.TMath.Exp(-L2/decay_length)
    #prob = (L2-L1)/decay_length

    return prob

def form_factor(m_a):
    # calculate the form factor
    m_rho     = 782
    width_rho = 146
    fv_list   = [0.616, 0.223, -0.339]
    
    form_factor = 0 + 0j
    for i in range(len(fv_list)):
        form_factor += (fv_list[i]*m_rho**2)/(m_rho**2-m_a**2-(m_rho*width_rho*1j))
    
    m_omega     = 770
    width_omega = 8.5
    
    fv_list   = [1.011, -0.881, 0.369]
    
    for i in range(len(fv_list)):
        form_factor += (fv_list[i]*m_omega**2)/(m_omega**2-m_a**2-(m_omega*width_omega*1j))
    
    return np.real(form_factor * np.conj(form_factor))


# get the elastic pp cross section from the data and interpolate
pp_el_xsec_data = np.genfromtxt('pp_inelastic_xsec.csv', delimiter=',').T
pp_el_xsec      = interp1d(pp_el_xsec_data[0], pp_el_xsec_data[1], kind='cubic')


def calculate_n_expected(m_a, epsilon, lumi, theta_max):

    data        = []
    dpt         = 0.1
    dtheta      = 1.0e-4
    pt_range    = np.arange(0.01, 20.0, dpt)
    theta_range = np.arange(0.1e-3, 10e-3, dtheta)
    m_p         = 0.938 # GeV

    form_factor_2 = form_factor(m_a)

    for pt in pt_range:
    
        for theta in theta_range:
        
            p      = pt / np.sin(theta)
            dp     = dpt / np.sin(theta)
            E_a    = (m_a**2 + p**2)**0.5
            dz     = dp * np.cos(theta) / m_p
        
            dsigma = form_factor_2 * pp_el_xsec(sqrt_s_prime(theta, pt, m_a))                     * w_z(theta, pt, m_a, epsilon) * get_decay_prob(theta, pt, m_a, epsilon, 1)                     * 2 * pt * dpt * dz
                
            condition_1 = (q_min(theta, pt, m_a)**2 < 250**2)
            condition_2 = (not np.isnan(dsigma))
            condition_3 = (E_a > 100 and E_a < 13000)
        
            if  condition_1 and condition_2 and condition_3 :
                data.append(np.array([theta, theta+dtheta, p, p+dp, dsigma]))

    data = np.asarray(data)
    mask = data[:, 1] < theta_max

    return lumi*np.sum(data[:,4][mask])


print("%.2f" %(calculate_n_expected(1.0e-1, 1.0e-5, 3.0e+9, 0.5e-3)))

#df = calculate_n_expected(3.0e-1, 3.0e-6, 3.0e+9)
#np.savetxt("dark_photon_bremsstrahlung_data.csv", data, delimiter=",",fmt="%.4e", header= "\
#pp-> p+A+X bremsstrahlung cross sections by FWW approximation\n \
#p and theta values given below belong to A (dark photon)\n \
#theta_min,theta_max,p_min[GeV],p_max[GeV],xsec[pb]")

#df = pd.read_csv("dark_photon_bremsstrahlung_data.csv", skiprows=2)