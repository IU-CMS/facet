import ROOT
import numpy as np
import pandas as pd
from itertools import product
import matplotlib.pyplot as plt
from tqdm import tqdm

from scipy.interpolate import interp1d

ROOT.ROOT.EnableImplicitMT(8)

def decay_mother(mother, daugther1_mass, daugther2_mass):
    
    event  = ROOT.TGenPhaseSpace()
    masses = np.array([daugther1_mass, daugther2_mass], dtype='float64')
    
    mother = ROOT.TLorentzVector(mother.X(), mother.Y(), mother.Z(), mother.T())
    event.SetDecay(mother, 2, masses)
    
    weight     = event.Generate()
    daugther1  = event.GetDecay(0)
    daugther2  = event.GetDecay(1)
    
    return mother, daugther1, daugther2

def get_prob(tree_reader, i_event, epsilon, ad_mass, daugther2_mass, br):
    
    treeReader.SetEntry(i_event)
    
    mother, ad, gamma = decay_mother(pi0_p4, ad_mass, 0)

    decay_length = 80 * br * (1e-5/epsilon)**2 * (ad.E()/1e+3) * (100/(ad_mass*1e+3))**2 # https://arxiv.org/pdf/1708.09389.pdf - Eq.7

    L1 = 96
    L2 = 116
    
    prob = ROOT.TMath.Exp(-L1/decay_length)-ROOT.TMath.Exp(-L2/decay_length)
    
    return mother, ad, prob

# aD->e+e- branching ratios from https://arxiv.org/pdf/1505.07459.pdf Fig.2
ad2ee_pair_br = np.genfromtxt('ad2ee_pair_br.csv', delimiter=',').T
br_e          = interp1d(ad2ee_pair_br[0], ad2ee_pair_br[1], kind='linear')

ROOT.gSystem.Load("libPhysics.so")

det_radius_in  = 0.125
det_radius_out = 0.500

ROOT.gRandom.SetSeed(123)

aD_p4 = ROOT.TLorentzVector(0., 0., 0., 0.)

file_name = "data_14TeVpp_eta.root"

file = ROOT.TFile.Open(file_name)

treeReader = ROOT.TTreeReader("tree", file)

outfile_name = 'pion2ad_test.csv'
branch_name = 'pi0_p4'
mother_mass = 139.0e-3 # in GeV
br_mother2gammagamma = 0.99 

if 'eta' in file_name:
    branch_name  = 'eta_p4'
    outfile_name = 'eta2ad_test.csv'
    mother_mass = 547.862e-3 # in GeV
    br_mother2gammagamma = 0.39
    
pi0_p4 = ROOT.TTreeReaderValue(ROOT.TLorentzVector)(treeReader, branch_name)
    
epsilon_exponents = np.arange(-7,-3, 0.05)
epsilon_range = [10**i for i in epsilon_exponents]
ad_mass_range = np.linspace(10.0e-3, mother_mass, 50)
print(epsilon_range, ad_mass_range)

lumi = 300.0e+16 # 3ab^-1
total_xsec = 110.69e-3 # 110.69 mb
n_generated = 1e+6/100

df = pd.DataFrame(columns=['epsilon','mass', 'n_expected'])

df = pd.DataFrame(columns=['epsilon','mass', 'n_expected'])

for ad_mass in ad_mass_range:
    
    for epsilon in epsilon_range:
        
        sum_prob = 0.0
        br_adgamma = 2 * epsilon**2 * (1-(ad_mass/mother_mass)**2)**3 * br_mother2gammagamma # https://arxiv.org/pdf/1708.09389.pdf - Eq.9
        print("running for eps=%.2E  ad_mass=%.2E" %(epsilon, ad_mass))
        
        for i_event in range(int(n_generated)):
            mother, ad, prob = get_prob(treeReader, i_event, epsilon, ad_mass, 0.000e-3, br_e(ad_mass))
            
            if ad.Theta()<np.arctan(0.5/96) and ad.Theta()>np.arctan(0.125/96) and ad.E()>100:
                sum_prob += (lumi) * (total_xsec)* (1/n_generated) * br_adgamma * prob
                
        print("n_expected", sum_prob)
        df = df.append(pd.Series([epsilon, ad_mass, sum_prob], index=df.columns), ignore_index=True)

df.to_csv(outfile_name)