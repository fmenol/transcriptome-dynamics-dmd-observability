import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from copy import deepcopy

from pydmd import DMD
from pydmd.bopdmd import BOPDMD
from pydmd.plotter import plot_eigs, plot_summary

_  = pd.read_csv('cdii.csv').drop(columns=['Unnamed: 0'])
Xn = _.values

t=np.linspace(0,230,24)

import pickle

with open('df_result_columns.pickle', 'rb') as file:
    gene_names = pickle.load(file)
    
bopdmd = BOPDMD(
svd_rank=24,
num_trials=0,
trial_size=0.8,
varpro_opts_dict={"verbose":False, "tol":0.8},
eig_constraints={"imag"},
use_proj=False,
mode_regularizer="l0",
compute_A=True
)


bopdmd.fit(Xn, t)