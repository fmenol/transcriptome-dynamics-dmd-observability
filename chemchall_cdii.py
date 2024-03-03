import numpy as np
import pandas as pd
import pickle
import seaborn as sns
import scipy as sp
from scipy.signal import savgol_filter as savgol 
from copy import deepcopy    
from random import choices
from scipy.stats import ttest_ind
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colors as mcolors
from matplotlib import cm
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib
npg_colors = ['#00A087FF','#E64B35FF','#3C5488FF','#4DBBD5FF','#F39B7FFF','#8491B4FF','#91D1C2FF','#DC0000FF','#7E6148FF','#B09C85FF']
npg_color_list = []
for this_color in npg_colors: 
    npg_color_list.append(list(int(this_color[1:][i:i+2], 16)/255 for i in (0, 2, 4)))
# plt.style.use('dark_background')
# dark_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  
plt.rcParams.update({'font.size':25});
plt.rcParams.update({'axes.linewidth':1.5})
plt.rc('lines',linewidth=3);
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False 
ncolors = 10
plt.rc('axes', axisbelow=True)
figsize=(9.32*3*0.6,3*3.74*0.6)

from preprocess import *
from dmd import *
from sensor_placement import *
# from query_uniprot import *
from regions import *
from plot_ellipsoid import *


#Load the matrix with all the data
data_full = pd.read_csv("./data_chem/features.csv")
#Load induction time profiles
u_ts = pd.read_csv("./data_chem/targets.csv")
#Load genes
genes = pd.read_csv("./data_chem/loading_record.csv")


#List and print all experiments
exp_id = data_full.exp_id.unique()
# print(exp_id)

#Select the experiment
exp_of_interest = 1 #this selects the first experiment only

#Subselecting data for the experiment of interest
df_exp = data_full[data_full.exp_id == exp_id[exp_of_interest-1]].drop(columns=["Unnamed: 1", "exp_id"])

#Renaming columns
df_exp.columns = genes.gene_name

#Filtering out spots with no strain averaging duplicate genes
df_exp=df_exp.iloc[:,df_exp.columns.notna()]



# This section extracts the input profile
# Display the unique values in the "induction_state" column
u_ts = deepcopy(u_ts[u_ts.exp_id == exp_id[exp_of_interest-1]])

unique_induction_states = u_ts["induction_state"].unique()

# Define a mapping for encoding of u_ts
encoding_map = {
    "('MilliQ_ddH2O',)": 0,
}

# Assign subsequent integers to other unique values
for idx, state in enumerate(unique_induction_states, start=0):
    if state not in encoding_map:
        encoding_map[state] = idx

u_ts["induction_state_encoded"] = u_ts["induction_state"].map(encoding_map)


# Create a dictionary to store start and end indices for each induction state occurrence
induction_intervals = {}

# Initialize variables to keep track of the current state and its start index
current_state = None
start_index = None

# Loop through rows in the DataFrame
for idx, row in u_ts.iterrows():
    state = row["induction_state"]
    
    # If we encounter a new state or reach the end of the DataFrame
    if current_state != state or idx == u_ts.index[-1]:
        # If we were tracking a state, save its interval (except for the first iteration)
        if current_state is not None:
            # Adjust end index for the last entry in the DataFrame
            end_index = idx if idx == u_ts.index[-1] and current_state == state else idx - 1
            if current_state in induction_intervals:
                induction_intervals[current_state].append((start_index, end_index))
            else:
                induction_intervals[current_state] = [(start_index, end_index)]
        
        # Start tracking the new state
        current_state = state
        start_index = idx

induction_intervals


#In the original dataset there are several "biological replicates", i.e. spots of the chip where the same gene is tracked
#This snippet of code averages all these spots to obtain a final dataset with one column per gene
duplicated_columns = df_exp.columns[df_exp.columns.duplicated(keep=False)].unique()

for i in duplicated_columns:
    print("The gene ",i, " has ", df_exp[i].shape[1], " spots associated to it (biological replicates)")

# Create a new DataFrame to store the results
df_result = pd.DataFrame()

# Process non-duplicated columns
for column in df_exp.columns:
    if column not in duplicated_columns:
        df_result[column] = df_exp[column]

# Process duplicated columns
for dup_col in duplicated_columns:
    # Take the average of duplicated columns and assign to the result DataFrame
    df_result[dup_col] = df_exp[dup_col].mean(axis=1)

df_result.head()


df = deepcopy(df_result) #renaming the array for simplicity (and because I'm lazy :P)

print("There are "+str(len(induction_intervals))+" in this experiment.")
for i,x in enumerate(induction_intervals):
    print(i+1,". ",x)
# stimulus_of_interest=int(input("Which input are you interested in?"))-1
stimulus_of_interest = 1


start_time = [i[0] for i in induction_intervals[u_ts.induction_state.unique()[stimulus_of_interest]]]
end_time = [i[1] for i in induction_intervals[u_ts.induction_state.unique()[stimulus_of_interest]]]
replicates = len(induction_intervals[u_ts.induction_state.unique()[stimulus_of_interest]])

input_label = u_ts.induction_state.unique()[stimulus_of_interest][2:-3]

print("\n You selected: "+u_ts.induction_state.unique()[stimulus_of_interest])


#This the maximum length of the time-series (i.e. the number of rows) in the 3D matrix, given that  
max_w = np.min([(end_time[i] - start_time[i]) for i in range(len(start_time))])

#Creating the 3D matrix
data = np.empty(shape=(max_w+1, len(df.columns), replicates))


#Assembling the 3D array, based on the stimulus selected (note that the 3rd dimension might be a singleton)
for  i in range(replicates):
    data[:,:,i] = df.iloc[start_time[i]:start_time[i]+max_w+1,:].values


from scipy.signal import savgol_filter as savgol 

data = data.transpose(1,0,2)
data_f = savgol(data, window_length=5, polyorder=2, axis=1)


from preprocess import standardize_time_series
data_fs = standardize_time_series(data_f, data_f.shape[1], data_f.shape[2])

data = deepcopy(data_fs)

Xp = data[:,:-1].reshape(len(data),(data.shape[1]-1)*data.shape[2],order='F') 
Xf = data[:,1:].reshape(len(data),(data.shape[1]-1)*data.shape[2],order='F')

r = 23
U,s,Vh = np.linalg.svd(Xp)
U_r = U[:,0:r] 
s_r = s[0:r]
Vh_r = Vh[0:r,:]

Atilde = U_r.T @ Xf @ Vh_r.T @ np.diag(1/s_r)
A = Xf @ np.linalg.pinv(Xp)

data_red = np.zeros((r,data.shape[1],data.shape[2]))
data_red[:,:,0] = np.dot(U_r.T ,data[:,:,0])
data_red[:,:,1] = np.dot(U_r.T ,data[:,:,1])

## This snippet of code uses the projections in the space of the modes to calculate the R^2 of the reconstruction
ntimepts = data_red.shape[1]
nreps = data_red.shape[2]

# Create the array with the data, reshaped to flatten the repeats
X = data_red[:, :ntimepts].reshape(len(data_red), (ntimepts)*nreps, order='F')
X_pred = np.zeros((Atilde.shape[0], ntimepts*nreps))

# Make the predictions taking the initial conditions and projecting from there
count = 0
for i in range(0, nreps):
    x_test_ic = X[:,i*(ntimepts):i*(ntimepts)+1]
    for j in range(0, ntimepts):
        X_pred[:, count:count+1] = np.dot(np.linalg.matrix_power(Atilde, j), x_test_ic)
        count += 1

feature_means = np.mean(X, axis=1).reshape(len(X), 1)

cd = 1 - ((np.linalg.norm(X-X_pred, ord=2))**2)/(np.linalg.norm(X-feature_means, ord=2)**2)

X_pred_red = X_pred.reshape(len(data_red), data_red.shape[1], data_red.shape[2])

L,W = np.linalg.eig(Atilde)

