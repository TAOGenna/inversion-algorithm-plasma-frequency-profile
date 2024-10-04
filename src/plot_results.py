import matplotlib.pyplot as plt
from src.get_ionogram import get_ionogram
import os

params = {'mew': 0, 'mec': 'none', 'alpha': 0.4}

def plot_results(ori_frq, ori_vh, fp, rh, batch_name, date, i):
    marker = 6
    plt.figure()
    s_f, s_vh = get_ionogram(rh, fp, fp)
    plt.plot(fp, rh, linewidth = 3.5, label='reconstructed fp profile', color='orange',alpha=0.9)
    plt.plot(s_f, s_vh, 's', markersize = marker, label='reconstructed ionogram', mew=0.5, mec='black', mfc = 'none')
    plt.plot(ori_frq, ori_vh, 'o', color='magenta', markersize = marker, label='original ionogram', **params)
    plt.legend()
    plt.title(date)
    currect_directory = os.getcwd()
    os.makedirs(currect_directory+'/results/'+batch_name, exist_ok=True)
    plt.savefig(currect_directory + '/results/' + batch_name + '/' + str(i) + '.pdf')
    plt.close()