import matplotlib.pyplot as plt
from src.get_ionogram import get_ionogram

params = {'mew': 0, 'mec': 'none', 'alpha': 0.4}

def plot_results(ori_frq, ori_vh, fp, rh, filename):
    marker = 6
    plt.figure()
    s_f, s_vh = get_ionogram(rh, fp, fp)
    plt.plot(fp, rh, linewidth = 3.5, label='reconstructed fp profile', color='orange',alpha=0.9)
    plt.plot(s_f, s_vh, 's', markersize = marker, label='reconstructed ionogram', mew=0.5, mec='black', mfc = 'none')
    plt.plot(ori_frq, ori_vh, 'o', color='magenta', markersize = marker, label='original ionogram', **params)
    plt.legend()
    plt.savefig(filename+'.pdf')
    plt.close()