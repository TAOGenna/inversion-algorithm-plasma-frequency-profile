import matplotlib.pyplot as plt
from src.get_ionogram import get_ionogram

params = {'mew': 0, 'mec': 'none', 'alpha': 0.7}

def plot_results(ori_frq, ori_vh, fp, rh, filename):
    plt.figure()
    plt.plot(fp, rh, '^', markersize = 3, label='reconstructed fp profile', **params)
    s_f, s_vh = get_ionogram(rh, fp, fp)
    plt.plot(s_f, s_vh, 's', markersize = 2.5, label='reconstructed ionogram', mew=0.5, mec='black', mfc = 'none')
    plt.plot(ori_frq, ori_vh, 'o', color='magenta', markersize = 2, label='original ionogram', **params)
    plt.legend()
    plt.savefig(filename+'.pdf')
    plt.close()