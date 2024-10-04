from src.inversion_algorithm import InversionAlgorithm, check_E_layer_existance, plot_inversion_algorithm
from src.handle_e_layer import handle_e_layer
from src.handle_f_layer import handle_f_layer
from src.plot_results import plot_results
from src.sao_reader import Ne_prof, Ne_map
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt

# PLAN: 
# pass through the .sao file through the Ne_prof function to fetch all the scaled data 
# pass each of the ionograms through the InversionAlgorithm
# output a single image: the autoscaled ionogram | the reconstructed fp profile and ionogram

currect_directory = os.getcwd()
sh_file = os.path.join(currect_directory, 'scripts', 'find_sao_files.sh')
subprocess.run(['bash', sh_file], check=True)

with open('sao_files_list.txt', 'r') as file_list:
    filenames = file_list.readlines()

# Process each .SAO file
for name in filenames:
    filename = f'sao_files/{name.strip()}'  # Construct the full path
    map_ionograms, dates = Ne_prof(filename)
    #Ne_map(filename)
    for i, date in enumerate(dates):
        frq, vh, foE = map_ionograms[date]
        frq, vh = np.array(frq), np.array(vh)
        
        if foE is None: continue
        # plt.figure()
        # plt.plot(frq,vh,'o', label='original ionogram')
        # plt.savefig('original_ionogram.png')
        # plt.close()
        flag_E_layer = check_E_layer_existance(frq,vh,foE)

        if flag_E_layer:
            aE, bE, rmE, fp_E, rh_E, index = handle_e_layer(frq=frq, vh=vh, foE=foE)
            QP = {'plasma_frequency': fp_E, 'real_height': rh_E, 
                  'numt': index, 'a_0': aE, 'r_m0': rmE, 'b_0': bE,
                  'f_c0': np.sqrt(aE)}
            plot_results(frq, vh, QP['plasma_frequency'], QP['real_height'], filename='cmp E layer')
            #print('------------------------------------------------------------------------')
            #print(QP)
            #print('------------------------------------------------------------------------')
            print('FINISHED E LAYER, STARTING F LAYER')
            QP = handle_f_layer(QP=QP.copy(), frq=frq, vh=vh)
            
            plot_results(frq, vh, QP['plasma_frequency'], QP['real_height'], filename='final comparison')
            
            break # we work with just one example 
            # profile = InversionAlgorithm(data_r = vh , data_f = frq, foE = foE)
            # path = os.getcwd() + '/' + filename[:-4]
            # plot_inversion_algorithm(frq, vh, profile.plasma_frequency, profile.real_height, profile.frequency, profile.virtual_height, path, date, foE, i)
            # break
