#from src.inversion_algorithm import InversionAlgorithm, check_E_layer_existance, plot_inversion_algorithm
from src.handle_e_layer import handle_e_layer
from src.handle_f_layer import handle_f_layer
from src.plot_results import plot_results
from src.sao_reader import Ne_prof, Ne_map
from src.utils import check_E_layer_existance, print_centered, print_yellow
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt

currect_directory = os.getcwd()
sh_file = os.path.join(currect_directory, 'scripts', 'find_sao_files.sh')
subprocess.run(['bash', sh_file], check=True)

with open('sao_files_list.txt', 'r') as file_list:
    filenames = file_list.readlines()

# Process each .SAO file 
for name in filenames:
    filename = f'sao_files/{name.strip()}'  # Construct the full path
    filename_fig = f'{name.strip()}'
    map_ionograms, dates = Ne_prof(filename)
    #Ne_map(filename)
    for i, date in enumerate(dates):
        frq, vh, foE = map_ionograms[date]
        frq, vh = np.array(frq), np.array(vh)

        if foE is None: continue
        flag_E_layer = check_E_layer_existance(frq,vh,foE)
        if flag_E_layer:
            print()  
            print()   
            print_centered('SOLVING FOR IONOGRAM '+date+' FROM BATCH '+name[:-4])
            print_yellow('SOLVING E LAYER')
            QP = handle_e_layer(frq=frq, vh=vh, foE=foE)
            print_yellow('SOLVING F LAYER')
            QP = handle_f_layer(QP=QP.copy(), frq=frq, vh=vh)  
            plot_results(frq, vh, QP['plasma_frequency'], QP['real_height'], batch_name=filename_fig[:-4], date=date, i=i)
  