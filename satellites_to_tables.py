import numpy as np
import pandas as pd
import yt
import ytree
import numpy as np
import yt.utilities.physical_constants as constants
from yt.mods import *
from yt.units.yt_array import YTQuantity
from yt.units import *
from matplotlib import pyplot as plt
import heapq
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
from matplotlib import rcParams
from matplotlib.patches import Circle, PathPatch
import itertools
import matplotlib.colors as mcolors
import json
import warnings
warnings.filterwarnings('ignore')
path_halos = "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/Rockstar_GARROTXA/"
path_rockstar_data = "/home/bego/GARROTXA/satellites_rockstar_analysis/results/"
path_rockstar_tables = "/home/bego/GARROTXA/satellites_rockstar_analysis/satelites_tables/"
#path_rockstar = "/home/bego/GARROTXA/satellites_rockstar_analysis/"
path_csv = "/media/temp/bego/snapshots_resim/"
path_datos = "/home/bego/GARROTXA_copia/datos_GARROTXA_resim/"
path_results = "/home/bego/GARROTXA/GalaDyn/results/"
path_crossmatch = "/home/bego/GARROTXA/satelites_crossmatch/"
path_figures_acceleration = "/home/bego/GARROTXA/aceleration_figures/"
path_figures = "/home/bego/GARROTXA/acceleration_figures/"
path_acceleration = "/home/bego/GARROTXA/acceleration/"
path_disk = "/home/bego/GARROTXA/disco/"
seconds_to_Myr = 3.15576e+16


snapshots_analysis = [520,523,525, 527,530,532,535, 537,539,541,
   543, 545,547, 550, 553, 555,557, 
   560, 563, 565, 567,570,573, 575, 577, 580,
   583, 585,587,590, 592,594,
   596,598,
  600, 602, 604, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 
 629, 630, 632, 634, 636, 639, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 
 664, 666, 668,670, 672, 674, 676, 679, 681, 682, 684, 687, 689,
 690, 692, 694, 698, 704,  706, 708,711, 712,714, 716,
 718, 720, 722, 724, 726, 728, 731, 732, 734, 736, 739, 740, 742, 744, 746, 748, 751,752,
 755, 756, 758, 761,763, 764, 766, 768, 770, 772, 774, 776, 778, 780, 
 782, 784, 786, 788, 790, 792, 794, 797, 798, 
802, 805, 806, 808, 810, 812, 814, 816,
818, 820, 822, 824, 826, 828, 830, 832, 834, 836, 839, 840, 842, 844, 846, 848, 850,
853, 855, 856, 858, 860, 862, 864, 867, 870, 872, 875, 877, 879, 881, 883, 884, 888,
890, 892, 894, 898,
 900, 902, 904, 907, 908, 910, 912, 915, 916, 918, 921, 922, 924, 927, 929, 
 930, 932, 934, 937,
 939, 941,942, 944, 946, 948, 950, 952, 954,956, 
 958, 961, 963, 965, 966, 968, 970, 972, 974, 976, 979,
 980, 982, 984, 989, 990, 993, 994, 996, 999]
datos_edades =  pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)
lookback = []
for i,snapshot in enumerate(snapshots_analysis):
    lb = datos_edades.loc[datos_edades['Snapshot'] == snapshot, 'Lookback'].iloc[0]
    lookback.append(lb)


class Rockstar:
        def __init__(self, name):
            self.name = name
            self.data = None
            self.list_halos = None
            def read_data():
            #    print(f"Initializing {self.name}")
                self.data= pd.read_csv(path_rockstar_data + f"DF_halos_{self.name}.csv", index_col = 0)
                with open(path_rockstar_data + f"coordenadas_de_halos_{self.name}.json") as json_file:
                    self.list_halos=json.load(json_file)
                for key in self.data["ID"]:
                   # print(np.mean(data_805[f"{key}_R"]))
                    if np.mean(self.list_halos[f"{key}_R"]) < 2:
                        self.data= self.data[self.data.ID != key]

                    if np.isnan(self.list_halos[f"{key}_R"]).all() == True:
                        self.data = self.data[self.data.ID != key]


            read_data()

            
datos_crossmatch = pd.read_csv( path_rockstar_data + "crossmatch_of_halos.csv", sep = ",")
list_of_halos  = pd.read_csv( path_rockstar_data + "list_of_halos.csv", sep = ",", index_col = 0)
rocksana_snapshots = [620, 689, 740, 805, 850,910, 999]

for halo, corresp in zip(list_of_halos["ID"], list_of_halos["Snapshot"]):

    print(f"analysing {halo} in{corresp}")
    csv_dato_halo =  pd.DataFrame([], columns=['Snapshot',  "Lookback", "X","Y","Z","R", "Mass", "Mstars" ])
    ind = np.digitize(snapshots_analysis, bins = rocksana_snapshots) -1
    ind[ind < 0] = 0
    halo_ref = Rockstar (int(corresp))
    crossmatch = datos_crossmatch.loc[datos_crossmatch[f'{int(corresp)}'] == halo].iloc[0]

    for i,name in enumerate(snapshots_analysis):

        #print(name)
        halo_ind = rocksana_snapshots[ind[i]]
        halo_ind_loaded= Rockstar(int(halo_ind))

        #key_crossmatch = crossmatch.loc[crossmatch[f'{int(halo_ind)}']].iloc[0]
       # print(key_crossmatch)
        key_crossmatch = crossmatch[f"{int(halo_ind)}"]
        if np.isnan(key_crossmatch) == True:
            cont = 0
            #print(cont)
            while np.isnan(key_crossmatch) == True:
                cont = cont +1
             #   print(cont)
                try:
                    key_crossmatch = crossmatch[f"{int(rocksana_snapshots[ind[i+cont]])}"]

                except:
                    break

        #    print(f"No hay datos en estos snapshots!")

            try:    
                halo_ind_lolo = rocksana_snapshots[ind[i+cont]]
                halo_ind_loaded_lolo= Rockstar(int(halo_ind_lolo))
              #  print(f"calculating mass! with {key_crossmatch} in {halo_ind_lolo}")
                #print(halo_ind_loaded_lolo.data)
                mass = halo_ind_loaded_lolo.data.loc[halo_ind_loaded_lolo.data['ID'] == key_crossmatch].iloc[0]["Mass"]
                mstars = halo_ind_loaded_lolo.data.loc[halo_ind_loaded_lolo.data['ID'] == key_crossmatch].iloc[0]["Mstars"]
            #    print(mass, mstars)
            except:
                mass= np.nan
                mstars =np.nan


        else:
            #print(f"calculating mass! with {key_crossmatch}")
            mass = halo_ind_loaded.data.loc[halo_ind_loaded.data['ID'] == key_crossmatch].iloc[0]["Mass"]
            mstars = halo_ind_loaded.data.loc[halo_ind_loaded.data['ID'] == key_crossmatch].iloc[0]["Mstars"]

        #print(mass, mstars)
        x = halo_ref.list_halos[f"{halo}_X"][i]
        y = halo_ref.list_halos[f"{halo}_Y"][i]
        z = halo_ref.list_halos[f"{halo}_Z"][i]
        r = halo_ref.list_halos[f"{halo}_R"][i]

        new_row = {'Snapshot':name,  "Lookback":lookback[i], "X":x,"Y":y,"Z":z,"R":r, "Mass":mass, "Mstars":mstars}
           # print (new_row)
          #  new_pd = pd.DataFrame(new_row, index = name)
           # csv_dato_halo = pd.concat([csv_dato_halo, new_pd])
        csv_dato_halo = csv_dato_halo.append(new_row, ignore_index = True)
    #rocksana_snapshots_re = rocksana_snapshots[::-1]

 #   a = 1.3273e11*csv_dato_halo["Mass"]/(3.086e16**2*csv_dato_halo["R"]**2)
 #   ax.plot(lookback,a,label = f"{mass:.1e}",ls  = "--", alpha =0.6)
    csv_dato_halo.to_csv(path_rockstar_tables + f"{halo}.csv", index = 0)