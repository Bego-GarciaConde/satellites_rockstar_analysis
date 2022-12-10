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
from config import*


datos_edades =  pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)
lookback = []
for i,snapshot in enumerate(snapshots_analysis):
    lb = datos_edades.loc[datos_edades['Snapshot'] == snapshot, 'Lookback'].iloc[0]
    lookback.append(lb)

rocksana_snapshots = [620, 689, 720, 740, 805, 850,900, 910, 999]
rocksana_snapshots_string = ["ID","620", "689","720" "740", "805", "850","900","910", "999"]
#for i in rocksana_snapshots:
#   rocksana_snapshots_string.append(str(i))
rocksana_snapshots_re = rocksana_snapshots[::-1]  #reverso order from present to past
print(rocksana_snapshots_re)

class Rockstar:
    def __init__(self, name):
        self.name = name
        self.data = None
        self.list_halos = None
        self.id_list = None
        def read_data():
            print(f"Initializing {self.name}")
            self.data= pd.read_csv(f"results/DF_halos_{self.name}.csv", index_col = 0)
            with open(f"results/coordenadas_de_halos_{self.name}.json") as json_file:
                self.list_halos=json.load(json_file)
            with open(f"results/halos_IDs_list_{self.name}.json") as json_file:
                self.id_list=json.load(json_file)
            for key in self.data["ID"]:
                # print(np.mean(data_805[f"{key}_R"]))
                if np.mean(self.list_halos[f"{key}_R"]) < 2:
                    self.data= self.data[self.data.ID != key]

                if np.isnan(self.list_halos[f"{key}_R"]).all() == True:
                    self.data = self.data[self.data.ID != key]


        read_data()
def create_list():         



    halos_lista = [] #we will modify this list
    halos_lista_ref = []
    for comp in rocksana_snapshots_re:
        halos_comp = Rockstar(comp)
        halos_lista.append(np.array(halos_comp.data["ID"]))
        halos_lista_ref.append(np.array(halos_comp.data["ID"]))

    dict_halos = {}
    for i, name in enumerate(rocksana_snapshots_re):
        dict_halos[f"{name}"]= i
    print(dict_halos)


    contador = 0
    rocksana_i = np.delete(np.array(rocksana_snapshots_re), np.where(np.array(rocksana_snapshots_re) ==rocksana_snapshots_re[0]))
    for j,comp in enumerate(rocksana_snapshots_re):
        df1 =  Rockstar(comp)
        if comp ==530:
            break
        halos_comp = Rockstar(comp)
        rocksana_i = np.delete(np.array(rocksana_i), np.where(np.array(rocksana_i) ==comp))
        print(f"for snapshot {comp} ",rocksana_i)
        
        for i, name in enumerate (rocksana_i):
            #df1 =rocksana[0]
        # list_of_pairs[f"530"] = np.nan
            df2 =  Rockstar(name)
            for key in halos_lista_ref[j]:
                start = 2
        
                while df1.list_halos[f"{key}_R"][start] is np.nan:
                    start = start + 10

                x = df1.list_halos[f"{key}_X"][start]
                y = df1.list_halos[f"{key}_Y"][start]
                z = df1.list_halos[f"{key}_Z"][start]
                r =  np.array(df1.list_halos[f"{key}_R"])
                mass_stars = halos_comp.data.loc[halos_comp.data['ID'] == float(key), "Mstars"].iloc[0]
                for kk in df2.data["ID"]:
                    
                # print(f"comparing {key} in 530 and {kk} in 999 ")
                #  if df2.list_halos[f"{kk}_R"][start] is not np.nan:
                    xx = df2.list_halos[f"{kk}_X"][start]
                    yy = df2.list_halos[f"{kk}_Y"][start]
                    zz = df2.list_halos[f"{kk}_Z"][start]
                    rr = np.array(df2.list_halos[f"{kk}_R"])
                    
                    
                    if mass_stars ==0:
                        reference = df1.id_list[f"{key}"]
                        comparison = df2.id_list[f"{kk}"]
                        res = (np.in1d(reference,comparison)==True)
                        
                        if res.sum()/len(res)>0.65:
                            halos_lista[dict_halos[f"{name}"]] = np.delete(np.array(halos_lista[dict_halos[f"{name}"]]), np.where(np.array(halos_lista[dict_halos[f"{name}"]]) ==kk))
                    elif mass_stars >0:
                        reference = df1.id_list[f"{key}"]
                        comparison = df2.id_list[f"{kk}"]
                        res = (np.in1d(reference,comparison)==True)
                        if  res.sum()/len(res)>0.65:
                            halos_lista[dict_halos[f"{name}"]] = np.delete(np.array(halos_lista[dict_halos[f"{name}"]]), np.where(np.array(halos_lista[dict_halos[f"{name}"]]) ==kk))

    dato_vacio = []
    vetados = [433.0, 2050.0, 2055.0]
    list_of_halos = pd.DataFrame(dato_vacio, columns=['ID', 'Snapshot', "Mass","Mstars"])
    for i, halo in enumerate(halos_lista):
        halos_comp = Rockstar(rocksana_snapshots_re[i])
        for key in halo:
            print(f"saving halo {key} which is in snapshot {rocksana_snapshots_re[i]}")
            mass_stars = halos_comp.data.loc[halos_comp.data['ID'] == float(key), "Mstars"].iloc[0]
            mass = halos_comp.data.loc[halos_comp.data['ID'] == float(key), "Mass"].iloc[0]
        #  new_row= [key,rocksana_snapshots_re[i], mass, mass_stars ]
        # print(new_row)
            if key in vetados:
                pass
            else:
                new_row = {'ID':key, 'Snapshot':rocksana_snapshots_re[i], "Mass":mass,"Mstars":mass_stars }#
                list_of_halos = list_of_halos.append(new_row, ignore_index = True)

    list_of_halos.to_csv("results/list_of_halos.csv", sep = ",")


    lista_definitiva = np.array(list_of_halos["ID"])
    IDs_ref = np.array(list_of_halos["ID"])
    print(len(lista_definitiva))
    #IDs_ref =[435.0, 1418.0]
    for key in IDs_ref:
    #key = 685.0

        tabla_halo1 = pd.read_csv(path_rockstar_tables + f"{key}.csv")
        x = np.array(tabla_halo1["X"])
        x[np.isnan(x)] = 200
        #x = x[~np.isnan(x)]
        y = np.array(tabla_halo1["Y"])
        y[np.isnan(y)] = 200
    #  y = y[~np.isnan(y)]
        z = np.array(tabla_halo1["Z"])
        z[np.isnan(z)] = 200
    # z = z[~np.isnan(z)]
        for kk in IDs_ref:
            if key== kk:
                pass
            else:

                tabla_halo2 = pd.read_csv(path_rockstar_tables + f"{kk}.csv")
                x2 = np.array(tabla_halo2["X"])
            #  x2 = x2[~np.isnan(x2)]
                x2[np.isnan(x2)] = 200
                y2 = np.array(tabla_halo2["Y"])
            #  y2 = y2[~np.isnan(y2)]
                y2[np.isnan(y2)] = 200
                z2 = np.array(tabla_halo2["Z"])
                z2[np.isnan(z2)] = 200
                #z2 = z2[~np.isnan(z2)]
                #if len(x) == len(x2):
                if (np.mean(np.abs(x -x2))< 1) &(np.mean(np.abs(y -y2))< 1) &(np.mean(np.abs(z -z2))< 1):
                    print(f"coincidencia de {key} con {kk}!")
                #  print(np.nanmean(tabla_halo2["Mass"]),np.nanmean(tabla_halo1["Mass"]) )
                    if np.nanmean(tabla_halo2["Mass"])> np.nanmean(tabla_halo1["Mass"]):
                        print(f"eliminar {key}")
                        if key in lista_definitiva:
                            lista_definitiva = lista_definitiva[lista_definitiva != key]
                            print(len(lista_definitiva))
                        #  lista_definitiva.delete(key)
                    elif np.nanmean(tabla_halo2["Mass"])< np.nanmean(tabla_halo1["Mass"]):
                        print(f"eliminar {kk}")
                        if kk in lista_definitiva:
                            lista_definitiva = lista_definitiva[lista_definitiva != kk]
                            print(len(lista_definitiva))
                        # lista_definitiva = lista_definitiva.delete(kk)
                    elif np.nanmean(tabla_halo2["Mass"])== np.nanmean(tabla_halo1["Mass"]):
                        if kk<key:
                            if kk in lista_definitiva:
                                print(f"eliminar {kk}")
                                lista_definitiva = lista_definitiva[lista_definitiva != kk]
                                print(len(lista_definitiva))
                        else:
                            if key in lista_definitiva:
                                print(f"eliminar {key}")
                                lista_definitiva = lista_definitiva[lista_definitiva != key]
                                print(len(lista_definitiva))
                            
                            
                else:
                    
                    pass



    print(f"Definitive list of halos is {lista_definitiva}")

def create_tables_coordinates(list_of_halos, datos_crossmatch):
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

            csv_dato_halo = csv_dato_halo.append(new_row, ignore_index = True)

        csv_dato_halo.to_csv(path_rockstar_tables + f"{halo}.csv", index = 0)
