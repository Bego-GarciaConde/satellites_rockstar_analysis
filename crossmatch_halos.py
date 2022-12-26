import numpy as np
import pandas as pd


import heapq

import json
from config import*


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



def  crossmatch_between_data ( df1, df2,namedf1,namedf2, df_halos):
    """This function identifies the same halo in different snapshots by crossmatching particles"""
    for key in df1.data["ID"]:
        start = 2
        mass_stars =  df1.data.loc[ df1.data['ID'] == float(key), "Mstars"].iloc[0]

        while df1.list_halos[f"{key}_R"][start] is np.nan:
            start = start + 10

        for kk in df2.data["ID"]:

            if mass_stars ==0:
                reference = df1.id_list[f"{key}"]
                comparison = df2.id_list[f"{kk}"]
                res = (np.in1d(reference,comparison)==True)
                if res.sum()/len(res)>0.65:
                    df_halos.loc[key,f"{namedf2}"] = kk
                    break
                else:
                    df_halos.loc[key,f"{namedf2}"] = np.nan
   
            elif mass_stars >0:
                reference = df1.id_list[f"{key}"]
                comparison = df2.id_list[f"{kk}"]
                res = (np.in1d(reference,comparison)==True)
                if  res.sum()/len(res)>0.65:
                    df_halos.loc[key,f"{namedf2}"] = kk
                    break
                else:
                    df_halos.loc[key,f"{namedf2}"] = np.nan


                    df_halos.loc[key,f"{namedf2}"] = np.nan

def crossmatch_halos():
    """Takes rocksana snapshots one by one and compares"""

    halos_recompilation = []

    rocksana_snapshots = [620, 689, 720, 740, 805, 850,900, 910, 999]
    rocksana_snapshots_string = ["ID","620", "689","720", "740", "805", "850","900","910", "999"]
    #for i in rocksana_snapshots:
    #   rocksana_snapshots_string.append(str(i))
    rocksana_snapshots_re = rocksana_snapshots[::-1]
    print(rocksana_snapshots_re)
    for comp in rocksana_snapshots_re:
        #comp = 999
        print(comp)
        dato_vacio = {}
        comparison = comp
        halos_comp = Rockstar(comparison)
        for snapshot in rocksana_snapshots_string:
            dato_vacio [f"{snapshot}"] = np.array(halos_comp.data["ID"])

        df_halos = pd.DataFrame(dato_vacio, columns=rocksana_snapshots_string)
        df_halos.set_index('ID',inplace = True)

        rocksana_i = np.delete(np.array(rocksana_snapshots), np.where(np.array(rocksana_snapshots) ==comp))
        print(f"for snapshot {comp} ",rocksana_i)
        for i, name in enumerate (rocksana_i):

            df2 =  Rockstar(name)
            crossmatch_between_data ( halos_comp, df2, comparison, name, df_halos)
        halos_recompilation.append(df_halos)



    recomp_definitive = pd.concat(halos_recompilation)
    #recomp_definitive.reset_index(drop=True)
    recomp_def = recomp_definitive.drop_duplicates(keep = "first")
    recomp_def.to_csv("results/crossmatch_of_halos.csv", sep = ",")