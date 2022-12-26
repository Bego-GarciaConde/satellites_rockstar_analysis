import numpy as np
import pandas as pd
import heapq
import itertools
import json
from config import*


datos_edades =  pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)
lookback = []
for i,snapshot in enumerate(snapshots_analysis):
    lb = datos_edades.loc[datos_edades['Snapshot'] == snapshot, 'Lookback'].iloc[0]
    lookback.append(lb)

rocksana_snapshots = [620, 689, 720, 740, 805, 850,900, 910, 999]
rocksana_snapshots_string = ["ID","620", "689","720" "740", "805", "850","900","910", "999"]

rocksana_snapshots_re = rocksana_snapshots[::-1]  #reverse order from present to past
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
    """ Creates a table in which we can find the ID of halo and the snapshot were we found it
    """     
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

    rocksana_i = np.delete(np.array(rocksana_snapshots_re), np.where(np.array(rocksana_snapshots_re) ==rocksana_snapshots_re[0]))
    for j,comp in enumerate(rocksana_snapshots_re):
        df1 =  Rockstar(comp)
        if comp ==530:
            break
        halos_comp = Rockstar(comp)
        rocksana_i = np.delete(np.array(rocksana_i), np.where(np.array(rocksana_i) ==comp))
        print(f"for snapshot {comp} ",rocksana_i)
        
        for i, name in enumerate (rocksana_i):
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
    vetados = [433.0, 2050.0, 2055.0] #some halos that give problems
    list_of_halos = pd.DataFrame(dato_vacio, columns=['ID', 'Snapshot', "Mass","Mstars"])

    for i, halo in enumerate(halos_lista):
        halos_comp = Rockstar(rocksana_snapshots_re[i])
        for key in halo:
            print(f"saving halo {key} which is in snapshot {rocksana_snapshots_re[i]}")
            mass_stars = halos_comp.data.loc[halos_comp.data['ID'] == float(key), "Mstars"].iloc[0]
            mass = halos_comp.data.loc[halos_comp.data['ID'] == float(key), "Mass"].iloc[0]

            if key in vetados:
                pass
            else:
                new_row = {'ID':key, 'Snapshot':rocksana_snapshots_re[i], "Mass":mass,"Mstars":mass_stars }#
                list_of_halos = list_of_halos.append(new_row, ignore_index = True)

    list_of_halos.to_csv("results/list_of_halos.csv", sep = ",")


    lista_definitiva = np.array(list_of_halos["ID"])
    IDs_ref = np.array(list_of_halos["ID"])
    print(len(lista_definitiva))


    #To get the final list, we filter by orbit, since filtering by particle ID is not 100% accurate,
    #since some halos are found twice with difference in Rvir. We take those with higher mass for the same halo
    for key in IDs_ref:

        tabla_halo1 = pd.read_csv(path_rockstar_tables + f"{key}.csv")
        x = np.array(tabla_halo1["X"])
        x[np.isnan(x)] = 200

        y = np.array(tabla_halo1["Y"])
        y[np.isnan(y)] = 200

        z = np.array(tabla_halo1["Z"])
        z[np.isnan(z)] = 200

        for kk in IDs_ref:
            if key== kk:
                pass
            else:

                tabla_halo2 = pd.read_csv(path_rockstar_tables + f"{kk}.csv")
                x2 = np.array(tabla_halo2["X"])
                x2[np.isnan(x2)] = 200    #If the halo is too high in Rvir, we set its R to 200 kpc
                y2 = np.array(tabla_halo2["Y"])
                y2[np.isnan(y2)] = 200
                z2 = np.array(tabla_halo2["Z"])
                z2[np.isnan(z2)] = 200
 
                if (np.mean(np.abs(x -x2))< 1) &(np.mean(np.abs(y -y2))< 1) &(np.mean(np.abs(z -z2))< 1):
                    print(f"coincidencia de {key} con {kk}!")

                    if np.nanmean(tabla_halo2["Mass"])> np.nanmean(tabla_halo1["Mass"]):
                        print(f"eliminar {key}")
                        if key in lista_definitiva:
                            lista_definitiva = lista_definitiva[lista_definitiva != key]
                            print(len(lista_definitiva))
    
                    elif np.nanmean(tabla_halo2["Mass"])< np.nanmean(tabla_halo1["Mass"]):
                        print(f"eliminar {kk}")
                        if kk in lista_definitiva:
                            lista_definitiva = lista_definitiva[lista_definitiva != kk]
                            print(len(lista_definitiva))

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
    with open('results/lista_definitiva.txt', 'w') as f:
        for line in lista_definitiva:
            f.write(f"{line}\n")

def create_tables_coordinates(list_of_halos, datos_crossmatch):
    """Creates a table of coordinates comparing the same halo with different snapshots in which it was found
    """
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
                    halo_ind_iteration = rocksana_snapshots[ind[i+cont]]
                    halo_ind_loaded_iteration= Rockstar(int(halo_ind_iteration))
                #  print(f"calculating mass! with {key_crossmatch} in {halo_ind_iteration}")
                    #print(halo_ind_loaded_iteration.data)
                    mass = halo_ind_loaded_iteration.data.loc[halo_ind_loaded_iteration.data['ID'] == key_crossmatch].iloc[0]["Mass"]
                    mstars = halo_ind_loaded_iteration.data.loc[halo_ind_loaded_iteration.data['ID'] == key_crossmatch].iloc[0]["Mstars"]
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
