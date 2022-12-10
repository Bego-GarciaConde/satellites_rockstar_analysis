
from config import *
import yt
from yt import YTArray
from useful_functions import*

from snapshot_definition import *

import json


class rockstar_snapshot:
    def __init__(self, number,name):
        self.number = number
        self.name = name       
        self.center = None
        self.mA = None
        self.mB = None
     #   self.yt_data = None
        self.lb = None
        self.halo_data = None
        self.dm_ID_list = None

      #  self.dm_ID_list = None
        



        def read_lb():
            self.lb = datos_edades.loc[datos_edades['Snapshot'] == self.name, 'Lookback'].iloc[0]
        def read_center():
            centro = np.loadtxt(path_datos +f'center_{self.name}.txt')
            center = YTArray([centro[0], centro[1], centro[2]], "cm")
            self.center = center


        def find_rotation_matrices ():
            self.mA =  np.loadtxt(path_datos + f"rotation_matrix/{name}_mA.txt")
            self.mB =  np.loadtxt(path_datos + f"rotation_matrix/{name}_mB.txt")

        read_lb()
        read_center()
        find_rotation_matrices ()

    def read_halos_data(self, snapshot):

        #Read halos with yt
        halos_ds = yt.load(path_halos + f"halos_{self.number}.0.bin")

        sp2=halos_ds.sphere(self.center,(200, "kpc"))

        x_pos = np.array(sp2["halos", "particle_position_x"].in_units("kpc") - self.center[0].in_units("kpc") )
        y_pos = np.array(sp2["halos", "particle_position_y"].in_units("kpc") - self.center[1].in_units("kpc") )
        z_pos = np.array(sp2["halos", "particle_position_z"].in_units("kpc") - self.center[2].in_units("kpc") )
        mass = np.array(sp2["halos", "particle_mass"].in_units("Msun"))
        Rvir = np.array(sp2["halos", "virial_radius"].in_units("kpc"))
        ID = np.array(sp2["halos", "particle_identifier"])
    
        ind = np.where((mass>1e6) &(mass<1e11))

        #Apply transformation matrices which situates in reference system of the galactic plane at z = 0
        x_re, y_re, z_re = apply_transformation_matrix(self.mA, x_pos[ind],y_pos[ind],z_pos[ind])
        x_rre, y_rre, z_rre = apply_transformation_matrix(self.mB, x_re,y_re,z_re)

        data = {"ID": ID[ind],"X":x_rre, "Y":y_rre, "Z":  z_rre, "Mass": mass[ind], "Rvir": Rvir[ind] }
        df_halos = pd.DataFrame(data = data)

        #Some halos are redundant, we discard the repeated ones
        dato_vacio = []
        df_vacio = pd.DataFrame(dato_vacio, columns=["ID","X", "Y", "Z", "Mass", "Rvir"])
        #Assign new IDS
        for i,halo_id in enumerate(df_halos['ID']):
            x_model, y_model, z_model  =df_halos.loc[df_halos['ID'] == halo_id, ['X', "Y", "Z"]].iloc[0]
            rvir_model = df_halos.loc[df_halos['ID'] == halo_id, "Rvir"].iloc[0]
                
            subgroup = np.where(np.sqrt((df_halos["X"]-x_model)**2  +(df_halos["Y"]-y_model)**2+(df_halos["Z"]-z_model)**2)<0.2*rvir_model)
            sg = df_halos.iloc[subgroup]
   
            if len(sg) >=2:
                df_halos.loc[subgroup[0],"new_ID"]=halo_id
  
            else:
                df_halos.loc[i,"new_ID"]=halo_id
                
        test_list = list(set(df_halos["new_ID"]))
        print(len(df_halos["new_ID"]),len(test_list))
       # df_index = []
        #Discard the redundant halos
        for i,halo_id in enumerate(test_list):
    
            sg = df_halos.loc[df_halos['new_ID'] == halo_id, ]
            if len(sg)>2:
                sg = sg.loc[sg['Mass'].idxmax(), :]   #If the halo is repeated, take the one with higher mass
            df_vacio.append(sg)
            
            df_vacio = df_vacio.append(sg)
            df_vacio = df_vacio.sort_index( ascending = True)
        #    df_index.append(i)

        df_vacio["index"]=np.arange(start=0, stop=len(df_vacio))
        df_vacio.set_index("index", inplace=True)

        self.halo_data = df_vacio

    #def find_IDs_particles (self): 

        ID_dm_list={}
        ID_stars_list={}
        #star_mass = 
        star_mass =[]
        df_index = []
        for i,halo_id in enumerate(self.halo_data['ID']):
            print(i, halo_id)
            x_model, y_model, z_model  =self.halo_data.loc[self.halo_data['ID'] == halo_id, ['X', "Y", "Z"]].iloc[0]
            rvir_model = self.halo_data.loc[self.halo_data['ID'] == halo_id, "Rvir"].iloc[0]
            
            ind = np.where(np.sqrt((snapshot.dm["X"]-x_model)**2  +(snapshot.dm["Y"]-y_model)**2+(snapshot.dm["Z"]-z_model)**2)<0.2*rvir_model)
            #print(len(ind))
            model = snapshot.dm.iloc[ind]
            print("This halo is located at ",x_model, y_model, z_model)
            print("number of particles, ",len(model))
            ID_dm_list[f"{halo_id}"]=np.array(model["ID"])
            
            ind = np.where(np.sqrt((snapshot.stars["X"]-x_model)**2  +(snapshot.stars["Y"]-y_model)**2+(snapshot.stars["Z"]-z_model)**2)<0.2*rvir_model)
            model = snapshot.stars.iloc[ind]
            ID_stars_list[f"{halo_id}"]=np.array(model["ID"])
            star_mass.append(np.sum(model["Mass"]))
            df_index.append(i)
            
        self.halo_data["Mstars"]=np.array(star_mass)
        self.halo_data["index"]=np.array(df_index)
        self.halo_data.set_index("index", inplace=True)
        self.dm_ID_list = ID_dm_list
    
    def save_halos_info (self):
        self.halo_data.to_csv(f"results/DF_halos_{self.name}.csv")

    def save_ID_info (self):
        for key, value in self.dm_ID_list.items():
            self.dm_ID_list[key]=value.tolist()
        with open(f"results/halos_IDs_list_{self.name}.json", "w") as outfile:
            json.dump(self.dm_ID_list, outfile)

    def find_coordinates_by_IDs(self):
        coord = ["X", "Y", "Z", "R"]
        coordinates={}
        for key, value in self.dm_ID_list.items():
            for co in coord:
                coordinates[f"{key}_{co}"] = []


        for name in snapshots_analysis:
            print(name)
            snapshot = Snapshot(name)
            snapshot.load_stars()
            snapshot.load_dm()
            snapshot.dm = pd.read_csv(path_csv + f"{name}_dm_Rvir.csv", sep = ",")
            for key, value in self.dm_ID_list.items():
                dfA =snapshot.dm[snapshot.dm['ID'].isin(value)]
                if len(dfA)> 0.5*len(value):
                    x= np.median(dfA["X"])
                    y = np.median(dfA["Y"])
                    z = np.median(dfA["Z"])
                    r = np.sqrt(x*x + y*y + z*z)

                else:
                    x = np.nan
                    y = np.nan
                    z = np.nan
                    r = np.nan
                #   print(f"Satellite {key} at {name} not found")

                coordinates[f"{key}_X"].append(x)
                coordinates[f"{key}_Y"].append(y)
                coordinates[f"{key}_Z"].append(z)
                coordinates[f"{key}_R"].append(r)
                
        with open(f"results/coordenadas_de_halos_{self.name}.json", "w") as outfile:
            json.dump(coordinates, outfile)
                
   