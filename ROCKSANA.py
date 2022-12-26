#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 05/10/2022
@author: B. García-Conde
"""


#from fourier_analysis_functions import *
#from accelerations import * 
from config import *
from satellite_definition import *
from snapshot_definition import *

from snapshot_definition import Snapshot

from create_list import*
from crossmatch_halos import*
#import ConfigParser


        ##################################################
        #             ROCKSTAR's satellites ANLYSIS                 #
    #
        ##################################################

def main ():
        if process_rockstar_snapshot ==1:
                snapshot = Snapshot(620) #We indicate the snaphot
                snapshot.load_stars()
                snapshot.load_dm()
                sat = rockstar_snapshot(43, 620)
                sat.read_halos_data(snapshot)
                sat.save_halos_info()
                sat.save_ID_info()
                sat.find_coordinates_by_IDs ()
        if process_all_satellites ==1:
                create_list()
                crossmatch_halos()
                list_of_halos =  pd.read_csv("satellites_rockstar_analysis/results/list_of_halos.csv")
                datos_crossmatch = pd.read_csv( path_rockstar_data + "crossmatch_of_halos.csv", sep = ",")
                create_tables_coordinates(list_of_halos, datos_crossmatch)

if __name__ == "__main__":
    main()
    