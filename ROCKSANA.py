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
#import ConfigParser


        ##################################################
        #             ROCKSTAR's satellites ANLYSIS                 #
    #
        ##################################################

def main ():
        snapshot = Snapshot(530)
        snapshot.load_stars()
        snapshot.load_dm()
        sat = rockstar_snapshot(15, 530, snapshot)
        sat.read_halos_data()
        sat.save_halos_info()
        sat.save_ID_info()
        sat.find_coordinates_by_IDs


if __name__ == "__main__":
    main()