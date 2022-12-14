
import pandas as pd
import numpy as np

path_halos = "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/Rockstar_GARROTXA/"
path_yt = "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/yt_snapshots/"
path_datos = "/home/bego/GARROTXA_copia/datos_GARROTXA_resim/"
path_csv =  "/mnt/usb-TOSHIBA_EXTERNAL_USB_20220124010088F-0:0-part2/snapshots_resim_new/"
datos_edades = pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)

path_rockstar_data = "/home/bego/GARROTXA/satellites_rockstar_analysis/results/"
path_rockstar_tables = "/home/bego/GARROTXA/satellites_rockstar_analysis/satelites_tables/"



#----This code operates in two modes----
process_rockstar_snapshot = 0 #Takes a new snapshot of rockstar and gives the satellites data
process_all_satellites = 0 #Takes all processed rockstar snapshots and identifies satellites, crossmatches them through all processed snapshots


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

satelites_rockstar = [275, 359, 375, 400, 425, 450, 475, 500, 505, 510, 515, 520, 523, 525, 527, 530,532,
535, 537, 539,541,543, 545,547, 550, 553, 555, 557, 560,563, 565, 567, 570, 573, 575, 577, 580, 583, 585, 587, 
590, 600, 610, 620, 629, 630, 639, 640, 650, 660, 670, 679, 687, 689, 690, 720, 739, 740, 755,
763, 770, 780,790,797, 805, 810, 820, 830, 840, 850, 853, 855, 860, 867, 870, 875, 877, 883, 890, 900, 910,
929, 930, 939, 950, 963, 970, 979, 980, 989, 990, 993, 999]