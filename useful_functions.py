


def apply_transformation_matrix(matriz_trans, X,Y,Z):
    X_re= X*matriz_trans[0,0]+ Y*matriz_trans[0,1] + Z*matriz_trans[0,2]       
    Y_re= X*matriz_trans[1,0]+ Y*matriz_trans[1,1] + Z*matriz_trans[1,2]
    Z_re= X*matriz_trans[2,0]+ Y*matriz_trans[2,1] + Z*matriz_trans[2,2]
    return X_re, Y_re, Z_re