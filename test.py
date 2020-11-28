import numpy as np
import scipy.spatial
import sys
import math
import csv
import random
import json 
import time
import copy 



def main():
    #-- read the needed parameters from the file 'params.json' (must be in same folder)
    try:
        jparams = json.load(open('params.json'))
    except:
        print("ERROR: something is wrong with the params.json file.")
        sys.exit()
    #-- store the input 3D points in list
    list_pts_3d = []
    with open(jparams['input-file']) as csvfile:
        r = csv.reader(csvfile, delimiter=' ')
        header = next(r)
        for line in r:
            p = list(map(float, line)) #-- convert each str to a float
            assert(len(p) == 3)
            list_pts_3d.append(p)
    gridsize= jparams['nn']['cellsize']
    list_pts = copy.copy(list_pts_3d)
    x = []
    y = []
    z = []
    sample_size = 0
    for point in list_pts_3d:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
        sample_size += 1 
    for pt in list_pts:
        pt.pop(2)
    kd = scipy.spatial.KDTree(list_pts)
    convex_hull = scipy.spatial.Delaunay(list_pts)
    print(convex_hull.find_simplex)
    ncols = math.ceil((max(x)-min(x)+(0.5 * gridsize))/gridsize)
    nrows = math.ceil((max(y)-min(y)+(0.5* gridsize))/gridsize)
    print(sample_size)
    yrange = reversed(range((int(min(y))),(int(max(y))+(gridsize)),gridsize))
    xrange = (range(int(min(x)),int(max(x)+(gridsize)),gridsize))
    coordinates = [[i, j] for j in yrange for i in xrange]
    for i in coordinates:
        if convex_hull.find_simplex(i) == -1:
            i.append(-9999)
        else:
            d, i_nn = kd.query(i,k=1)
            i.append(z[i_nn])
    row_num = 0
    col_num = 0
    with open('tas_test.asc', 'w') as fh:
        fh.writelines('NCOLS {}\n'.format(ncols))
        fh.writelines('NROWS {}\n'.format(nrows))
        fh.writelines('XLLCENTER {}\n'.format(min(x)))
        fh.writelines('YLLCENTER {}\n'.format(min(y)))
        fh.writelines('CELLSIZE {}\n'.format(jparams['nn']['cellsize']))
        fh.writelines('NODATA_VALUE {}\n'.format(-9999))
        for i in coordinates:
            fh.write(str(i[-1])+' ')
            col_num += 1
            if col_num  == ncols:
                col_num = 0
                row_num+= 1
                fh.write('\n')
3