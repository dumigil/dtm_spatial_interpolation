import numpy
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
    gridsize= int(jparams['tin']['cellsize'])
    list_pts = copy.deepcopy(list_pts_3d)
    x = []
    y = []
    z = []
    for point in list_pts_3d:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    for pt in list_pts:
        pt.pop(2)
    convex_hull = scipy.spatial.Delaunay(list_pts)
    ncols = math.ceil((max(x)-min(x))/gridsize)
    nrows = math.ceil((max(y)-min(y))/gridsize)

    yrange = reversed(numpy.arange(min(y) + 0.5 * gridsize, max(y) + 0.5 * gridsize, gridsize))
    xrange = numpy.arange(min(x) + 0.5 * gridsize, max(x) + 0.5 * gridsize, gridsize)

    coordinates = [[i, j] for j in yrange for i in xrange]
    for i in coordinates:
        if convex_hull.find_simplex(i) == -1:
            i.append(-9999)
        else:
            v1, v2, v3 = convex_hull.simplices[convex_hull.find_simplex(i)]
            p1 = list_pts[v1]
            p2 = list_pts[v2]
            p3 = list_pts[v3]
            if xy_distance(i, p1) == 0:
                i.append(z[v1])
            elif xy_distance(i, p2) == 0:
                i.append(z[v2])
            elif xy_distance(i, p3) == 0:
                i.append(z[v3])
            else:
                weight_1 = 1/xy_distance(i, p1)
                weight_2 = 1/xy_distance(i, p2)
                weight_3 = 1/xy_distance(i, p3)
                i.append(((weight_1 * z[v1]) + (weight_2 * z[v2]) + (weight_3 * z[v3])) / (weight_1 + weight_2 + weight_3))

    
    row_num = 0
    col_num = 0
    with open('tas_tin_test.asc', 'w') as fh:
        fh.writelines('NCOLS {}\n'.format(ncols))
        fh.writelines('NROWS {}\n'.format(nrows))
        fh.writelines('XLLCENTER {}\n'.format(min(x)))
        fh.writelines('YLLCENTER {}\n'.format(min(y)))
        fh.writelines('CELLSIZE {}\n'.format(jparams['tin']['cellsize']))
        fh.writelines('NODATA_VALUE {}\n'.format(-9999))
        for i in coordinates:
            fh.write(str(i[-1])+' ')
            col_num += 1
            if col_num  == ncols:
                col_num = 0
                row_num+= 1
                fh.write('\n')

def xy_distance(point1, point2):
    dx = point1[0] - point2[0]
    dy = point1[1] - point2[1]
    distance = math.sqrt(dx **2 + dy **2)
    return distance


if __name__ == '__main__':
    main()