import numpy 
import scipy.spatial
import sys
import math
import csv
import random
import json 
import time
import copy 
import startin


numpy.set_printoptions(threshold=sys.maxsize)

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

    clean_points_list = []
    for point1 in list_pts_3d:
	    repeated = False
	    for point2 in clean_points_list:
		    if point1[0] == point2[0] and point1[1] == point2[1]:
			    repeated = True
	    if repeated == False:
		    clean_points_list.append(point1)
	    else:
		    print("Repeated point: " + str(point1[0]) + " " + str(point1[1]))
    
    list_pts = (clean_points_list)
    
    gridsize= int(jparams['kriging']['cellsize'])
    radius = jparams['kriging']['radius']
    var_range = 260
    var_nugget = 18
    var_sill = 1300

    

    def distance(point1, point2):
	    return math.sqrt((point2[0]-point1[0])*(point2[0]-point1[0])+(point2[1]-point1[1])*(point2[1]-point1[1]))

    def variogram(h):
        r = var_range
        n = var_nugget
        s = var_sill
        h_theor = n + s * (1 - math.exp(((-9 * (h * h) ))/(r * r)))
        return h_theor
    

    
    #copy 3d list before splitting it
    #list_pts = copy.deepcopy(list_pts_3d)
    x = []
    y = []
    z = []

    #splitting list
    for point in list_pts_3d:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    #removing z values
    for pt in list_pts:
        pt.pop(2)
    #storing xy vals in kdtree
    kd = scipy.spatial.KDTree(list_pts)
    convex_hull = scipy.spatial.Delaunay(clean_points_list)
    #determining grid size based on 
    ncols = math.ceil((max(x)-min(x))/gridsize)
    nrows = math.ceil((max(y)-min(y))/gridsize)
    
    yrange = reversed(numpy.arange(min(y) + 0.5 * gridsize, max(y) + 0.5 * gridsize, gridsize))
    xrange = numpy.arange(min(x) + 0.5 * gridsize, max(x) + 0.5 * gridsize, gridsize)
    
    coordinates = [[i, j] for j in yrange for i in xrange]
    for i in coordinates:
        if convex_hull.find_simplex(i) == -1:
            i.append(-9999)
        else:
            radius_points = kd.query_ball_point(i, radius)
            sample_points = []
            query_distances = []
            z_vals = []
            for points in radius_points:
                pointe = [x[points],y[points]]
                sample_points.append(pointe)
            distance_matrix  = (scipy.spatial.distance_matrix(sample_points, sample_points))
            
            for ii in range(len(distance_matrix)):
                for jj in range((ii)):
                    if distance_matrix[ii][jj] != 0:
                        distance_matrix[ii][jj] = variogram(distance_matrix[ii][jj])
            c = numpy.matrix(distance_matrix)
            for points in radius_points:
                dist = distance(i, [x[points],y[points]])
                query_distances.append(dist)
                z_vals.append(z[points])
            c_trans = numpy.linalg.inv(c)
            weight = numpy.dot(c_trans, query_distances)
            z_new = numpy.dot(weight, z_vals)
            i.append(z_new)                        

                    
                    
            
    row_num = 0
    col_num = 0



    with open(j_idw['output-file'], 'w') as fh:
        fh.writelines('NCOLS {}\n'.format(ncols))
        fh.writelines('NROWS {}\n'.format(nrows))
        fh.writelines('XLLCENTER {}\n'.format(min(x)))
        fh.writelines('YLLCENTER {}\n'.format(min(y)))
        fh.writelines('CELLSIZE {}\n'.format(j_idw['cellsize']))
        fh.writelines('NODATA_VALUE {}\n'.format(-9999))
        for i in coordinates:
            fh.write(str(i[-1])+' ')
            col_num += 1
            if col_num  == ncols:
                col_num = 0
                row_num+= 1
                fh.write('\n')


if __name__ == '__main__':
    main()

