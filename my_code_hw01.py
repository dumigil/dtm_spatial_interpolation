#-- my_code_hw01.py
#-- hw01 GEO1015.2020
#-- Michiel de Jong
#-- 4376978
#-- [YOUR NAME]
#-- [YOUR STUDENT NUMBER] 


#-- import outside the standard Python library are not allowed, just those:
import math
import numpy
import scipy.spatial
import startin 
import copy
#-----
def xy_distance(point1, point2):
        dx = point1[0] - point2[0]
        dy = point1[1] - point2[1]
        distance = math.sqrt(dx **2 + dy **2)
        return distance


def nn_interpolation(list_pts_3d, j_nn):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with nearest neighbour interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_nn:        the parameters of the input for "nn"
    Output:
        returns the value of the area
 
    """  
    # print("cellsize:", j_nn['cellsize'])

    #-- to speed up the nearest neighbour us a kd-tree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
    # kd = scipy.spatial.KDTree(list_pts)
    # d, i = kd.query(p, k=1)
    gridsize= int(j_nn['cellsize'])
    
    #copy 3d list before splitting it
    list_pts = copy.deepcopy(list_pts_3d)
    x = []
    y = []
    z = []
    
    #splitting the list in xyz
    for point in list_pts_3d:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    #removing z coordinates
    for pt in list_pts:
        pt.pop(2)
    
    #adding xy points to KD tree
    kd = scipy.spatial.KDTree(list_pts)
    
    #creating convex hull
    convex_hull = scipy.spatial.Delaunay(list_pts)
    
    #defining the number of rows and columns, taking into account that the grid is based on the centerpoints.
    ncols = math.ceil((max(x)-min(x))/gridsize)
    nrows = math.ceil((max(y)-min(y))/gridsize)

    #defining the range for the grid, reversing the y axis because the grid has to be initialised bottom left 
    yrange = reversed(numpy.arange(min(y) + 0.5 * gridsize, max(y) + 0.5 * gridsize, gridsize))
    xrange = numpy.arange(min(x) + 0.5 * gridsize, max(x) + 0.5 * gridsize, gridsize)
    
    #creating the grid array
    coordinates = [[i, j] for j in yrange for i in xrange]

    #interpolation actually happens here, first checking if point is inside convex hull, otherwise it will be no_data
    for i in coordinates:
        if convex_hull.find_simplex(i) == -1:
            i.append(-9999)
        else:
            d, i_nn = kd.query(i,k=1)
            i.append(z[i_nn])
    
    row_num = 0
    col_num = 0

    #writing file
    with open(j_nn['output-file'], 'w') as fh:
        fh.writelines('NCOLS {}\n'.format(ncols))
        fh.writelines('NROWS {}\n'.format(nrows))
        fh.writelines('XLLCENTER {}\n'.format(min(x)))
        fh.writelines('YLLCENTER {}\n'.format(min(y)))
        fh.writelines('CELLSIZE {}\n'.format(j_nn['cellsize']))
        fh.writelines('NODATA_VALUE {}\n'.format(-9999))
        for i in coordinates:
            fh.write(str(i[-1])+' ')
            col_num += 1
            if col_num  == ncols:
                col_num = 0
                row_num+= 1
                fh.write('\n')
    print("File written to", j_nn['output-file'])




def idw_interpolation(list_pts_3d, j_idw):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with IDW
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_idw:       the parameters of the input for "idw"
    Output:
        returns the value of the area
 
    """  
    # print("cellsize:", j_idw['cellsize'])
    # print("radius:", j_idw['radius'])

    #-- to speed up the nearest neighbour us a kd-tree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
    # kd = scipy.spatial.KDTree(list_pts)
    # i = kd.query_ball_point(p, radius)
    
    #pulling parameters from json
    gridsize= int(j_idw['cellsize'])
    radius = j_idw['radius']
    power = j_idw['power']

    #copy 3d list before splitting it
    list_pts = copy.deepcopy(list_pts_3d)
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
    convex_hull = scipy.spatial.Delaunay(list_pts)
    
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
            d, i_nn = kd.query(i,k=1)
            i_idw = kd.query_ball_point(i, radius)
            weight_sum = 0
            z_sum = 0
            for idx in i_idw:
                distance = math.sqrt(((x[idx]-i[0])**2)+((y[idx]-i[1])**2))
                if distance != 0:
                    weight = 1/(distance ** power)
                    z_val = z[idx] * weight
                    weight_sum += weight
                    z_sum += z_val
            i.append((z_sum/weight_sum))
            
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
    print("File written to", j_idw['output-file'])


def tin_interpolation(list_pts_3d, j_tin):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with linear in TIN interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_tin:       the parameters of the input for "tin"
    Output:
        returns the value of the area
 
    """  
    #-- example to construct the DT with scipy
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html#scipy.spatial.Delaunay
    # dt = scipy.spatial.Delaunay([])

    #-- example to construct the DT with startin
    # minimal docs: https://github.com/hugoledoux/startin_python/blob/master/docs/doc.md
    # how to use it: https://github.com/hugoledoux/startin_python#a-full-simple-example
    # you are *not* allowed to use the function for the tin linear interpolation that I wrote for startin
    # you need to write your own code for this step
    # but you can of course read the code [dt.interpolate_tin_linear(x, y)]
    
    gridsize= int(j_tin['cellsize'])
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
    with open(j_tin['output-file'], 'w') as fh:
        fh.writelines('NCOLS {}\n'.format(ncols))
        fh.writelines('NROWS {}\n'.format(nrows))
        fh.writelines('XLLCENTER {}\n'.format(min(x)))
        fh.writelines('YLLCENTER {}\n'.format(min(y)))
        fh.writelines('CELLSIZE {}\n'.format(j_tin['cellsize']))
        fh.writelines('NODATA_VALUE {}\n'.format(-9999))
        for i in coordinates:
            fh.write(str(i[-1])+' ')
            col_num += 1
            if col_num  == ncols:
                col_num = 0
                row_num+= 1
                fh.write('\n')

    print("File written to", j_tin['output-file'])


def kriging_interpolation(list_pts_3d, j_kriging):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with ordinary kriging interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_kriging:       the parameters of the input for "kriging"
    Output:
        returns the value of the area
 
    """  
    
    
    print("File written to", j_kriging['output-file'])
