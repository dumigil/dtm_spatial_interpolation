import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import scipy.spatial
import sys
import math
import csv
import random
import json 
import time

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
    print(gridsize)
    x = []
    y = []
    z = []
    sample_size = 0
    for point in list_pts_3d:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
        sample_size += 1 
    kd = scipy.spatial.KDTree(list_pts_3d)
    d, i = kd.query([1,20,75], k=1)
    print(d, i)
    """
    # data coordinates and values

    # target grid to interpolate to
    xi = yi = np.arange(min(x),max(y)+2,10)
    xi,yi = np.meshgrid(xi,yi)

    # set mask
    #mask = (xi > 0.5) & (xi < 0.6) & (yi > 0.5) & (yi < 0.6)

    # interpolate
    zi = griddata((x,y),z,(xi,yi),method='nearest')

    # mask out the field
    #zi[mask] = np.nan

    # plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    plt.contourf(xi,yi,zi,np.arange(min(z),max(z),1))
    #plt.plot(x,y,'k.')
    plt.xlabel('xi',fontsize=16)
    plt.ylabel('yi',fontsize=16)
    #plt.savefig('interpolated.png',dpi=100)
    plt.show()
    """

if __name__ == '__main__':
    main()