# -*- coding: utf-8 -*-
'''
Some simple functions

'''

import subprocess as sp
from os import path, system
import numpy as np

def printlog(str, logfile):
    print(str)
    logfile.write(str)
    logfile.write("\n")

def area(x1, y1, x2, y2, x3, y3):
    #  calculate the triangula area that given by coordiinate prameters
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1)
                + x3 * (y1 - y2)) / 2.0)


def isInside(x1, y1, x2, y2, x3, y3, x, y):
    #  check whether point P(x, y) lies inside the triangle formed by  A(x1, y1), B(x2, y2) and C(x3, y3)

    # area of base triangle
    a = area(x1, y1, x2, y2, x3, y3)
    # area of 1st triangle
    a1 = area(x, y, x2, y2, x3, y3)
    # area of 2nd triangle
    a2 = area(x1, y1, x, y, x3, y3)
    # area of 3rd triangle
    a3 = area(x1, y1, x2, y2, x, y)
    # Check if sum of A1, A2 and A3 is same as A
    if round(a, 9) == round(a1 + a2 + a3, 9):
        return True
    else:
        return False


def triangCoversPoint(edgeLenghts, centroidCoords, vertexCoords, xi, yi):
    """ determines triangle that point (xi,yi)  is inside
      returns the;
        indice or triange
        centroid coordinates of triangle as 2-d array
        vertex coordinates of tiangle as 2x3 matix
        """

    tIndice = 0
    tCentroid = 0
    tVertex = 0
    maxElen = np.amax(edgeLenghts)

    for i in xrange(0, len(centroidCoords), 1):
        if xi > centroidCoords[i][0] + maxElen:
            continue
        elif xi < centroidCoords[i][0] - maxElen:
            continue
        elif yi > centroidCoords[i][1] + maxElen:
            continue
        elif yi < centroidCoords[i][1] - maxElen:
            continue
        elif xi > vertexCoords[3 * i][0] and \
                xi > vertexCoords[3 * i + 1][0] and \
                xi > vertexCoords[3 * i + 2][0]:
            continue
        elif xi < vertexCoords[3 * i][0] and \
                xi < vertexCoords[3 * i + 1][0] and \
                xi < vertexCoords[3 * i + 2][0]:
            continue
        elif yi > vertexCoords[3 * i][1] and \
                yi > vertexCoords[3 * i + 1][1] and \
                yi > vertexCoords[3 * i + 2][1]:
            continue
        elif yi < vertexCoords[3 * i][1] and \
                yi < vertexCoords[3 * i + 1][1] and \
                yi < vertexCoords[3 * i + 2][1]:
            continue
        elif isInside(vertexCoords[3 * i][0], vertexCoords[3 * i][1], vertexCoords[3 * i + 1][0],
                      vertexCoords[3 * i + 1][1], vertexCoords[3 * i + 2][0], vertexCoords[3 * i + 2][1], xi, yi):
            tIndice = i
            tCentroid = centroidCoords[i]
            tVertex = [vertexCoords[3 * i], vertexCoords[3 * i + 1], vertexCoords[3 * i + 2]]
            print "centorid:" + str(centroidCoords[i]) + "--vertex1:" + str(vertexCoords[3 * i]) + "--vertex2:" \
                  + str(vertexCoords[3 * i + 1]) + "vertex3:" + str(vertexCoords[3 * i + 2])
            break
    return tIndice, tCentroid, tVertex
