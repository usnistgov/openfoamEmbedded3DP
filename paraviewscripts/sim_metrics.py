#!/usr/bin/env python
'''Functions for collecting data from simulations'''
# RG

from interfacePlots import *
import cv2 as cv
from pylab import *
import numpy as np

CONV = 4.221/260 # converts pixels to millimeters (bath height/image length)

#------------- repeat-use functions ------------------

def threshold(im:np.array) -> Union[float, np.array]:
    '''tresholds images'''
    im = removeNozzle(im)
    imgray = cv.cvtColor(im, cv.COLOR_BGR2GRAY)
    ret, thresh = cv.threshold(imgray, 10, 255, cv.THRESH_BINARY) # first number 1-170 seems acceptable
    return ret, thresh

def removeNozzle(im:np.array) -> np.array:
    '''if the nozzle is in the image, remove it'''
    ret, thresh = cv.threshold(im, 220, 255, cv.THRESH_BINARY)
    im[thresh==255] = 0
    return im

#------------- fusion metrics ------------------

def getPerimeter(im:np.array, n:int) -> float:
    '''measure perimeter of object
    n is the number of lines
    componentMask is a binarized image of just one segment'''
    im = im[240:500, 480:740] # crop image
    cm = im.copy()
    ret, thresh = threshold(cm)
    contours, hierarchy = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
    perimeters = [-1]*len(contours)
    for i,cnt in enumerate(contours):
        perimeters[i] = cv.arcLength(cnt,True)*CONV
        cv.drawContours(cm, cnt, -1, (110, 245, 209), 6)
    if sum(perimeters)==0:
        return -1
    # imshow(cm)
    # cv.imshow('Perimeter', cm)
    ideal = ((n-1)/2+1)*np.pi*0.603 # ideal perimeter of n lines side-to-side
    return sum(perimeters)/ideal-1
    # return perimeters/ideal
    
    
def getRoughness(im:np.array) -> List[float]:
    '''measure roughness as perimeter of object / perimeter of convex hull. 
    componentMask is a binarized image of just one segment'''
    im = im[240:500, 480:740]
    cm = im.copy()
    ret, thresh = threshold(cm)
    contours, hierarchy = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key=lambda x: cv.contourArea(x), reverse=True) # select the largest contour
    cnt = contours[0]
    perimeter = cv.arcLength(cnt,True)
    if perimeter==0:
        return -1
    hull = cv.convexHull(cnt)
    hullperimeter = cv.arcLength(hull,True)
    roughness = perimeter/hullperimeter-1  # how much extra perimeter there is compared to the convex hull

    # annotate image
    cv.drawContours(cm, [hull], -1, (110, 245, 209), 6)
    cv.drawContours(cm, cnt, -1, (186, 6, 162), 6)
    # imshow(cm)
    return roughness


def asymmetry(im:np.array, vert:bool) -> Union[float, float]:
    '''measures how symmetric the filament(s) is/are. 
    vert True means look at vertical symmetry
    vert False means look at horizontal symmetry'''
    im = im[240:500, 480:740]
    if not vert:
        im = cv.rotate(im, cv.ROTATE_90_CLOCKWISE)
    cm = im.copy()
    ret, thresh = threshold(cm)

    lens = []
    for row in range(cm.shape[0]): # go through each row and record length between first and last ink pixel in row
        r = thresh[row]
        left = np.argmax(r)
        r = np.flip(r)
        right = cm.shape[1] - np.argmax(r)
        if left != 0 and right != 0:
            lens.append(right-left)
            if row % 10 == 0:
                cv.line(cm, (left, row), (right, row), (0,0,255), 1)
    lenslength = int(len(lens)/2)
    lenst = lens[0:lenslength]
    lensb = lens[-1:lenslength:-1]
    diff = [abs(1-t/b) for t,b in zip(lenst,lensb)] # indicator of how symmetrical the cross section is (ideal 0)
    avgdiff = sum(diff)/len(diff)
    sd = np.std(diff) # how much this symmetry varies (ideal 0)
    if not vert:
        cm = cv.rotate(cm, cv.ROTATE_90_COUNTERCLOCKWISE)
    # imshow(cm)
    return avgdiff, sd


def fused(im) -> bool:
    '''check if the filament is fused
    ignores tiny specks which may be blur artifacts'''
    im = im[240:500, 480:740]
    cm = im.copy()
    ret, thresh = threshold(cm)
    contours, hierarchy = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
    j = 0
    for i in contours:
        if cv.arcLength(i,True) > 30: # arbitrary small perimeter
            j+=1
    f = j==1
    return f

#----------- no extrusion metrics --------------------

def contoursNoLeak(thresh:np.array) -> Union[List[tuple], np.array]: # verify this works for the leaks
    '''gets contours of an image with no extrusion, ignoring leakage from nozzle'''
    c, h = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
    for i in c:
        if cv.arcLength(i,True) > 30:
            del i
    return c, h


def dispY(im1, im2, lines) -> int:
    '''Measures horizontal displacement of filament due to nozzle
    im1 is the cross section after the nozzle
    im2 is the cross section before the nozzle
    lines is the number of lines in the simulation'''
    im1 = im1[240:500, 480:740] # crop image
    cm1 = im1.copy()
    ret1, thresh1 = threshold(cm1)
    contours1, hierarchy1 = contoursNoLeak(thresh1)
    contours1 = sorted(contours1, key=lambda x: cv.contourArea(x), reverse=True) # select the largest contour
    contour1 = contours1[0]
    M1 = cv.moments(contour1)
    Y = int(M1["m10"] / M1["m00"])*CONV
    
    im2 = im2[240:500, 480:740] # crop image
    cm2 = im2.copy()
    ret2, thresh2 = threshold(cm2)
    contours2, hierarchy2 = cv.findContours(thresh2, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
    contours2 = sorted(contours2, key=lambda x: cv.contourArea(x), reverse=True) # select the largest contour
    contour2 = contours2[0]
    M2 = cv.moments(contour2)
    
    if lines==1: 
        idealY = int(M2["m10"]/M2["m00"])*CONV
        return idealY-Y
    
    if lines==2:
        adjY = int(M2["m10"]/M2["m00"])*CONV
        nozY = 0
        idealY = (adjY+nozY)/2
        return idealY-Y
    
    raise Exception('Valid number of lines is 1 or 2')
                    

def dispZ(im1, im2, lines) -> int:
    '''Measures horizontal displacement of filament due to nozzle'''
    im1 = im1[240:500, 480:740] # crop image
    cm1 = im1.copy()
    ret1, thresh1 = threshold(cm1)
    contours1, hierarchy1 = contoursNoLeak(thresh1)
    contours1 = sorted(contours1, key=lambda x: cv.contourArea(x), reverse=True) # select the largest contour
    contour1 = contours1[0]
    M1 = cv.moments(contour1)
    Z = int(M1["m01"] / M1["m00"])*CONV
    
    im2 = im2[240:500, 480:740] # crop image
    cm2 = im2.copy()
    ret2, thresh2 = threshold(cm2)
    contours2, hierarchy2 = cv.findContours(thresh2, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
    contours2 = sorted(contours2, key=lambda x: cv.contourArea(x), reverse=True) # select the largest contour
    contour2 = contours2[0]
    M2 = cv.moments(contour2)
                    
    if lines==1:
        idealZ = int(M2["m10"]/M2["m00"])*CONV
        return idealZ-Z
    
    if lines==2:
        adjZ = int(M2["m10"]/M2["m00"])*CONV
        nozZ = 0.30015
        idealZ = (adjY+nozY)/2
        return idealZ-Z
    
    raise Exception('Valid number of lines is 1 or 2')