from pylab import *
import sys, shutil, os
from mpl_toolkits.mplot3d import Axes3D

def read_2D(title,plottitle):
    infile = open(title)
    a = infile.readline().split()
    for i in range(len(a)):
        a[i] = float(a[i])
    
    u = zeros([len(a),len(a)])
    u[0,:] = a
    j = 1
    for line in infile:
        sline = line.split()
        for i in range(len(sline)):
            a[i] = float(sline[i])
        u[j,:] = a
        j += 1

    fig = figure()
    ax = Axes3D(fig)
    X = linspace(0,1,len(a))
    Y = linspace(0,1,len(a))
    X,Y = meshgrid(X,Y)
    ax.plot_surface(X,Y,u, rstride=1, cstride=1,cmap='cool')
    
    #ax.contourf(X, Y, u, zdir='z', offset=-10000, cmap=cm.hot)
    
    xlabel("x"); ylabel("y"); ax.set_zlabel("u(x,y)")
    ax.set_title(plottitle)
    show()

def read_1D_MC(title):

    infile = open(title)
    X = []
    for line in infile:
        if float(line) >= 0:
            X.append(float(line))
    infile.close()

    Y = [1.]
    j = 0
    for i in range(1,len(X)-1):
        if X[i+1] == X[i]:
            Y[j] += 1.
        else:
            j += 1
            Y.append(1.)

    print Y[0]
    Y = array(Y)
    Y = Y/Y[0]
    x = linspace(0,1,len(Y))

    plot(x,Y,'-')
    show()


read_1D_MC('build-main-Desktop-Debug/RandomWalk_a.txt')
read_2D("build-main-Desktop-Debug/Explicit.txt","Explicit sheme for 2+1D diffusion equation")
read_2D("build-main-Desktop-Debug/Implicit.txt","Implicit sheme for 2+1D diffusion equation")

