#!/usr/bin/env python
"""
Created on Thu Nov 19 09:26:59 2015

@author: zhang
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import h5py
import glob


def get_vars_from_hdf5(filepath):
    
    f = h5py.File(filepath, 'r')
    
    d={}
    
    d['rho'] = f.get('rho')[...]
    d['u'] = f.get('u')[...]
    d['v'] = f.get('v')[...]
    d['w'] = f.get('w')[...]
    d['p'] = f.get('p')[...]
    d['T'] = f.get('T')[...]
    
    d['x'] = f.get('x')[...]
    d['y'] = f.get('y')[...]
    d['z'] = f.get('z')[...]

    d['time'] = np.float(f.get('time')[...])
    d['time_step'] = np.int(f.get('time_step')[...])
    
    d['Ni'] = f.get('Ni')[...]
    d['Nj'] = f.get('Nj')[...]
    d['Nk'] = f.get('Nk')[...]
    
    f.close()
    
    return d


filenames = glob.glob('./*h5')
filenames.sort()

Nt = len(filenames)

for count,file in enumerate(filenames):
    
    filepath = './' + file
    data = get_vars_from_hdf5(filepath)
    print ('Acquired data from: '+file)


    ################################################
    #Pick the data on the center surface (center on z direction) for contour plot
    #figureName = 'ContourPlots_'+str(file)
    X,Y = np.meshgrid(data['x'],data['y'])
    Nk = np.int(data['Nk'])
    Nk_center = Nk/2
    U_center = np.transpose(np.array(data['u'])[:,:,Nk_center])
    V_center = np.transpose(np.array(data['v'])[:,:,Nk_center])
    T_center = np.transpose(np.array(data['T'])[:,:,Nk_center])
    P_center = np.transpose(np.array(data['p'])[:,:,Nk_center])


    ### Plotting
    
    fontSize     = 30
    gcaFontSize  = 20
    figwidth  = 18;
    figheight = 8;
    
    lineWidth = 1
    iskip_plt = 1
    plot_markers = False
    markerSize = 0.5
    
   
    xmin = np.min(X)
    xmax = np.max(X)
    ymin = np.min(Y)
    ymax = np.max(Y)
    #pmin = np.min(P_center)
    #pmax = np.max(P_center)
    pmin = 2.0
    pmax = 4.0
    p_levels = np.linspace(pmin,pmax,50)
    #Tmin = np.min(T_center)
    #Tmax = np.max(T_center)
    Tmin = 1.0
    Tmax = 5.0
    T_levels = np.linspace(Tmin,Tmax,50)
    #Umin = np.min(U_center)
    #Umax = np.max(U_center)
    Umin = 0.0
    Umax = 6.0
    U_levels = np.linspace(Umin,Umax,50)
    #Vmin = np.min(V_center)
    #Vmax = np.max(V_center)
    Vmin = -1.0
    Vmax = 1.0
    V_levels = np.linspace(Vmin,Vmax,50)
    
    figureName = 'Contour_P_T_' + str(count).zfill(5)
    fig = plt.figure(0, figsize=(figwidth,figheight))
    
    ax0   = fig.add_subplot(2,1,1)
    ax1   = fig.add_subplot(2,1,2)
    
    ax = ax0
    plt.axes(ax)
    plt.contourf(X,Y,P_center,p_levels,interpolation='nearest', cmap=plt.cm.jet, extend="both")
    plt.colorbar()
    plt.title("Pressure")

    ax = ax1
    plt.axes(ax)
    plt.contourf(X,Y,T_center,T_levels,interpolation='nearest', cmap=plt.cm.jet, extend="both")
    plt.colorbar()
    plt.title("Temperature")

    for ax in [ax0,ax1]:
        plt.axes(ax)
        ax.set_xlabel(r"$x$",fontsize=fontSize)
        ax.set_ylabel(r"$y$",fontsize=fontSize) # (\% of Atm)
    
        plt.setp(ax.get_xticklabels(),fontsize=gcaFontSize)
        plt.setp(ax.get_yticklabels(),fontsize=gcaFontSize)
    
        plt.axis('equal')
        plt.axis('tight')
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
    
    print ("Saving figure: " + figureName)
    plt.tight_layout()
    plt.savefig(figureName,dpi=100)
    plt.close()

    figureName = 'Contour_U_V_' + str(count).zfill(5)
    fig = plt.figure(0, figsize=(figwidth,figheight))
    
    ax0 = fig.add_subplot(2,1,1)
    ax1 = fig.add_subplot(2,1,2)

    ax = ax0
    plt.axes(ax)
    plt.contourf(X,Y,U_center,U_levels,interpolation='nearest', cmap=plt.cm.jet, extend="both")
    plt.colorbar()
    plt.title("U")

    ax = ax1
    plt.axes(ax)
    plt.contourf(X,Y,V_center,V_levels,interpolation='nearest', cmap=plt.cm.jet, extend="both")
    plt.colorbar()
    plt.title("V")


    for ax in [ax0,ax1]:
        plt.axes(ax)
        ax.set_xlabel(r"$x$",fontsize=fontSize)
        ax.set_ylabel(r"$y$",fontsize=fontSize) # (\% of Atm)
    
        plt.setp(ax.get_xticklabels(),fontsize=gcaFontSize)
        plt.setp(ax.get_yticklabels(),fontsize=gcaFontSize)
        
        plt.axis('equal')
        plt.axis('tight')
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        
    print ("Saving figure: " + figureName)
    plt.tight_layout()
    plt.savefig(figureName,dpi=100)
    plt.close()

########################################









