import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
plt.ion()
 
# Define a function to make the ellipses
def ellipse(ra,rb,ang,x0,y0,Nb=100):
    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang
    co,si=np.cos(an),np.sin(an)
    the=linspace(0,2*np.pi,Nb)
    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
    return X,Y


def make2Dplot(fig,Data_BEAMX,Data_BEAMY,Nbinx,Nbiny):

    # Define the x and y data 
    # For example just using random numbers
	
    x = Data_BEAMX
    y = Data_BEAMY
 
    # Set up default x and y limits
    xlims = [min(x),max(x)]
    ylims = [min(y),max(y)]
  
    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
 
    # Set up the geometry of the three plots
    rect_beam = [left, bottom, width, height]  # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25]  # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height]  # dimensions of y-histogram
 
    # Set up the size of the figure
    #fig = plt.figure(1, figsize=(9.5,9))
	
    # Make the three plots
    axBeam = fig.add_axes(rect_beam) # beam plot
    mean_DataX=str(round(Data_BEAMX.mean(),3))
    std_DataX=str(round(Data_BEAMX.std(),3))
    mean_DataY =str(round(Data_BEAMY.mean(),3))
    std_DataY =str(round(Data_BEAMY.std(),3))

    axHistx = fig.add_axes(rect_histx) # x histogram
    axHistx.set_ylabel("Counts")
    axHistx.set_title('Mean : '+mean_DataX + ' std : '+std_DataX)
    axHistx.grid(True)

    axHisty = fig.add_axes(rect_histy) # y histogram
    axHisty.set_xlabel("Counts")
    axHisty.set_title('Mean : '+mean_DataY + ' std : '+std_DataY,rotation=270,x=1.08,y=0.75)
    axHisty.grid(True)
	
    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
	
    # Find the min/max of the data
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(ylims)
	
    # Make the 'main' beam plot
    # Define the number of bins
    nxbins = Nbinx
    nybins = Nbiny
	
    xbins = linspace(start = xmin, stop = xmax, num = nxbins)
    ybins = linspace(start = ymin, stop = ymax, num = nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0
    aspectratio = 1.0*(xmax - 0)/(1.0*ymax - 0)
	
    H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins))
    X = xcenter
    Y = ycenter
    Z = H
 
    # Plot the beam data
    cax = (axBeam.imshow(H, extent=[xmin,xmax,ymin,ymax],
    interpolation='nearest', origin='lower',aspect='auto'))

    print('xmin : '+str(xmin) +' xmax : ' +str(xmax) + ' ymin: ' + str(ymin) + ' ymax: '+str(ymax))
    # Plot the beam plot contours
    contourcolor = 'white'
    xcenter = np.mean(x)
    ycenter = np.mean(y)
    ra = np.std(x)
    rb = np.std(y)
    ang = 0 ##To change for rotated ellipse : call georges.phys
	
    X,Y=ellipse(ra,rb,ang,xcenter,ycenter)
    axBeam.plot(X,Y,"k:",ms=1,linewidth=2.0)
    axBeam.annotate('$1\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
				textcoords='offset points', horizontalalignment='right',
				verticalalignment='bottom',fontsize=25)
 
    # X,Y=ellipse(2*ra,2*rb,ang,xcenter,ycenter)
    # axBeam.plot(X,Y,"k:",color = contourcolor,ms=1,linewidth=2.0)
    # axBeam.annotate('$2\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
    #             textcoords='offset points',horizontalalignment='right',
    #             verticalalignment='bottom',fontsize=25, color = contourcolor)
 
    X,Y=ellipse(3*ra,3*rb,ang,xcenter,ycenter)
    axBeam.plot(X,Y,"k:",color = contourcolor, ms=1,linewidth=2.0)
    axBeam.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
                       textcoords='offset points',horizontalalignment='right',
                       verticalalignment='bottom',fontsize=25, color = contourcolor)

    # Plot the axes labels
    axBeam.set_xlabel(Data_BEAMX.name, fontsize=25)
    axBeam.set_ylabel(Data_BEAMY.name, fontsize=25)

    #Make the tickmarks pretty
    ticklabels = axBeam.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family('serif')
 
    ticklabels = axBeam.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family('serif')
 
    #Set up the plot limits
    axBeam.set_xlim(xlims)
    axBeam.set_ylim(ylims)

    #Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax-xmin)/nxbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)/nybins)
	
    #Plot the histograms
    axHistx.hist(x, bins=xbins, color = 'blue', histtype='step', normed=True)
    axHisty.hist(y, bins=ybins, orientation='horizontal', color = 'red', histtype='step', normed=True)
    #Set up the histogram limits
    axHistx.set_xlim( min(x), max(x) )
    axHisty.set_ylim( min(y), max(y) )
 
    #Make the tickmarks pretty
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')
 
    #Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')
 
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))

    #Show the plot
    #ax.draw()

