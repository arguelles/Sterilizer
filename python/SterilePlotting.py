import numpy
import pylab

def Draw2DAbsHistogram(data, bins,title="",SaveAs="",clim=[0,5]):
    pylab.rc('font', family='serif', size=27)
    pylab.rc('legend', fontsize = 20)
    pylab.rc('axes', labelsize = 27)
    pylab.rc('axes', titlesize = 27)
    pylab.rcParams['xtick.major.pad']='8'
    pylab.rcParams['ytick.major.pad']='8'

    yticks=[3,4]
    ylabels=[r"$10^{3}$",r"$10^{4}$"]
    
    ebins=bins[0]
    zbins=bins[1]
    zbincenters=(zbins[1:]+zbins[:-1])/2.
    ebincenters=(ebins[1:]+ebins[:-1])/2.
    ev,cv=numpy.meshgrid(numpy.log10(ebincenters),zbincenters)
    ev_flat=ev.flatten()
    cv_flat=cv.flatten()
    data_flat=data.flatten()

    fig=pylab.figure(figsize=(9,8))
    pylab.hist2d(cv_flat,ev_flat,weights=numpy.log10(data_flat), bins=(zbins,numpy.log10(ebins)),cmap=pylab.get_cmap('Blues'))
    for i in range(0,len(data_flat)):
        if(not numpy.isnan(data_flat[i]) and not numpy.isinf(data_flat[i])):
            pylab.annotate("{:.0f}".format(data_flat[i]),xy=(cv_flat[i],ev_flat[i]),size=12,color='black',rotation='vertical')
    pylab.clim(clim[0],clim[1])

    if(title!=""):
        pylab.title(title, y=1.02)    
    pylab.yticks(yticks,ylabels)
    pylab.xlabel(r"cos$\theta^{reco}_{\mu,z}}$")
    pylab.ylabel(r"$E^{reco}_{\mu,proxy}/ GeV$")

    pylab.tight_layout()
    if(SaveAs!=""):
        print("saving as ",SaveAs)
        pylab.savefig(SaveAs,dpi=300,bbox_inches='tight')
    pylab.show()


def Draw2DRatioHistogram(data1, data2, bins,title="",SaveAs="",clim=[-20,20]):
    pylab.rc('font', family='serif', size=27)
    pylab.rc('legend', fontsize = 20)
    pylab.rc('axes', labelsize = 27)
    pylab.rc('axes', titlesize = 27)
    pylab.rcParams['xtick.major.pad']='8'
    pylab.rcParams['ytick.major.pad']='8'

    yticks=[3,4]
    ylabels=[r"$10^{3}$",r"$10^{4}$"]

    ebins=bins[0]
    zbins=bins[1]
    zbincenters=(zbins[1:]+zbins[:-1])/2.
    ebincenters=(ebins[1:]+ebins[:-1])/2.
    ev,cv=numpy.meshgrid(numpy.log10(ebincenters),zbincenters)
    ev_flat=ev.flatten()
    cv_flat=cv.flatten()
    data_flat=(((data1-data2)/data1)*100).flatten()

    fig=pylab.figure(figsize=(9,8))
    pylab.hist2d(cv_flat,ev_flat,weights=data_flat, bins=(zbins,numpy.log10(ebins)),cmap=pylab.get_cmap('RdBu'))
    for i in range(0,len(data_flat)):
        if(not numpy.isnan(data_flat[i]) and not numpy.isinf(data_flat[i])):
            pylab.annotate("{:.1f}".format(data_flat[i]),xy=(cv_flat[i]-0.01,ev_flat[i]),size=12,color='black',rotation='vertical')
    pylab.clim(clim[0],clim[1])
    pylab.tick_params(width=1, length = 8)
    pylab.colorbar(format='%1.1f',label = r"Deviation (%)")

    if(title!=""):
        pylab.title(title, y=1.02)
    pylab.yticks(yticks,ylabels)
    pylab.xlabel(r"cos$\theta^{reco}_{\mu,z}}$")
    pylab.ylabel(r"$E^{reco}_{\mu,proxy}/ GeV$")

    pylab.tight_layout(pad=0.2,h_pad=0.2,w_pad=0.2)
    if(SaveAs!=""):
        print("saving as ",SaveAs)
        pylab.savefig(SaveAs,dpi=300,bbox_inches='tight')
    pylab.show()



def Draw2DPullPlot(data, expectation, bins,title="",SaveAs="",clim=[-5,5]):
    pylab.rc('font', family='serif', size=27)
    pylab.rc('legend', fontsize = 20)
    pylab.rc('axes', labelsize = 27)
    pylab.rc('axes', titlesize = 27)
    pylab.rcParams['xtick.major.pad']='8'
    pylab.rcParams['ytick.major.pad']='8'

    yticks=[3,4]
    ylabels=[r"$10^{3}$",r"$10^{4}$"]

    ebins=bins[0]
    zbins=bins[1]
    zbincenters=(zbins[1:]+zbins[:-1])/2.
    ebincenters=(ebins[1:]+ebins[:-1])/2.
    ev,cv=numpy.meshgrid(numpy.log10(ebincenters),zbincenters)
    ev_flat=ev.flatten()
    cv_flat=cv.flatten()
    data_flat=((data-expectation)/numpy.sqrt(expectation)).flatten()

    fig=pylab.figure(figsize=(9,8))
    pylab.hist2d(cv_flat,ev_flat,weights=data_flat, bins=(zbins,numpy.log10(ebins)),cmap=pylab.get_cmap('RdBu'))
    for i in range(0,len(data_flat)):
        if(not numpy.isnan(data_flat[i]) and not numpy.isinf(data_flat[i])):
            pylab.annotate("{:.1f}".format(data_flat[i]),xy=(cv_flat[i]-0.01,ev_flat[i]),size=12,color='black',rotation='vertical')
    pylab.clim(clim[0],clim[1])
    pylab.tick_params(width=1, length = 8)
    pylab.colorbar(format='%1.1f',label = r"Pull per bin ($\sigma$)")

    if(title!=""):
        pylab.title(title, y=1.02)
    pylab.yticks(yticks,ylabels)
    pylab.xlabel(r"cos$\theta^{reco}_{\mu,z}}$")
    pylab.ylabel(r"$E^{reco}_{\mu,proxy}/ GeV$")

    pylab.tight_layout(pad=0.2,h_pad=0.2,w_pad=0.2)
    if(SaveAs!=""):
        print("saving as ",SaveAs)
        pylab.savefig(SaveAs,dpi=300,bbox_inches='tight')
    pylab.show()

