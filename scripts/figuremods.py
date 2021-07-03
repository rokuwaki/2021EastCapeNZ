import warnings
warnings.filterwarnings('ignore')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.image as mpimg
from mpl_toolkits.basemap import Basemap
from xml.dom import minidom
import numpy as np
import os
import subprocess
from cmcrameri import cm
import shapefile
import pandas as pd
from netCDF4 import Dataset
from scipy.interpolate import griddata
from pyrocko.plot import beachball
from pyrocko import moment_tensor as pmt
import obspy
from obspy import read, read_inventory
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84
import utils
figsize = (5.6, 5.6)

# Data and model directory
datarootdir = '../materials/'


def fig1(figname):
    j = 13
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    fig = plt.figure(figsize=figsize)    

    # basemap and bathymetry
    axpxloc, axpyloc, axpwidth = 0.1, 0.1, 0.6
    tickintx, tickinty, tickformat = 1, 1, 0
    lonmin, lonmax, latmin, latmax = lonminb, lonmaxb, latminb, latmaxb = 178.4, 181, -39.1, -36    
    m, ax, ax2 = utils.drawbasemap(fig,lonmin,lonmax,latmin,latmax,axpxloc,axpyloc,axpwidth,tickintx,tickinty,tickformat)
    utils.drawbathy(datarootdir+'work/bathymetry/bathy_cut.nc', fig, ax, m)
    
    # trench with indentation
    xmldoc = minidom.parse(datarootdir+'HikurangiTrench.kml')
    coordinates = xmldoc.getElementsByTagName('coordinates')[0].firstChild.data.replace('\n', '').replace('\t', '').split()
    coordinates = np.array([coordinates[i].split(',') for i in range(len(coordinates))]).astype(float)
    for n in range(len(coordinates)):
        if coordinates[n,0] < 0: coordinates[n,0] += 360
    x, y = m(coordinates[:,0], coordinates[:,1])
    ax.plot(x, y, color=cm.lapaz(255), lw=0.75, linestyle='--', zorder=1, alpha=0.5)

    # beachballs
    utils.drawbeachball(datarootdir+'work/SPUD_QUAKEML_bundle.xml',datarootdir,lonmin,lonmax,latmin,latmax,ax,m)

    # epicenter
    model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
    elat, elon = model_para.lat[0], model_para.lon[0]
    x, y=m(elon, elat)
    sc=ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=100, 
                  path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])

    # draw some labels and plate motion vectors
    utils.drawlabels(elat, elon, ax, m)
    
    # regional map
    axp = ax.get_position()
    axpxloc, axpyloc, axpwidth = axp.x1+0.005, axp.y0, 0.18
    tickintx, tickinty, tickformat = 10, 10, 0
    latmin, latmax, lonmin, lonmax = -50, -15, 165, 187
    m, ax, ax2 = utils.drawbasemap(fig,lonmin,lonmax,latmin,latmax,axpxloc,axpyloc,axpwidth,tickintx,tickinty,tickformat,1,'l',0.2)
    ax2.set_xticks([]); ax2.set_yticks([])
    
    # epicenter
    x, y=m(elon, elat)
    sc=ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=101, 
                  path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])
    
    # plate boundaries
    src = shapefile.Reader(datarootdir+'work/tectonicplates/PB2002_boundaries.shp')
    for tmp in src.shapeRecords():
        x, y = [i[0] for i in tmp.shape.points[:]], [i[1] for i in tmp.shape.points[:]]
        for n in range(len(x)):
            if x[n] < 0: x[n] += 360
        x, y = m(x, y)
        ax.plot(x, y, color='k', lw=0.5, linestyle='--', zorder=0)

    # left-map region and plate names
    tmplon = [lonminb, lonmaxb, lonmaxb, lonminb, lonminb]
    tmplat = [latminb, latminb, latmaxb, latmaxb, latminb]
    x, y = m(tmplon, tmplat)
    ax.plot(x, y, color='k', lw=0.75, solid_capstyle='round', zorder=100)
    tmplon, tmplat, texts = [180,170,179], [-47,-30,-34], ['PA','AU','KE']
    for lon, lat, text in zip(tmplon, tmplat, texts):    
        x, y = m(lon, lat)
        ax.text(x, y, text, path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])
    
    plt.savefig(figname, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.show()
    
    
def fig2(figname, ptflag, j, rcmtflag=1, synflag=0):
    fig = plt.figure(figsize=figsize)
    ax = fig.gca(projection='3d', alpha=1)
    azim, elev = 195, 25
    
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
    elat, elon = model_para.lat[0], model_para.lon[0]
    data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/FFM_DCpreferred.txt', skiprows=1)
    lon, lat, depth, slip = data[:,2], data[:,3], data[:,4], data[:,1]
    for n in range(len(lon)):
        if lon[n] < 0: lon[n] += 360

    # 3D frame
    ax, axlabel = utils.drawbaseframe3D(model_para, fig, ax, azim, elev)    
    utils.drawslip3D(datarootdir, model_para, model[j], fig, ax, axlabel, lon, lat, depth, slip)
    
    # 2D basemap
    axpxloc, axpyloc, axpwidth = 1.02, 0.18, 0.4
    tickintx, tickinty, tickformat = 0.2, 0.2, 1
    lonmargin, latmargin = 0.05, 0.05
    lonmin=np.min(lon)-lonmargin; lonmax=np.max(lon)+lonmargin+0.06; latmin=np.min(lat)-latmargin; latmax=np.max(lat)+latmargin
    for axes in range(2):
        if axes == 0:
            axpxloc = axpxloc
        else:
            axp = ax.get_position()
            axpxloc = axp.x1+0.01
        m, ax, ax2 = utils.drawbasemap(fig,lonmin,lonmax,latmin,latmax,axpxloc,axpyloc,axpwidth,tickintx,tickinty,tickformat)
        ax.set_facecolor('whitesmoke')
        
        # plate boundaries
        src = shapefile.Reader(datarootdir+'work/tectonicplates/PB2002_boundaries.shp')
        for tmp in src.shapeRecords():
            x, y = [i[0] for i in tmp.shape.points[:]], [i[1] for i in tmp.shape.points[:]]
            for n in range(len(x)):
                if x[n] < 0: x[n] += 360
            x, y = m(x, y)
            ax.plot(x, y, color='k', lw=0.5, linestyle='--', zorder=0)
        x, y = m(179.875, -37.4)
        ax.text(x, y, 'Trench', ha='left', va='center', fontsize=8, color='k', rotation=65,
                path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])
        
        #slab2df = pd.read_csv(datarootdir+'work/ker_slab2_clp_02.24.18.csv', header=None, names=['lon', 'lat'])
        #x, y = m(slab2df['lon'], slab2df['lat'])
        #ax.plot(x, y, color='r', lw=0.75, linestyle='--')
        
        
        # epicenter
        x, y=m(elon, elat)
        sc=ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=101, 
                      path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])
                
        # beach balls and P axes
        utils.drawPaxis(fig, ax, ax2, m, axes, datarootdir, model[j], model_para)
        
        if rcmtflag == 1:
            # R-CMTs
            utils.drawRCMT(fig, ax, m, datarootdir, axes, model_para)
            
    if synflag == 1:
        if j == 13:
            fig.text(0.16, 0.81, 'Input model')
        elif j == 21:
            fig.text(0.16, 0.81, 'Reproduced model')

    
    plt.savefig(figname, pad_inches=0.1, dpi=300, transparent=False, bbox_inches=mpl.transforms.Bbox([[0.8, 0.6], [11, 4.65]]))
    plt.show()
    
    
    
def fig3(figname):
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    for j in np.arange(13, 14, 1):
        model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')

        data = np.loadtxt(datarootdir+'model_'+str(model[j])+'/tw_mec.dat')
        snap_t, str_total, dip_total, rake_total = data[:,0], data[:,9], data[:,11], data[:,13]
        m1, m2, m3, m4, m5, m6 = data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8]

        col_names = ['l','k','tw','x','y','slip',
                    'mrr', 'mtt', 'mff', 'mrt', 'mrf', 'mtf',
                    'strike1','strike2','dip1','dip2','rake1','rake2']
        df = pd.read_table(datarootdir+'model_'+str(model[j])+'/snap.dat', header=None, delim_whitespace=True, names=col_names)
        dep = [ model_para.depth[0] - np.sin( np.deg2rad(model_para.dip[0]) )*((df['y'][i])) for i in range(len(df)) ]
        df['dep'] = dep

        cmap = cm.bilbao
        fig = plt.figure(figsize=figsize)
        xmin, xmax = np.min(df['x'])-model_para.xx[0]/2, np.max(df['x'])+model_para.xx[0]/2
        ymin, ymax = 0, np.max(df['dep'])+model_para.yy[0]/2

        snapnum = 0
        for tw in np.arange(1, 7, 1):

            axwidth = 0.15
            axheigt=(ymax-ymin) / (xmax-xmin) * axwidth
            if snapnum == 0:
                axxloc = 0.1
                axyloc = 0.1
            else:
                axxloc = axpmain.x1+0.01
                axyloc = axpmain.y0

            ax = fig.add_axes([axxloc, axyloc, axwidth, axheigt])
            if snapnum == 0:
                axp0 = ax.get_position()

            axpmain = ax.get_position()
            snapint = 5
            fig.text(axpmain.x0+0.005, axpmain.y1+0.002, str(snapnum*snapint)+'–'+str((snapnum+1)*snapint)+' s', va='bottom',
                    path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])


            tmpdf = df[(df['tw'] == tw)].reset_index(drop=True)
            ax.scatter(0, model_para.depth[0], s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=10,
                      path_effects=[path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])

            tmpx, tmpy = tmpdf['x'], tmpdf['dep']
            xi=np.linspace(np.min(tmpx)-model_para.xx[0]/2, np.max(tmpx)+model_para.xx[0]/2, 100)
            yi=np.linspace(np.min(tmpy)-model_para.yy[0]/2, np.max(tmpy)+model_para.yy[0]/2, 100)

            X, Y=np.meshgrid(xi, yi)
            zi=griddata((tmpx, tmpy), tmpdf['slip'], (X, Y),'linear')
            interval=np.linspace(0, np.max(df['slip']), 11)
            sc=ax.contourf(X, Y, zi, interval, vmin=0, vmax=np.max(df['slip']), cmap=cmap)
            interval = [0.063 , 0.126]
            ax.contour(X, Y, zi, interval, colors='k', linewidths=0.75)

            ax.set_xlim(xmin, xmax)
            ax.set_xticklabels([])
            ax.set_ylim(ymin, ymax)
            ax.set_yticks(np.arange(0, 120, 20))
            ax.invert_yaxis()

            if snapnum == 0:
                ax.set_ylabel('Depth (km)')
            else:
                ax.set_yticklabels([])

            ax.axhline(7, linestyle='--', lw=0.75, color='k')

            axp = ax.get_position()

            snapnum += 1

        #############################################################################################
        #############################################################################################
        #############################################################################################
        axp = ax.get_position()
        axheight = axp.height*2+0.015
        axwidth = axheight * (xmax-xmin) / (ymax-ymin)
        ax = fig.add_axes([axp.x1+0.07, axp.y1-axheight, axwidth, axheight])
        axpsummary = ax.get_position()
        cmap = cm.nuuk

        snapnum = 0
        centroidloc = []
        for tw in np.arange(1, 7, 1):
            tmpdf = df[(df['tw'] == tw)].reset_index(drop=True)
            ax.scatter(0, model_para.depth[0], s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=10,
                      path_effects=[path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])

            tmpx, tmpy = tmpdf['x'], tmpdf['dep']
            xi=np.linspace(np.min(tmpx)-model_para.xx[0]/2, np.max(tmpx)+model_para.xx[0]/2, 100)
            yi=np.linspace(np.min(tmpy)-model_para.yy[0]/2, np.max(tmpy)+model_para.yy[0]/2, 100)

            X, Y=np.meshgrid(xi, yi)
            zi=griddata((tmpx, tmpy), tmpdf['slip'], (X, Y),'linear')
            interval = [0.063, 0.126]
            ax.contourf(X, Y, zi, interval, colors=mpl.colors.rgb2hex(cmap(snapnum/5)), alpha=np.max(tmpdf['slip'])/np.max(df['slip'])*0.75, zorder=int(snapnum/5*10))
            ax.contour(X, Y, zi, interval, colors='k', linewidths=0.75, zorder=int(snapnum/5*10))
            if snapnum <= 3:
                tmpx, tmpy = 22-3, model_para.depth[0]-55
                ax.plot([tmpx, tmpx+17], [tmpy+snapnum*4.5, tmpy+snapnum*4.5], color=mpl.colors.rgb2hex(cmap(snapnum/5)), clip_on=False,
                        alpha=np.max(tmpdf['slip'])/np.max(df['slip'])*0.75, lw=7)
                ax.text(tmpx, tmpy+0.35+snapnum*4.5, 'E'+str(snapnum+1)+': '+str(snapnum*snapint)+'–'+str((snapnum+1)*snapint)+' s', ha='left', va='center', color='k', fontsize=6,
                       path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])
                centroidloc.append([tmpdf.loc[tmpdf['slip'].idxmax(), 'x'], tmpdf.loc[tmpdf['slip'].idxmax(), 'y']])

                meclocs, mos, focmecsRCMT, labels = utils.loadmecs(datarootdir)
                if snapnum == 1 or snapnum == 2:
                    if snapnum == 1:
                        k = 4
                        cmtdepth = 52
                        timeshift = 6.8
                        color_t = mpl.colors.rgb2hex(cmap(snapnum/5))
                    elif snapnum == 2:
                        k = 5
                        cmtdepth = 12
                        timeshift = 14.6
                        color_t = mpl.colors.rgb2hex(cmap(snapnum/5))
                    tmp = beachball.plot_beachball_mpl(focmecsRCMT[k], ax, size=12, position=(xmin-5-0.5, cmtdepth),
                                                       beachball_type='dc', edgecolor='C7', color_t=color_t,
                                                       color_p='none', linewidth=0.75, alpha=0.75, zorder=100).set_clip_on(False)
                    #ax.axhline(cmtdepth, color='C7', linestyle='--', lw=0.5, zorder=-1)
                    ax.text(xmin-1-0.5, cmtdepth-4, 'R-CMT\n'+str(timeshift)+' s', ha='right', va='bottom', fontsize=6, color='k',
                           path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])
                    ax.text(xmin, cmtdepth, '$\leftarrow$', ha='left', va='center', fontsize=6, color='k',
                           path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])


            snapnum += 1

        ax.axhline(7, linestyle='--', lw=0.75, color='k')

        ax.set_xlabel('Strike (km)')
        ax.set_ylabel('Depth (km)')

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.invert_yaxis()
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.set_label_position('right')


        ax.text(3, model_para.depth[0]-17, 'E1', ha='center', va='center', color='k', fontsize=8,
               path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()], zorder=50)

        ax.text(23, model_para.depth[0]+18, 'E2', ha='center', va='center', color='k', fontsize=8,
               path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()], zorder=50)

        ax.text(-8, model_para.depth[0]-46, 'E3', ha='center', va='center', color='k', fontsize=8,
               path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()], zorder=50)

        ax.text(23, model_para.depth[0]-8, 'E3', ha='center', va='center', color='k', fontsize=8,
               path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()], zorder=50)

        ax.text(-5, model_para.depth[0]+10, 'E4', ha='center', va='center', color='k', fontsize=8,
               path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()], zorder=50)

        axp = ax.get_position()
        fig.text(axp.x1-0.05-0.015, axp.y1-0.15, 'Slip rate\n≥'+str(interval[0])+' m/s', ha='center', va='top', fontsize=6)
        fig.text(axp.x0+0.01, axp.y0+0.005, r'$\leftarrow$North', ha='left', va='bottom', fontsize=8)
        fig.text(axp.x1-0.01, axp.y0+0.005, r'South$\rightarrow$', ha='right', va='bottom', fontsize=8)

        fig.text(axp.x1-0.005-0.015, axp.y1-0.032, 'Plate interface', ha='right', va='center', color='k', fontsize=6,
               path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()], zorder=50)

        ax.set_facecolor('whitesmoke')

        axp = ax.get_position()
        cmap = cm.bilbao
        cax=fig.add_axes([axp.x1, axp.y0-0.2, 0.01, 0.1])
        norm=mpl.colors.Normalize(vmin=0, vmax=np.max(df['slip']))
        cb=mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, label='Slip rate (m/s)',
                                     ticks=np.linspace(0, np.max(df['slip']), 3), format='%.2f')

        axmrf=fig.add_axes([axp.x0, axp.y0-0.2, 0.18, 0.1])
        data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/st_'+str(model[j])+'.dat')
        mrft, mrfamp=data[:,0], data[:,1]
        axmrf.plot(mrft, mrfamp, color='k', lw=0.5)
        axmrf.fill(mrft, mrfamp, color='C7')
        axmrf.set_ylabel('Moment rate\n'+r'($\times 10^{18}$ Nm/s)')
        axmrf.set_xlabel('Time (s)')
        axmrf.set_ylim([0, max(mrfamp)+max(mrfamp)*0.1])
        axmrf.set_xticks(np.arange(0, 50, 10))
        axmrf.set_xlim(-2, 35)
        model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
        moment = model_para.moment[0]
        mw = model_para.mw[0]
        #print(model_para)
        axpmrf = axmrf.get_position()
        #text=fig.text(axpmrf.x1-0.005, axpmrf.y1-0.005,
        #        str('{:.2e}'.format(moment))+' Nm\n($M_{\mathrm{W}}$ '+str('{:.1f}'.format(mw))+')', size=6, va='top', ha='right')
        axmrf.yaxis.tick_right()
        axmrf.yaxis.set_label_position("right")


        #############################################################################################
        #############################################################################################
        #############################################################################################

        cmap = cm.bilbao
        snapnum = 0
        for tw in np.arange(1, 7, 1):
            axwidth = 0.15
            axheigt=(ymax-ymin) / (xmax-xmin) * axwidth
            if snapnum == 0:
                axxloc = axp0.x0
                axyloc = axp0.y0-axp0.height-0.015
            else:
                axxloc = axpmain.x1+0.01
                axyloc = axpmain.y0

            ax = fig.add_axes([axxloc, axyloc, axwidth, axheigt])
            axpmain = ax.get_position()

            tmpdf = df[(df['tw'] == tw)].reset_index(drop=True)
            ax.scatter(0, model_para.depth[0], s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=10,
                      path_effects=[path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])

            azi_list = []
            for i in range(len(tmpdf)):
                focmecs=[tmpdf['mrr'][i], tmpdf['mtt'][i], tmpdf['mff'][i], tmpdf['mrt'][i], tmpdf['mrf'][i], tmpdf['mtf'][i]]
                focmecs = utils.convertUSEtoNED(focmecs)
                focmecsDC=[tmpdf['strike1'][i], tmpdf['dip1'][i], tmpdf['rake1'][i]]
                mt = pmt.MomentTensor(
                    strike=tmpdf['strike1'][i],
                    dip=tmpdf['dip1'][i],
                    rake=tmpdf['rake1'][i])
                pnt = utils.get_plunge_azimuth(mt.m6_up_south_east())
                t_azi = pnt[5] # T-axis azimuth
                p_azi = pnt[1] # p-axis azimuth
                azi_list.append(p_azi)

                try:
                    tmp = beachball.plot_beachball_mpl(focmecs, ax, size=7*tmpdf['slip'][i]/np.max(df['slip']), position=(tmpdf['x'][i], tmpdf['dep'][i]),
                                             beachball_type='dc', edgecolor='C7', color_t=cmap(tmpdf['slip'][i]/np.max(df['slip'])),
                                             color_p='w', linewidth=0.2,
                                                       alpha=tmpdf['slip'][i]/np.max(df['slip']), zorder=(int(2+tmpdf['slip'][i]/np.max(df['slip'])*10)))
                except:
                    tmp = beachball.plot_beachball_mpl(focmecsDC, ax, size=7*tmpdf['slip'][i]/np.max(df['slip']), position=(tmpdf['x'][i], tmpdf['dep'][i]),
                                             beachball_type='dc', edgecolor='C7', color_t=cmap(tmpdf['slip'][i]/np.max(df['slip'])),
                                             color_p='w', linewidth=0.2,
                                                       alpha=tmpdf['slip'][i]/np.max(df['slip']), zorder=(int(2+tmpdf['slip'][i]/np.max(df['slip'])*10)))

                x = np.cos(np.deg2rad(90-p_azi)) * tmpdf['slip'][i]/np.max(df['slip']) * 15
                y = np.sin(np.deg2rad(90-p_azi)) * tmpdf['slip'][i]/np.max(df['slip']) * 15 * -1 # !! [-1] is required because y-axis is later inverted
                #x = np.cos(np.deg2rad(90-t_azi)) * tmpdf['slip'][i]/np.max(df['slip']) * 10
                #y = np.sin(np.deg2rad(90-t_azi)) * tmpdf['slip'][i]/np.max(df['slip']) * 10 * 1
                for pm in [-1, 1]:
                    ax.plot([tmpdf['x'][i], tmpdf['x'][i]+pm*x], [tmpdf['dep'][i], tmpdf['dep'][i]+pm*y],
                            color='k', lw=0.5, alpha=tmpdf['slip'][i]/np.max(df['slip']), zorder=1)

            dx, dy = model_para.xx[0], model_para.yy[0]
            for i in range(len(tmpdf)):
                tmpx, tmpy = tmpdf['x'][i], tmpdf['dep'][i]
                tmpxlist = [tmpx-dx/2, tmpx+dx/2, tmpx+dx/2, tmpx-dx/2, tmpx-dx/2]
                tmpylist = [tmpy-dy/2, tmpy-dy/2, tmpy+dy/2, tmpy+dy/2, tmpy-dy/2]
                ax.fill(tmpxlist, tmpylist, facecolor=cmap(tmpdf['slip'][i]/np.max(df['slip'])), edgecolor='none', zorder=0)

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_yticks(np.arange(0, 120, 20))
            ax.invert_yaxis()

            if snapnum == 0:
                ax.set_xlabel('Strike (km)')
                ax.set_ylabel('Depth (km)')

            else:
                ax.set_yticklabels([])

            ax.axhline(7, linestyle='--', lw=0.75, color='k')

            axp = ax.get_position()


            #############################
            #############################
            baxbase = fig.add_axes([axp.x0, axpmrf.y0, axp.width, axpmrf.height])
            baxbase.set_xticks([]); baxbase.set_yticks([])
            bax0 = fig.add_axes([axp.x0, axp.y0-axp.width*0.49-0.112, axp.width*0.49, axp.width*0.49])
            scale = np.max(tmpdf['slip'])/np.max(df['slip'])
            baxp = bax0.get_position()
            bax1 = fig.add_axes([baxp.x1+(baxp.width*(1-scale)/2)-0.005,
                                 baxp.y0+(baxp.width*(1-scale)/2),
                                 baxp.width*scale, baxp.height*scale], projection='polar')
            bax0.set_xticks([]); bax0.set_yticks([])
            bax0.set_axis_off()
            if snapnum == 5:
                baxp = bax1.get_position()
                fig.text(baxp.x0+baxp.width/2, baxp.y1+0.002, 'P-axis\nazimuth\nN', fontsize=6, ha='center', va='bottom',
                         path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])
                fig.text(baxp.x0+baxp.width/2, baxp.y0-0.005, 'S', fontsize=6, ha='center', va='top')
                fig.text(baxp.x1+0.004, baxp.y0+baxp.height/2, 'E', fontsize=6, ha='left', va='center')
                fig.text(baxp.x0-0.004, baxp.y0+baxp.height/2, 'W', fontsize=6, ha='right', va='center')

            focmecs=[m1[tw-1],m2[tw-1],m3[tw-1],m4[tw-1],m5[tw-1],m6[tw-1]] # total moment tensor for each snapshot
            focmecs = utils.convertUSEtoNED(focmecs)
            x, y = 0.5, 0.5
            tmp = beachball.plot_beachball_mpl(focmecs, bax0, size=20*np.max(tmpdf['slip'])/np.max(df['slip']), position=(x, y),
                                     beachball_type='deviatoric', edgecolor='none', color_t=cmap(np.max(tmpdf['slip'])/np.max(df['slip'])),
                                     color_p='w', linewidth=0, alpha=1, zorder=2)
            tmp = beachball.plot_beachball_mpl(focmecs, bax0, size=20*np.max(tmpdf['slip'])/np.max(df['slip']), position=(x, y),
                                     beachball_type='dc', edgecolor='C7', color_t='none',
                                     color_p='none', linewidth=0.5, alpha=1, zorder=2)

            tmp = np.array(azi_list)
            bin_edges = np.arange(-5, 366, 10)
            number_of_strikes, bin_edges = np.histogram(tmp, bin_edges, weights=tmpdf['slip']/np.max(df['slip']))

            number_of_strikes[0] += number_of_strikes[-1]
            half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
            two_halves = np.concatenate([half, half])

            bax1.bar(np.deg2rad(np.arange(0, 360, 10)), two_halves,
                   width=np.deg2rad(10), bottom=0.0, color='C7', edgecolor='k', linewidth=0.1)
            x = [ (bin_edges[i-1]+bin_edges[i])/2 for i in range(1, len(bin_edges)) ]
            bax1.set_theta_zero_location('N')
            bax1.set_theta_direction(-1) # direction in which theta increases: clockwise
            bax1.set_thetagrids(np.arange(0, 360, 45), labels=[])
            #ax.set_rgrids(np.arange(1, two_halves.max() + 1, 2000), angle=0, labels=[])
            bax1.set_rgrids(np.arange(0, two_halves.max(), 2000), angle=0, labels=[])
            bax1.set_xticklabels([])
            bax1.grid(linewidth=0.35)
            bax1.set_axisbelow(True)

            snapnum += 1

        plt.savefig(figname, bbox_inches="tight", pad_inches=0.1, dpi=300)
        plt.show()

        
def figSynFFMmerge(figname, inputimage, outputimage):
    fig = plt.figure(figsize=(10, 10))

    image = mpimg.imread(inputimage)
    ax = fig.add_axes([0, 0, image.T.shape[1]*1e-3, image.T.shape[2]*1e-3])
    im = ax.imshow(image)
    ax.axis('off')

    ax = fig.add_axes([0, 0-image.T.shape[2]*1e-3, image.T.shape[1]*1e-3, image.T.shape[2]*1e-3])
    image = mpimg.imread(outputimage)
    im = ax.imshow(image)

    ax.axis('off')

    plt.savefig(figname, pad_inches=0.1, dpi=50, transparent=False,bbox_inches='tight')
    plt.close(fig)
    plt.show()

    
def fit(figname):
    def aziequi(ax, stalist, distancetextsize):
        data = np.loadtxt(stalist, usecols=(5, 4))
        d, a = (data[:, 0], 90-data[:,1])
        x, y=(d*np.cos(a*np.pi/180.0), d*np.sin(a*np.pi/180.0))
        sc=ax.scatter(x, y, s=15, marker='^', edgecolor='k', facecolor='ivory', alpha=1, lw=0.75, zorder=10)
        #sc=ax.scatter(x, y, s=8, marker='^', edgecolor='k', facecolor='none', alpha=1, lw=0.5, zorder=10)
        stalist=np.loadtxt(stalist, usecols=(0), dtype=str)
        #print(stalist)
        #for (i, staname) in enumerate(stalist):
        #    text = ax.text(x[i], y[i], staname, size=6, va='center')
        #    text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='w', alpha=0.85), path_effects.Normal()])
        ax.scatter(0, 0, s=100, marker='*', edgecolor='k', facecolor='none', lw=0.75)
        theta=np.linspace(0, 360, 360)
        for i in [30, 90]:
            x, y=(i*np.cos(theta*np.pi/180.0), i*np.sin(theta*np.pi/180.0))
            ax.plot(x, y, color='C7', zorder=0, solid_capstyle='round', lw=1, linestyle='--')
            x, y=(i*np.cos(-90*np.pi/180.0), i*np.sin(-90*np.pi/180.0))
            text = ax.text(x, y, str(i)+'$\degree$', size=distancetextsize, va='center', ha='center')
            text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])

        x, y=(100*np.cos(theta*np.pi/180.0), 100*np.sin(theta*np.pi/180.0))
        ax.plot(x, y, color='k', solid_capstyle='round', lw=1)
        ax.fill(x, y, edgecolor='none', facecolor='w', zorder=0)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_yticks([])


    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    for j in np.arange(13, 14, 1):
        modelid = model[j]
        alpha=1
        axw=0.25
        axh=0.075

        fig=plt.figure(figsize=figsize)

        stalist=np.loadtxt(datarootdir+'model_'+str(modelid)+'/station_'+str(modelid)+'.list', usecols=(0), dtype=str)
        stadata=np.loadtxt(datarootdir+'model_'+str(modelid)+'/station_'+str(modelid)+'.list', usecols=(4, 5))
        azi, dis=(stadata[:,0], stadata[:,1])

        datalen = np.loadtxt(datarootdir+'model_'+str(modelid)+'/'+str(modelid)+'.station.abic', skiprows=1, usecols=5)
        maxdatalen = max(datalen)
        for i in np.arange(0, len(stalist), 1):
        #for i in np.arange(0, 12, 1):
            num = 15
            mod = i // num
            if i == 0:
                ax=fig.add_axes([0.1, 0.1, axw, axh])
                axp0=ax.get_position()
                #fig.text(axp.x0, axp.y1+0.03, str(pwd)+', '+modelid, va='bottom', ha='left')
            elif i == num*(mod):
                ax=fig.add_axes([axp0.x1+axw*(mod-1)+0.08*(mod), axp0.y0, axw, axh])
            elif i > num*(mod) and i < num*(mod+1):
                ax=fig.add_axes([axp.x0, axp.y0-axh, axw, axh])
            axp=ax.get_position()

            if i == num*(mod+1)-1 or i == len(stalist)-1:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(True)
                ax.spines['left'].set_visible(False)
                ax.set_xlabel('Time (s)')
                ax.set_yticklabels([])
                ax.set_yticks([])
                axpbottom = ax.get_position()
            else:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.set_xticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_yticks([])

            sta=stalist[i]
            obsdata=datarootdir+'model_'+str(modelid)+'/obssyn/obs_'+sta+'.txt'
            syndata=datarootdir+'model_'+str(modelid)+'/obssyn/syn_'+sta+'.txt'

            data=np.loadtxt(obsdata, usecols=(0, 1))
            t, amp=(data[:,0]-5, data[:,1])
            ax.plot(t, amp, color='k', lw=1, alpha=alpha)
            data=np.loadtxt(syndata, usecols=(0, 1))
            t, amp=(data[:,0]-5, data[:,1])
            ax.plot(t, amp, color='r', lw=1, alpha=alpha)

            fig.text(axp.x0, axp.y1-(axp.height/2), sta+'\nAz.='+str(azi[i])+'$\degree$\nDel.='+str(dis[i])+'$\degree$', \
                                va='center', ha='right', size=6, color='k', clip_on=False, alpha=alpha)
            ax.set_ylim([-1.2, 1.2])
            ax.set_xlim([-5, 75])
            ax.set_xticks(np.arange(0, 80, 20))

        ax=fig.add_axes([axp.x0, axpbottom.y0-axp.width-0.1, axp.width, axp.width])
        axp=ax.get_position()
        sc = aziequi(ax, datarootdir+'model_'+str(modelid)+'/station_'+str(modelid)+'.list', 10)

        plt.savefig(figname, bbox_inches='tight', pad_inches=0.1, dpi=200)
        plt.show()

        
def dcpopulation(figname):
    
    models = (np.loadtxt(datarootdir+'modellist.txt', usecols=0, dtype=int)).astype(str)
    for model in models[13:14]:
        model_para = utils.load_fort40(os.path.join(datarootdir, 'model_'+model, 'fort.40'))
        fig=plt.figure(figsize=figsize)
        axx0, axy0 = 0.1, 0.1
        axw = 0.35

        aspect = (model_para.xx[0] * model_para.mn[0])/(model_para.yy[0] * model_para.nn[0])
        axh = axw / aspect

        ax = fig.add_axes([axx0, axy0, axw, axh])
        axp = ax.get_position()

        data = np.loadtxt(os.path.join(datarootdir, 'model_'+model, 'FFM_DCall.txt'), skiprows=1)
        lon,lat,dep,slip = data[:,1],data[:,2],data[:,3],data[:,4]
        strike0,dip0,rake0,strike1,dip1,rake1 = data[:,5],data[:,6],data[:,7],data[:,8],data[:,9],data[:,10]
        xloc,yloc = data[:,11],data[:,12]

        data=np.loadtxt(os.path.join(datarootdir, 'model_'+model, 'FFM_MT.txt'))
        m1,m2,m3,m4,m5,m6 = data[:,4],data[:,5],data[:,6],data[:,7],data[:,8],data[:,9]

        strike, dip, rake = [],[],[]
        for i in range(len(strike0)):
            tmpstr, tmpdip, tmprake = utils.selectplane(model_para.strike[0],model_para.dip[0],strike0[i],dip0[i],rake0[i],strike1[i],dip1[i],rake1[i])
            strike.append(tmpstr)
            dip.append(tmpdip)
            rake.append(tmprake)

        dx, dy = model_para.xx[0], model_para.yy[0]
        xmin, xmax = np.min(xloc)-dx/2, np.max(xloc)+dx/2
        ymin, ymax = np.min(yloc)-dy/2, np.max(yloc)+dy/2
        cmap = cm.bilbao

        for i in np.arange(len(slip)):
            tmpx, tmpy = xloc[i], yloc[i]
            tmpxlist = [tmpx-dx/2, tmpx+dx/2, tmpx+dx/2, tmpx-dx/2, tmpx-dx/2]
            tmpylist = [tmpy-dy/2, tmpy-dy/2, tmpy+dy/2, tmpy+dy/2, tmpy-dy/2]
            ax.fill(tmpxlist, tmpylist, facecolor=cmap(slip[i]/max(slip)), edgecolor='none', zorder=0)


        arrowflag = 'beach'
        x, y = xloc, yloc
        for i in range(len(strike)):
            if arrowflag == 'rake':
                a=np.deg2rad(rake[i])
                length=slip[i] / max(slip) * dx
                x1=np.cos(a)*length
                y1=np.sin(a)*length
                ax.plot([x[i], x1+x[i]], [y[i], y1+y[i]], color='C2', lw=1, solid_capstyle='round', clip_on=False)
                x2, y2 =np.cos(a-np.deg2rad(150))*length*0.2, np.sin(a-np.deg2rad(150))*length*0.2
                ax.plot([x1+x[i], x1+x2+x[i]], [y1+y[i], y1+y2+y[i]], color='C2', lw=1, solid_capstyle='round', clip_on=False)
                x2, y2 =np.cos(a+np.deg2rad(150))*length*0.2, np.sin(a+np.deg2rad(150))*length*0.2
                ax.plot([x1+x[i], x1+x2+x[i]], [y1+y[i], y1+y2+y[i]], color='C2', lw=1, solid_capstyle='round', clip_on=False)

            elif arrowflag == 'strike':
                a=np.deg2rad(-strike[i]+90)
                length=slip[i] / max(slip) * dx
                x1=np.cos(a)*length
                y1=np.sin(a)*length
                ax.plot([x[i]-x1, x1+x[i]], [y[i]-y1, y1+y[i]], color='C2', lw=1, solid_capstyle='round', clip_on=False)            

            elif arrowflag == 'dip':
                a=np.deg2rad(dip[i]-180)
                length=slip[i] / max(slip) * dy
                x1=np.cos(a)*length
                y1=np.sin(a)*length
                ax.plot([x[i]-x1, x1+x[i]], [y[i]-y1, y1+y[i]], color='C2', lw=1, solid_capstyle='round', clip_on=False)

            elif arrowflag == 'beach':
                focmec = [m1[i],m2[i],m3[i],m4[i],m5[i],m6[i]]
                focmec = utils.convertUSEtoNED(focmec)
                tmp = beachball.plot_beachball_mpl(focmec, ax, size=dx*1.5, position=(x[i], y[i]),
                                         beachball_type='deviatoric', edgecolor='none', color_t=cmap(slip[i]/max(slip)),
                                         color_p='w', linewidth=0.5, alpha=1, zorder=1, view='top')
                tmp = beachball.plot_beachball_mpl(focmec, ax, size=dx*1.5, position=(x[i], y[i]),
                                         beachball_type='dc', edgecolor='k', color_t='none',
                                         color_p='none', linewidth=0.5, alpha=1, zorder=1, view='top')

        ax.scatter(0, 0, marker='*', s=500, facecolor='none', edgecolor='k', 
                   path_effects = [path_effects.Stroke(linewidth=2, foreground='w'), path_effects.Normal()])

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel('Strike (km)')
        ax.set_ylabel('Dip (km)')
        if model_para.dip[0] > 0:
            ax2 = ax.twinx()
            depmin = model_para.depth[0]-np.sin(np.deg2rad(model_para.dip[0]))*ymin
            depmax = model_para.depth[0]-np.sin(np.deg2rad(model_para.dip[0]))*ymax
            ax2.set_ylim(depmin, depmax)
            ax2.set_ylabel('Depth (km)')
            cax = fig.add_axes([axp.x1+0.07, axp.y0, 0.01, axp.height/4])
        else:
            cax = fig.add_axes([axp.x1+0.005, axp.y0, 0.01, axp.height/4])

        norm = mpl.colors.Normalize(vmin=0, vmax=np.max(slip))
        mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, label='Slip (m)', ticks=np.linspace(0,np.max(slip),3), format='%.2f')


        #plt.savefig('slip_'+arrowflag+'_'+model+'.png', bbox_inches='tight', pad_inches=0.1, dpi=300)        
        #plt.show()    
    
    
    dep0, dep1 = 0, 50
    ratio_dc, ratio_dc_shallow, ratio_dc_deep = [], [], []
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    for j in np.arange(13, 14, 1):
        model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')

        data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/FFM_DCpreferred.txt', skiprows=1)
        lon, lat, slip, strike, dip, rake = data[:,2], data[:,3],data[:,1], data[:,7], data[:,8],data[:,9]
        depth = data[:,4]

        data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/FFM_MT.txt')
        m1,m2,m3,m4,m5,m6 = data[:,4],data[:,5],data[:,6],data[:,7],data[:,8],data[:,9]

        for i in range(len(depth)):

            mt = pmt.MomentTensor.from_values(
                [m1[i], m2[i], m3[i], m4[i], m5[i], m6[i]])
            focmec = mt.m6()
            decon = mt.standard_decomposition()
            ratio_dc.append(decon[1][1]*100)

            if depth[i] <= 50:
                ratio_dc_shallow.append(decon[1][1]*100)
            elif depth[i] > 50:
                ratio_dc_deep.append(decon[1][1]*100)

    #fig = plt.figure(figsize=figsize)
    axp = ax.get_position()
    ax = fig.add_axes([axp.x1+0.22, axp.y0, axp.height, axp.height])

    bins=np.arange(0, 110, 10)
    ax.hist(ratio_dc_shallow, label='Shallow layers\n'+'($\leq50$ km depth)', bins=bins, density=True, facecolor=cmap(0.5/np.max(slip)), edgecolor='none', lw=0.5, alpha=0.55)
    ax.hist(ratio_dc_deep, label='Deep layers\n'+'($>50$ km depth)', bins=bins, density=True, facecolor=cmap(1./np.max(slip)), edgecolor='none', lw=0.5, alpha=0.55)
    ax.hist(ratio_dc_shallow, bins=bins, density=True, facecolor='none', edgecolor='w', lw=0.5)
    ax.hist(ratio_dc_deep, bins=bins, density=True, facecolor='none', edgecolor='w', lw=0.5)

    ax.set_xlim(-5, 105)

    ax.set_yticks([])
    ax.set_xlabel('DC ratio (%)')
    ax.set_ylabel('Normalised counts')
    ax.legend(loc='upper left', fontsize=8, handlelength=1, facecolor='none')

    plt.savefig(figname, pad_inches=0.1, dpi=150, bbox_inches="tight")
    plt.show()

    
def geometryFFM(figname, j, geomflag):
    fig = plt.figure(figsize=(5.4, 5.4))
    ax = fig.gca(projection='3d', alpha=1)
    azim, elev = 195, 25
    
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
    elat, elon = model_para.lat[0], model_para.lon[0]
    data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/FFM_DCpreferred.txt', skiprows=1)
    lon, lat, depth, slip = data[:,2], data[:,3], data[:,4], data[:,1]
    for n in range(len(lon)):
        if lon[n] < 0: lon[n] += 360
    lonmin=model_para.lon[0]-1*0.85
    lonmax=model_para.lon[0]+1*0.85
    latmin=model_para.lat[0]-1*0.5
    latmax=model_para.lat[0]+0.75*0.5

    # 3D frame
    ax, axlabel = utils.drawbaseframe3D(model_para, fig, ax, azim, elev)    
    
    # Model geometries
    vertslist = []
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0, comments='#')
    model_paras = pd.read_csv(datarootdir+'FFMmodelPara.csv')
    
    if geomflag == 'vertical':
        models = np.arange(10, 20, 1)
    else:
        models = [0,1,2,3,4,5,6,7,8,9,13,20]
    
    for j in models:
        #model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
        model_para = model_paras[model_paras['modelid'] == model[j]].reset_index()

        if geomflag == 'vertical':
            if j == 13:
                verts = utils.plotModelEdge3D(model_para, ax, 1)
                pc = Poly3DCollection(verts)
                pc.set_alpha(0.1)
                pc.set_facecolor('C5')
                ax.add_collection3d(pc)
                
        else:
            verts = utils.plotModelEdge3D(model_para, ax, 1)
            pc = Poly3DCollection(verts)
            pc.set_alpha(0.1)
            pc.set_zorder(1001)
            if j == 13:
                pc.set_facecolor('C5')
            else:
                pc.set_facecolor('C7')
            ax.add_collection3d(pc)
            

        if model_para.depth[0] == 72:
            ax.scatter(model_para.lon[0], model_para.lat[0], model_para.depth[0], marker='*', s=100, edgecolor='k', facecolor='C5', zorder=1000)
        else:
            ax.scatter(model_para.lon[0], model_para.lat[0], model_para.depth[0], marker='*', s=100, edgecolor='k', facecolor='w', zorder=1000)
    
    
    ax.plot([lonmin, lonmax], [latmin, latmin], [0, 0], color='k', lw=0.75)
    ax.plot([lonmin, lonmin], [latmin, latmin], [0, 115], color='k', lw=0.75)
    ax.plot([lonmax, lonmax], [latmin, latmin], [0, 115], color='k', lw=0.75)
    ax.plot([lonmin, lonmin], [latmax, latmax], [0, 115], color='k', lw=0.75)
    ax.plot([lonmax, lonmax], [latmax, latmax], [0, 115], color='k', lw=0.75)
    
    
    plt.savefig(figname, pad_inches=0, dpi=300, transparent=False, bbox_inches=mpl.transforms.Bbox([[0.9, 0.7], [4.7, 4.3]]))
    #plt.close(fig)
    plt.show()
    
    
def varianceFFM(figname):
    deplist = []
    varlist = []
    vrlist = []

    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    model_paras = pd.read_csv(datarootdir+'FFMmodelPara.csv')
    
    for j in np.arange(0, 21, 1):
        #model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
        model_para = model_paras[model_paras['modelid'] == model[j]].reset_index()
        
        
        if model_para.dip.values == 0.0 and model_para.strike.values == 200.0:
            #print(model_para.depth.values, model_para.variance.values)
            deplist.append(model_para.depth.values[0])
            varlist.append(model_para.variance.values[0])
            vrlist.append(model_para.vr.values[0])

    tmp = np.argsort(np.array(deplist))

    deplist = [ deplist[i] for i in tmp ]
    varlist = [ varlist[i] for i in tmp ]

    fig=plt.figure(figsize=figsize)
    ax = fig.add_axes([0.1, 0.1, 0.5, 0.5])
    ax.plot(deplist, varlist, 'o-', color='k', markersize=6, label='Horizontal models')
    ax.set_xlabel('Hypo depth (km)')
    ax.set_ylabel('Variance')

    axp = ax.get_position()
    fig.text(axp.x0+0.01, axp.y1-0.01, 'Model plane geometry', ha='left', va='top', fontsize=10)

    #model_para = utils.load_fort40(datarootdir+'model_'+str(model[20])+'/fort.40')
    model_para = model_paras[model_paras['modelid'] == model[20]].reset_index()
    ax.scatter(model_para.depth.values, model_para.variance.values, facecolor='k', s=30, marker='D', edgecolor='w', lw=0.75,
               label='Vertical model '+r'($\phi_{\rm{m}}=110\degree$, $\delta_{\rm{m}}=90\degree$)', zorder=100)

    #model_para = utils.load_fort40(datarootdir+'model_'+str(model[13])+'/fort.40')
    model_para = model_paras[model_paras['modelid'] == model[13]].reset_index()
    ax.scatter(model_para.depth.values, model_para.variance.values, facecolor='C5', s=40, lw=1, edgecolor='k',
               marker='o', label='Optimum model '+r'($\phi_{\rm{m}}=200\degree$, $\delta_{\rm{m}}=90\degree$)')

    ax.legend(loc='lower left', fontsize=10, facecolor='none', handletextpad=0.1, borderpad=0.1, labelspacing=0.1)

    ax.set_ylim(0.2, 0.4)
    ax.set_yticks(np.arange(0.2, 0.5, 0.05))
    ax.set_xlim(0, 120)
    #############################################################
    #############################################################
    #############################################################

    deplist = []
    varlist = []
    vrlist = []

    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    for j in np.arange(0, 21, 1):
        #model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
        model_para = model_paras[model_paras['modelid'] == model[j]].reset_index()
        if model_para.dip.values == 90.0 and model_para.strike.values == 200.0:
            #print(model_para.depth.values, model_para.variance.values)

            deplist.append(model_para.depth.values[0])
            varlist.append(model_para.variance.values[0])
            vrlist.append(model_para.vr.values[0])

    tmp = np.argsort(np.array(deplist))

    deplist = [ deplist[i] for i in tmp ]
    varlist = [ varlist[i] for i in tmp ]

    axp = ax.get_position()
    ax = fig.add_axes([axp.x0, axp.y0-axp.height-0.1, axp.width, axp.height])
    ax.plot(deplist, varlist, 'o-', color='k', markersize=6, label='Vertical model '+r'($\phi_{\rm{m}}=200\degree$, $\delta_{\rm{m}}=90\degree$)')

    #model_para = utils.load_fort40(datarootdir+'model_'+str(model[13])+'/fort.40')
    model_para = model_paras[model_paras['modelid'] == model[13]].reset_index()
    ax.scatter(model_para.depth.values, model_para.variance.values, facecolor='C5', s=40, lw=1, edgecolor='k', zorder=10,
               marker='o', label='Optimum model '+r'($\phi_{\rm{m}}=200\degree$, $\delta_{\rm{m}}=90\degree$)')

    axp = ax.get_position()
    fig.text(axp.x0+0.01, axp.y1-0.01, 'Nucleation depth', ha='left', va='top', fontsize=10)
    ax.set_yticks(np.arange(0.25, 0.28, 0.01))
    ax.set_ylim(0.25, 0.29)
    ax.set_xlim(0, 120)
    ax.legend(loc='lower left', fontsize=10, facecolor='none', handletextpad=0.1, borderpad=0.1, labelspacing=0.1)
    ax.set_xlabel('Hypo depth (km)')
    ax.set_ylabel('Variance')

    plt.savefig(figname, bbox_inches="tight", pad_inches=0, dpi=300, transparent=False)
    #plt.close(fig)
    plt.show()
    
    
def sensitivityMerge(figname, modelgeometry0, modelgeometry1, variance):
    fig = plt.figure(figsize=(10, 10))

    image = mpimg.imread(modelgeometry0)
    ax = fig.add_axes([0, 0.05, image.T.shape[1]*1e-3, image.T.shape[2]*1e-3])
    im = ax.imshow(image)
    ax.axis('off')

    ax = fig.add_axes([0, 0.05-image.T.shape[2]*1e-3, image.T.shape[1]*1e-3, image.T.shape[2]*1e-3])
    image = mpimg.imread(modelgeometry1)
    im = ax.imshow(image)
    ax.axis('off')

    image = mpimg.imread(variance)
    ax = fig.add_axes([1.2, -1, image.T.shape[1]*1e-3, image.T.shape[2]*1e-3])
    im = ax.imshow(image)
    ax.axis('off')

    plt.savefig(figname, pad_inches=0.1, dpi=50, transparent=False,bbox_inches="tight")
    plt.close(fig)
    plt.show()

    
def slab2(figname, ptflag):
    fig = plt.figure(figsize=figsize)
    cmap = cm.bilbao
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    for j in np.arange(13, 14, 1):
        model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')

        axpxloc, axpyloc, axpwidth = 1.02, 0.18, 0.4
        tickintx, tickinty, tickformat = 1, 1, 0
        lonmargin, latmargin = 1, 1.1

        data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/FFM_DCpreferred.txt', skiprows=1)
        lon, lat, slip, depth = data[:,2], data[:,3],data[:,1], data[:,4]

        data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/FFM_DCall.txt', skiprows=1)
        stk0, dip0, rake0 = data[:,5], data[:,6],data[:,7]
        stk1, dip1, rake1 = data[:,8], data[:,9],data[:,10]

        strike, dip, rake = [],[],[]
        for i in range(len(stk0)):
            tmpstr, tmpdip, tmprake = utils.selectplane(model_para.strike[0], model_para.dip[0], stk0[i], dip0[i], rake0[i], stk1[i], dip1[i], rake1[i])
            strike.append(tmpstr)
            dip.append(tmpdip)
            rake.append(tmprake)

        data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/FFM_MT.txt')
        m1,m2,m3,m4,m5,m6 = data[:,4],data[:,5],data[:,6],data[:,7],data[:,8],data[:,9]


        for n in range(len(lon)):
            if lon[n] < 0:
                lon[n] = lon[n] + 360

        lonmin=np.min(lon)-lonmargin; lonmax=np.max(lon)+lonmargin; latmin=np.min(lat)-latmargin; latmax=np.max(lat)+latmargin

        m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
                  rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')

        x, y=m([lonmin, lonmax], [latmin, latmax])
        aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

        mapheight=axpwidth/aspect

        for axes in range(2):
            if axes == 0:
                axpxloc = axpxloc
            else:
                axp = ax.get_position()
                axpxloc = axp.x1+0.01

            ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
            axp=ax.get_position()
            #print(axp)
            m.fillcontinents(color='C7', zorder=0, alpha=0.1)

            cmap_slab2 = cm.batlow_r
            fh = Dataset(datarootdir+'work/ker_slab2_dep_02.24.18.grd', mode='r')
            lons = fh.variables['x'][:]; lats = fh.variables['y'][:]; tmax = fh.variables['z'][:] * -1
            fh.close()
            vmin, vmax = 0, 700
            x, y = m(lons, lats)
            CS = ax.contour(x, y, tmax, levels=np.arange(vmin, vmax+5, 5), colors='k', linewidths=0)
            tmpx = np.linspace(lonmin, lonmax, 11)
            tmpy = np.ones(len(tmpx))*-36.5
            manual_locations = [ (tmpx[i],tmpy[i]) for i in range(len(tmpx)) ]

            clabel = ax.clabel(CS, fontsize=0, inline=True, inline_spacing=0, manual=manual_locations, fmt='%d')
            CS = ax.contour(x, y, tmax, levels=np.arange(vmin, vmax+5, 5), colors='k', linewidths=0.75)

            for text in clabel[1:7]:
                #print(text)
                x, y = m(text.get_position()[0], text.get_position()[1])
                label = text.get_text()
                ax.text(x, y, label, ha='center', fontsize=8,
                        path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])


            slab2df = pd.read_csv(datarootdir+'work/ker_slab2_clp_02.24.18.csv', header=None, names=['lon', 'lat'])
            x, y = m(slab2df['lon'], slab2df['lat'])
            ax.plot(x, y, color='k', lw=0.75, linestyle='--')

            x, y = m(model_para.lon[0], model_para.lat[0])
            sc=ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=2000)
            sc.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])

            if axes == 0:
                dep0, dep1 = 0, 50
                axp = ax.get_position()
                if ptflag == 'P':
                    fig.text(axp.x1-0.01, axp.y0+0.01, 'P-axes: Shallow layers\n(0–50 km depth)', ha='right', va='bottom', fontsize=8,
                            path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])
                elif ptflag == 'T':
                    fig.text(axp.x1-0.01, axp.y0+0.01, 'T-axes: Shallow layers\n(0–50 km depth)', ha='right', va='bottom', fontsize=8,
                            path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])


            else:
                dep0, dep1 = 50, 1000
                axp = ax.get_position()
                if ptflag == 'P':
                    fig.text(axp.x1-0.01, axp.y0+0.01, 'P-axes: Deep layers\n(>50 km depth)', ha='right', va='bottom', fontsize=8,
                            path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])
                elif ptflag == 'T':
                    fig.text(axp.x1-0.01, axp.y0+0.01, 'T-axes: Deep layers\n(>50 km depth)', ha='right', va='bottom', fontsize=8,
                            path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])

            fig.text(axp.x1-0.01, axp.y1-0.01, 'Slab2 iso-depth (km)', ha='right', va='top', fontsize=8,
                    path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])
            for i in range(len(strike)):

                if depth[i] > dep0 and depth[i] <= dep1:
                    if lon[i] < 0:
                        lon[i] = lon[i] + 360
                    x, y = m(lon[i], lat[i])
                    focmecs=[strike[i], dip[i], rake[i]]
                    color=cmap(slip[i]/max(slip))

                    dx = model_para.xx[0]
                    focmecs=[m1[i],m2[i],m3[i],m4[i],m5[i],m6[i]]
                    focmecs = utils.convertUSEtoNED(focmecs)

                    mt = pmt.MomentTensor(
                        strike=strike[i],
                        dip=dip[i],
                        rake=rake[i])
                    pnt = utils.get_plunge_azimuth(mt.m6_up_south_east())
                    t_azi = pnt[5] # T-axis azimuth
                    p_azi = pnt[1] # p-axis azimuth
                    for pm in [-1, 1]:
                        if ptflag == 'P':
                            tmp = geod.Direct(lat[i], lon[i], p_azi, pm*20*1e3*slip[i]/max(slip))
                        elif ptflag == 'T':
                            tmp = geod.Direct(lat[i], lon[i], t_azi, pm*15*1e3*slip[i]/max(slip))
                        if tmp['lon1'] < 0: tmp['lon1'] += 360
                        if tmp['lon2'] < 0: tmp['lon2'] += 360
                        x0, y0 = m(tmp['lon1'], tmp['lat1'])
                        x1, y1 = m(tmp['lon2'], tmp['lat2'])
                        ax.plot([x0,x1], [y0,y1], color=cmap(slip[i]/max(slip)), lw=0.75,
                                zorder=int(slip[i]/max(slip)*100), alpha=slip[i]/max(slip))

            ax2 = utils.mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
            ax2.set_zorder(-1)
            if axes == 1:
                ax2.set_yticklabels([])

                axp = ax.get_position()
                cax=fig.add_axes([axp.x1+0.005, axp.y0, 0.01, axp.height/2])
                norm=mpl.colors.Normalize(vmin=0, vmax=max(slip))
                cb=mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, label='Slip (m)',
                                             ticks=np.linspace(0, max(slip), 5), format='%.2f', alpha=1)
            ax.set_facecolor('whitesmoke')

    plt.savefig(figname, bbox_inches="tight", pad_inches=0.1, dpi=300)
    plt.show()

    
    
def reloc(figname):
    # origin: OT 2021 03 04 13 27 35.705896 Lat -37.466117 Long 179.773495 Depth 71.972656
    elat, elon, edep = -37.466117, 179.773495, 71.972656
    xyloc = np.loadtxt(datarootdir+'hypocentreRelocation/vinti.20210304.132755.grid0.loc.scat.lonlat.XY')
    xzloc = np.loadtxt(datarootdir+'hypocentreRelocation/vinti.20210304.132755.grid0.loc.scat.lonlat.XZ')
    zyloc = np.loadtxt(datarootdir+'hypocentreRelocation/vinti.20210304.132755.grid0.loc.scat.lonlat.ZY')
    stacode = np.loadtxt(datarootdir+'hypocentreRelocation/vinti.sum.grid0.loc.stations', dtype=str, usecols=(0))

    data = np.loadtxt(datarootdir+'work/NZstation_info.txt', delimiter='|', usecols=(2,3), skiprows=1)
    lon, lat = data[:,1], data[:,0]
    code = np.loadtxt(datarootdir+'work/NZstation_info.txt', delimiter='|', usecols=(1), skiprows=1, dtype=str)
    df_station = pd.DataFrame(data = np.vstack([lon, lat, code]).T,
                       columns=['lon', 'lat', 'code'])

    lon, lat = xyloc[:,0], xyloc[:,1]
    for n in range(len(lon)):
        if lon[n] < 0: lon[n] += + 360
    lonmargin, latmargin = 1, 1
    lonmin=np.min(lon)-lonmargin; lonmax=np.max(lon)+lonmargin; latmin=np.min(lat)-latmargin; latmax=np.max(lat)+latmargin

    axpxloc, axpyloc = 0, 0
    tickintx,tickinty, tickformat = 1, 1, 0
    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
              rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')

    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

    axpwidth = 0.25
    mapheight=axpwidth/aspect

    fig=plt.figure(figsize=figsize)
    ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
    axpmap = ax.get_position()

    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.drawcoastlines(color='k', linewidth=0.5, zorder=10)

    x, y = m(lon, lat)
    ax.scatter(x, y, s=0.1, edgecolor='none', facecolor='k', alpha=1)

    x, y = m(elon, elat)
    sc=ax.scatter(x, y, s=10, marker='+', color='w', alpha=1, lw=0.5, zorder=100)

    ax.axhline(y, color='C7', lw=0.5, linestyle='--', zorder=0)
    ax.axvline(x, color='C7', lw=0.5, linestyle='--', zorder=0)

    src = shapefile.Reader(datarootdir+'work/tectonicplates/PB2002_boundaries.shp')
    for tmp in src.shapeRecords():
        x = [i[0] for i in tmp.shape.points[:]]
        y = [i[1] for i in tmp.shape.points[:]]
        for n in range(len(x)):
            if x[n] < 0:
                x[n] = x[n] + 360
        x, y = m(x, y)
        ax.plot(x, y, color='k', lw=0.3, linestyle='--', zorder=0)
    ax.plot([], [], color='k', lw=0.3, linestyle='--', zorder=0, label='Trench (Bird, 2003)')
    ax.legend(loc='upper left', fontsize=6)

    tmp = geod.Direct(-38, elon+0.05, 90, 50*1e3)
    if tmp['lon2'] < 0: tmp['lon2'] += 360
    x0, y0 = m(tmp['lon1'], tmp['lat1'])
    x1, y1 = m(tmp['lon2'], tmp['lat2'])
    ax.plot([x0, x1], [y0, y1], color='k', lw=0.75)
    x, y = m(tmp['lon1']+(tmp['lon2']-tmp['lon1'])/2, tmp['lat1']-0.05)
    ax.text(x, y, '50 km', ha='center', va='top', fontsize=8)

    ax2 = utils.mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
    #ax2.set_zorder(-1)
    ax2.xaxis.set_ticks_position('top')
    ax2.set_ylabel('Latitude')

    #print(geod.Inverse(latmin, lonmin, latmin, lonmax), axpmap)
    #print(geod.Inverse(latmin, elon, latmax, elon), axpmap)

    ###########################################################################################
    ###########################################################################################
    tmp = geod.Inverse(latmin, elon, latmax, elon)
    xzheight = (150/(tmp['s12']*1e-3)) * axpmap.height

    tmp = geod.Inverse(latmin, lonmin, latmin, lonmax)
    zywidth = (150 * axpmap.width) / (tmp['s12']*1e-3)
    #print(zywidth)


    ax=fig.add_axes([axpmap.x0, axpmap.y0-xzheight-0.01, axpmap.width, xzheight])
    lon, z = xzloc[:,0], xzloc[:,1]
    for n in range(len(lon)):
        if lon[n] < 0: lon[n] += + 360
    ax.scatter(lon, z, s=0.1, edgecolor='none', facecolor='k', alpha=1)
    sc=ax.scatter(elon, edep, s=10, marker='+', color='w', alpha=1, lw=0.5, zorder=100)
    ax.axhline(edep, color='C7', lw=0.5, linestyle='--', zorder=0)
    ax.axvline(elon, color='C7', lw=0.5, linestyle='--', zorder=0)

    ax.set_xlim(lonmin, lonmax)
    tmp = ax.get_xticks()
    label = [ str(int(tmp[i]))+'$\degree$E' for i in range(len(tmp)) ]
    ax.set_xticklabels(label)
    ax.set_yticks(np.arange(0, 200, 50))
    ax.set_ylim(0, 150)
    ax.invert_yaxis()
    ax.set_ylabel('Depth (km)')
    ax.set_xlabel('Longitude')
    ###########################################################################################
    ###########################################################################################
    ax=fig.add_axes([axpmap.x1+0.01, axpmap.y0, zywidth, axpmap.height])
    z, lat = zyloc[:,0], zyloc[:,1]
    ax.scatter(z, lat, s=0.1, edgecolor='none', facecolor='k', alpha=1)
    sc=ax.scatter(edep, elat, s=10, marker='+', color='w', alpha=1, lw=0.5, zorder=100)
    ax.axhline(elat, color='C7', lw=0.5, linestyle='--', zorder=0)
    ax.axvline(edep, color='C7', lw=0.5, linestyle='--', zorder=0)

    ax.yaxis.set_ticks_position('right')
    ax.set_xticks(np.arange(0, 200, 50))
    ax.set_xlim(0, 150)
    ax.set_ylim(latmin, latmax)
    tmp = ax.get_yticks()
    label = [ str(int(tmp[i]))+'$\degree$S' for i in range(len(tmp)) ]
    ax.set_yticklabels(label)
    ax.set_xlabel('Depth (km)')
    ######################################################################################
    ######################################################################################
    ######################################################################################
    lon, lat = xyloc[:,0], xyloc[:,1]
    for n in range(len(lon)):
        if lon[n] < 0: lon[n] += + 360

    lonmargin, latmargin = 10, 10

    axp = ax.get_position()
    axpxloc, axpyloc = axp.x1+0.1, axpmap.y0-xzheight-0.01
    tickintx,tickinty, tickformat = 5, 5, 0

    lonmin=np.min(lon)-lonmargin; lonmax=np.max(lon)+2; latmin=np.min(lat)-latmargin; latmax=np.max(lat)+latmargin

    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
              rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')

    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

    mapheight = axpmap.height+xzheight+0.01
    mapwidth = mapheight * aspect
    ax=fig.add_axes([axpxloc, axpyloc, mapwidth, mapheight])

    m.fillcontinents(color='C7', zorder=0, alpha=0.2)
    m.drawcoastlines(color='k', linewidth=0.75, zorder=10)

    #x, y = m(lon, lat)
    #ax.scatter(x, y, s=2, edgecolor='none', facecolor='C7', alpha=0.1)

    x, y = m(elon, elat)
    sc=ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=100, label='Relocated epicentre')
    sc.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])

    for j in range(len(stacode)):
        if stacode[j] == 'RAO':
            slon, slat = 177.929000, -29.245000
        elif stacode[j] == 'SNZO':
            slon, slat = 174.7046, -41.3101

        else:
            slat = float(df_station.loc[(df_station['code'] == stacode[j])]['lat'].values[0])
            slon = float(df_station.loc[(df_station['code'] == stacode[j])]['lon'].values[0])


        if slon < 0: slon += 360
        x, y = m(slon, slat)
        ax.scatter(x, y, marker='^', facecolor='ivory', edgecolor='k', lw=0.75, s=20, zorder=10)


    src = shapefile.Reader(datarootdir+'work/tectonicplates/PB2002_boundaries.shp')
    for tmp in src.shapeRecords():
        x = [i[0] for i in tmp.shape.points[:]]
        y = [i[1] for i in tmp.shape.points[:]]
        for n in range(len(x)):
            if x[n] < 0:
                x[n] = x[n] + 360
        x, y = m(x, y)
        ax.plot(x, y, color='k', lw=0.3, linestyle='--', zorder=0)

    ax2 = utils.mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
    #ax2.set_zorder(-1)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')

    plt.savefig(figname, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.show()

    
def relocAftershock(figname, original, relocated, starttime, endtime):
    fig = plt.figure(figsize=figsize)
    for num,inputcatalog in enumerate([original, relocated]):
        model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
        j = 13 # optimum model
        model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
        elon, elat, edep = model_para.lon[0], model_para.lat[0], model_para.depth[0]

        data = np.loadtxt(inputcatalog, usecols=(8,7,9))
        np.savetxt(datarootdir+'work/tmp.txt', data)
        print(len(data))

        command = ['gmt','project',datarootdir+'work/tmp.txt','-C'+str(elon)+'/'+str(elat),'-A20','-Q','>',
                   datarootdir+'work/tmp_proj_20.txt']
        process = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        data = np.loadtxt(inputcatalog, usecols=(1,2,3,4,5,6,10))
        year,month,day,hour,minute,sec,mag = data[:,0],data[:,1],data[:,2],data[:,3],data[:,4],data[:,5],data[:,6]

        for n in range(len(sec)):
            if sec[n] == 60.0: sec[n] = 59.99999

        tmp_time = np.array([ obspy.UTCDateTime(str(int(year[i]))+'-'+str(int(month[i]))+'-'+str(int(day[i]))+'T'+str(int(hour[i]))+':'+str(int(minute[i]))+':'+str(sec[i])) for i in range(len(year)) ])

        data = np.loadtxt(datarootdir+'work/tmp_proj_20.txt')
        lon, lat, dep, x, y = data[:,0],data[:,1],data[:,2],data[:,3],data[:,4]*-1
        for n in range(len(lon)):
            if lon[n] < 0: lon[n] += 360

        df = pd.DataFrame(data = np.vstack([tmp_time, lat, lon, dep, mag, x, y]).T,
                           columns=['origintime', 'latitude', 'longitude', 'depth', 'magnitude', 'projx', 'projy'])

        lonmin, lonmax, latmin, latmax=178.1, 181.3, -39.5, -36
        m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
                  rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')
        x, y=m([lonmin, lonmax], [latmin, latmax])
        aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

        axpxloc, axpyloc, axpwidth = 0.1, 0.1,0.5
        if num == 1:
            axpyloc = axp0.y0-axp.height-0.1
        mapheight=axpwidth/aspect

        ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
        axp=ax.get_position()
        if num == 0:
            axp0 = axp
            label = 'Original location (GeoNet)'
        elif num == 1:
            label = 'Relocation'
        fig.text(axp.x0+0.01, axp.y1-0.01, label, ha='left', va='top', size=8,
               path_effects=[path_effects.Stroke(linewidth=3, foreground='w', alpha=1), path_effects.Normal()])
        fig.text(axp.x0+0.01, axp.y0+0.01, 'From '+starttime[0:10]+' to '+endtime[0:10], ha='left', va='bottom', size=8,
               path_effects=[path_effects.Stroke(linewidth=3, foreground='w', alpha=1), path_effects.Normal()])

        x, y = m(model_para.lon[0], model_para.lat[0])
        ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=10, label='Relocated epicentre (This study)',
                  path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])

        tmpdf = df[(df['longitude'] >= lonmin) & (df['longitude'] <= lonmax) &
                   (df['origintime'] >= starttime) & (df['origintime'] <= endtime) &
                   (df['magnitude'] >= 0)]
        
        print(len(tmpdf))

        x, y = m(tmpdf['longitude'], tmpdf['latitude'])
        sc = ax.scatter(x, y, s=1, marker='o', facecolor='k', edgecolor='none', alpha=1, zorder=1)

        src = shapefile.Reader(datarootdir+'work/tectonicplates/PB2002_boundaries.shp')
        for tmp in src.shapeRecords():
            x = [i[0] for i in tmp.shape.points[:]]
            y = [i[1] for i in tmp.shape.points[:]]
            for n in range(len(x)):
                if x[n] < 0: x[n] = x[n] + 360
            x, y = m(x, y)
            ax.plot(x, y, color='C7', lw=0.5, linestyle='--', zorder=0)
        ax.plot([], [], color='C7', lw=0.5, linestyle='--', label='Trench (Bird, 2003)')
        ax2 = utils.mapTicksBasemap(fig,m,ax,1,1,lonmin,lonmax,latmin,latmax,0)

        label = ['S', 'N', 'W', 'E']
        num = 0
        for az in [20, 110]:
            for pm in [-1, 1]:
                tmp = geod.Direct(model_para.lat[0], model_para.lon[0], az, pm*100*1e3)
                if tmp['lon2'] < 0: tmp['lon2'] += 360
                x0, y0 = m(tmp['lon1'], tmp['lat1'])
                x1, y1 = m(tmp['lon2'], tmp['lat2'])
                ax.plot([x0, x1], [y0, y1], lw=0.75, color='k', linestyle='--')
                if az == 20 and pm == -1:
                    ha, va = 'right', 'top'
                    color='k'
                elif az == 20 and pm == 1:
                    ha, va = 'left', 'bottom'
                    color='k'
                elif az == 110 and pm == -1:
                    ha, va = 'right', 'center'
                    color='k'
                elif az == 110 and pm == 1:
                    ha, va = 'left', 'center'
                    color='k'
                ax.text(x1, y1, label[num], ha=ha, va=va, color=color,
                       path_effects=[path_effects.Stroke(linewidth=3, foreground='w', alpha=0.5), path_effects.Normal()])

                num += 1

        for xpanel in range(2):
            if xpanel == 0:
                axp = ax.get_position()
                ax = fig.add_axes([axp.x1+0.01, axp.y1-(axp.height/2.05), 0.3, axp.height/2.05])
                x = tmpdf['projx']
                rlabel = 'N'
                llabel = 'S'
            else:
                ax = fig.add_axes([axp.x1+0.01, axp.y0, 0.3, axp.height/2.05])
                x = tmpdf['projy']
                rlabel = 'E'
                llabel = 'W'

            ax.scatter(x, tmpdf['depth'], s=1, marker='o', facecolor='k', edgecolor='none', alpha=1, zorder=1)

            ax.scatter(0, edep, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=0, label='Relocated epicentre (This study)',
                          path_effects=[path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])

            ax.text(-100, 153, llabel, ha='center', va='bottom', color='k',
                   path_effects=[path_effects.Stroke(linewidth=3, foreground='w', alpha=0.5), path_effects.Normal()])
            ax.text(100, 153, rlabel, ha='center', va='bottom', color='k',
                   path_effects=[path_effects.Stroke(linewidth=3, foreground='w', alpha=0.5), path_effects.Normal()])

            ax.set_xlim(-120, 120)
            ax.set_ylim(0, 155)
            ax.invert_yaxis()
            ax.set_ylabel('Depth (km)')
            ax.yaxis.set_ticks_position('right')
            ax.yaxis.set_label_position('right')
            if xpanel == 1:
                ax.set_xlabel('Distance (km)')
            else:
                ax.set_xticklabels([])
                
                
            histflag = 0
            if histflag == 1:
                _tmpdf = df[(df['longitude'] >= lonmin) & (df['longitude'] <= lonmax) &
                           (df['origintime'] >= starttime) & (df['origintime'] <= endtime) &
                           (df['depth'] != 12.0) & (df['depth'] != 33.0 ) ]

                _axp = ax.get_position()
                axhist = fig.add_axes([_axp.x1+0.1, _axp.y0, 0.15, _axp.height])
                bins = np.arange(0, 1000, 5)
                axhist.hist(_tmpdf['depth'], bins=bins, orientation='horizontal', color='C7')
                axhist.set_ylim(0, 155)
                #axhist.set_xlim(0, 200)
                axhist.invert_yaxis()
                axhist.set_ylabel('Depth (km)')
                axhist.yaxis.set_ticks_position('right')
                axhist.yaxis.set_label_position('right')
                axhist.set_xlabel('Count')
                axhist.set_xscale('log')
            


    plt.savefig(figname, bbox_inches="tight", pad_inches=0.1, dpi=300)
    plt.show()

    
    
def seismicityFocalMech(figname):
    url = 'https://raw.githubusercontent.com/GeoNet/data/main/moment-tensor/GeoNet_CMT_solutions.csv'
    dfGeoNetCMTData = pd.read_csv(url)
    dfGeoNetCMTData.Date = pd.to_datetime(dfGeoNetCMTData.Date, format='%Y%m%d%H%M%S')
    dfGeoNetCMTData.Longitude[dfGeoNetCMTData.Longitude < 0] += 360

    fig = plt.figure(figsize=figsize)

    lonmin, lonmax, latmin, latmax=178.1, 181, -39.5, -36
    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
              rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')
    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

    axpxloc, axpyloc, axpwidth = 0.1, 0.1,0.4
    mapheight=axpwidth/aspect
    cmap = cm.batlow

    for j, catalogflag in enumerate(['geonet', 'gcmt']):
        if catalogflag == 'gcmt':
            ax=fig.add_axes([axp.x1+0.01, axp.y0, axpwidth, mapheight])
            axp=ax.get_position()
            fig.text(axp.x0+0.01, axp.y1-0.01, 'GCMT (Mw≥5, ≤2021-03-11)', ha='left', va='top', fontsize=8,
                    path_effects=[path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])
        else:
            ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
            axp=ax.get_position()
            fig.text(axp.x0+0.01, axp.y1-0.01, 'GeoNet (Mw≥5, ≤2021-03-11)', ha='left', va='top', fontsize=8,
                    path_effects=[path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])

        m.fillcontinents(color='C7', zorder=0, alpha=0.1)
        m.drawcoastlines(color='k', linewidth=1, zorder=10)

        model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
        j = 13
        model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')

        x, y = m(model_para.lon[0], model_para.lat[0])
        ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=10, label='Relocated epicentre (This study)',
                  path_effects=[path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])

        if catalogflag == 'geonet':

            for n in range(2):
                if n == 0:
                    df = dfGeoNetCMTData[(dfGeoNetCMTData['Date'] >= '1900-03-04') & (dfGeoNetCMTData['Date'] < '2021-03-04') &
                                         (dfGeoNetCMTData['Mw'] >= 5) & (dfGeoNetCMTData['CD'] <= 10000)]
                    color_t = 'C7'
                    zorder_base = 0
                else:
                    df = dfGeoNetCMTData[(dfGeoNetCMTData['Date'] >= '2021-03-04') & (dfGeoNetCMTData['Date'] <= '2021-03-11') &
                                         (dfGeoNetCMTData['Mw'] >= 5) & (dfGeoNetCMTData['CD'] <= 10000)]
                    color_t = 'C3'
                    zorder_base = 1000
                cmap=cm.batlow_r
                for i in range(df.shape[0]):
                    focmecs=[df.strike1.values[i],df.dip1.values[i],df.rake1.values[i]]
                    x, y = df.Longitude.values[i], df.Latitude.values[i]
                    focmecs = [df.Mxx.values[i],
                               df.Myy.values[i],
                               df.Mzz.values[i],
                               df.Mxy.values[i],
                               df.Mxz.values[i],
                               df.Myz.values[i]]
                    tmp = beachball.plot_beachball_mpl(focmecs, ax, size=6, position=(x, y),
                                             beachball_type='deviatoric', edgecolor='none', color_t=color_t,
                                             color_p='w', linewidth=0.5, alpha=1, zorder=int(100-df.Mw.values[i]*10)+zorder_base)
                    tmp = beachball.plot_beachball_mpl(focmecs, ax, size=6, position=(x, y),
                                             beachball_type='dc', edgecolor='k', color_t='none',
                                             color_p='none', linewidth=0.5, alpha=1, zorder=int(100-df.Mw.values[i]*10)+zorder_base)
                    if df.Mw.values[i] == 7.3:
                        tmp = beachball.plot_beachball_mpl(focmecs, ax, size=6, position=(x, y),
                                                 beachball_type='deviatoric', edgecolor='none', color_t='C5',
                                                 color_p='w', linewidth=0.5, alpha=1, zorder=int(100-df.Mw.values[i]*10)+zorder_base)
                        tmp = beachball.plot_beachball_mpl(focmecs, ax, size=6, position=(x, y),
                                                 beachball_type='dc', edgecolor='k', color_t='none',
                                                 color_p='none', linewidth=0.5, alpha=1, zorder=int(100-df.Mw.values[i]*10)+zorder_base)

        elif catalogflag == 'gcmt':
            # GCMT
            #'''
            cat = obspy.read_events(datarootdir+'work/SPUD_QUAKEML_bundle.xml')
            for event in cat.events:
                mag = event.magnitudes[0].mag
                lon, lat, dep = event.origins[1].longitude, event.origins[1].latitude, event.origins[1].depth
                if lon < 0: lon += 360
                if (event.origins[1].time >= obspy.UTCDateTime('1900-03-04')) and (event.origins[1].time <= obspy.UTCDateTime('2021-03-11')) and \
                mag >= 5:

                    if event.origins[1].time < obspy.UTCDateTime('2021-03-04'):
                        color_t = 'C7'
                        zorder_base = 0
                    else:
                        color_t = 'C3'
                        zorder_base = 1000

                    focmecs=[event.focal_mechanisms[0].moment_tensor.tensor.m_rr,
                             event.focal_mechanisms[0].moment_tensor.tensor.m_tt,
                             event.focal_mechanisms[0].moment_tensor.tensor.m_pp,
                             event.focal_mechanisms[0].moment_tensor.tensor.m_rt,
                             event.focal_mechanisms[0].moment_tensor.tensor.m_rp,
                             event.focal_mechanisms[0].moment_tensor.tensor.m_tp]
                    focmecs = utils.convertUSEtoNED(focmecs)
                    x, y = m(lon, lat)

                    size = 6

                    tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                             beachball_type='deviatoric', edgecolor='none', color_t=color_t,
                                             color_p='w', linewidth=0.5, alpha=1, zorder=int(100-mag*10)+zorder_base)
                    tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                             beachball_type='dc', edgecolor='k', color_t='none',
                                             color_p='none', linewidth=0.5, alpha=1, zorder=int(100-mag*10)+zorder_base)
                    if mag == 7.2:
                        #print(lon, lat)
                        #print(event)
                        #print(focmecs)
                        tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                                 beachball_type='deviatoric', edgecolor='none', color_t='C5',
                                                 color_p='w', linewidth=0.5, alpha=1, zorder=int(100)+zorder_base+1000)
                        tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                                 beachball_type='dc', edgecolor='k', color_t='none',
                                                 color_p='none', linewidth=0.5, alpha=1, zorder=int(100)+zorder_base+1000)
            #'''


        src = shapefile.Reader(datarootdir+'work/tectonicplates/PB2002_boundaries.shp')
        for tmp in src.shapeRecords():
            x = [i[0] for i in tmp.shape.points[:]]
            y = [i[1] for i in tmp.shape.points[:]]
            for n in range(len(x)):
                if x[n] < 0:
                    x[n] = x[n] + 360
            x, y = m(x, y)
            ax.plot(x, y, color='C7', lw=0.5, linestyle='--', zorder=0)
        ax.plot([], [], color='C7', lw=0.5, linestyle='--', label='Trench (Bird, 2003)')
        ax2 = utils.mapTicksBasemap(fig,m,ax,1,1,lonmin,lonmax,latmin,latmax,0)

        ax.legend(loc='lower right', fontsize=6).set_zorder(1000)
        if catalogflag == 'gcmt':
            ax2.set_yticklabels([])

    plt.savefig(figname, bbox_inches="tight", pad_inches=0.1, dpi=300)
    plt.show()
    
    
def wphase(figname):
    # /Users/ryo/GoogleDrive/Work/2021NewZealand/2021-03-04 New Zealand event shared folder/WphaseCMT/report/gfz2021ekhv/wphase_cmt_gfz2021ekhv/plots/fits_waveform/default
    # grep td_tele.wphase  fits_waveform.default.plot_group.yaml | grep - | sed -e "s/- td_tele.wphase.//g" > station_used.txt
    datadir=datarootdir+'WphaseCMT/report/gfz2021ekhv/wphase_cmt_gfz2021ekhv/plots/fits_waveform/default'
    invdir = datarootdir+'WphaseCMT/stations'
    data = np.loadtxt(datadir+'/station_used.txt', dtype=str)


    stationlist = []
    for j in range(len(data)):
        network = data[j].split('.')[0]
        stationcode = data[j].split('.')[1]
        inv = read_inventory(invdir+'/'+network+'.'+stationcode+'.xml')
        lon, lat = inv[0][0].longitude, inv[0][0].latitude
        stationlist.append([stationcode, lat, lon])
    #print(len(stationlist))


    stationlist_unique = utils.get_unique_list(stationlist)
    stationlist_unique

    from pyrocko.model import load_events, load_stations, Channel
    datadir=datarootdir+'WphaseCMT'
    res = load_events(datadir+"/ensemble.yaml")  # Read in pyrocko event file containing ensemble events
    depths = [e.depth/1000 for e in res]  # Make list of depths, convert to m->km

    mean = load_events(datadir+"/mean_wphase.yaml")[0]  # Read in event file containing best solution
    mean_depth = mean.depth/1000

    fig=plt.figure(figsize=figsize)
    ax=fig.add_axes([0.1, 0.1, 0.3, 0.5])
    bins = np.arange(0, 102, 2)
    ax.hist(depths, orientation='horizontal', bins=bins, density=True,
            color='C7', edgecolor='w', linewidth=0.5, label='Sample (2-km bin)')

    ax.axhline(mean_depth, color='k', zorder=1, linestyle='--', lw=1, label='Mean')
    ax.set_ylim(0, 100)
    ax.invert_yaxis()
    ax.set_xlabel('Normed counts')
    ax.set_ylabel('Depth (km)')
    ax.legend(fontsize=8, loc='upper left')


    meclocs, mos, focmecs, labels = utils.loadmecs(datarootdir)

    lonmargin, latmargin = 0.2, 0.2
    tickintx,tickinty, tickformat = 0.5, 0.5, 1
    lonmin=179.5; lonmax=180.2; latmin=-37.75; latmax=-37
    #print(lonmin, lonmax, latmin, latmax)

    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
              rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')

    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

    axp = ax.get_position()
    axpwidth = axp.width*0.75
    mapheight=axpwidth/aspect
    #mapheight = axp.height
    #axpwidth = mapheight * aspect

    axpxloc, axpyloc = axp.x1+0.12, axp.y1-mapheight

    ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.drawcoastlines(color='k', linewidth=0.5, zorder=10)

    src = shapefile.Reader(datarootdir+'work/tectonicplates/PB2002_boundaries.shp')
    for tmp in src.shapeRecords():
        x = [i[0] for i in tmp.shape.points[:]]
        y = [i[1] for i in tmp.shape.points[:]]
        for n in range(len(x)):
            if x[n] < 0:
                x[n] = x[n] + 360
        x, y = m(x, y)
        ax.plot(x, y, color='C7', lw=0.5, linestyle='--')

    x, y = m(179.875, -37.4)
    ax.text(x, y, 'Trench', ha='left', va='center', fontsize=8, color='gray', rotation=65,
            path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])

    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    model_para = utils.load_fort40(datarootdir+'model_'+str(model[13])+'/fort.40')

    x, y = m(model_para.lon[0], model_para.lat[0])
    sc=ax.scatter(x, y, s=200, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=100)
    sc.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])


    focmec = focmecs[2]
    x, y = m(meclocs[2][0], meclocs[2][1])
    tmp = beachball.plot_beachball_mpl(focmec, ax, size=20, position=(x, y),
                             beachball_type='deviatoric', edgecolor='none', color_t='C5',
                             color_p='w', linewidth=0, alpha=1, zorder=2)
    tmp.set_clip_on(False)
    tmp = beachball.plot_beachball_mpl(focmec, ax, size=20, position=(x, y),
                             beachball_type='dc', edgecolor='k', color_t='none',
                             color_p='none', linewidth=0.75, alpha=1, zorder=2)

    ax.set_zorder(11)
    ax2 = utils.mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
    ax2.set_zorder(10)

    ax=fig.add_axes([axp.x1+0.1, axp.y0-0.06, axp.width*0.85, axp.width*0.85])
    axp=ax.get_position()
    sc = utils.aziequi_wphase(ax, stationlist_unique, meclocs[2], 10)

    plt.savefig(figname, bbox_inches="tight", pad_inches=0.1, dpi=300)
    plt.show()

    
    
def rcmtLF(figname):
    # awk 'NR >= 4 && NR <= 38 {print $0}' < inv1.dat > inv1_extracted.txt
    data = np.loadtxt(datarootdir+'R-CMT/LowFreq_SingleSubEvent_LineDepths/sources.gmt', skiprows=1)
    dep = data[:,2]

    data = np.loadtxt(datarootdir+'R-CMT/LowFreq_SingleSubEvent_LineDepths/inv1_extracted.txt')
    corr, dc, str1, dip1, rake1, ishift = data[:,2], data[:,4], data[:,5], data[:,6], data[:,7], data[:,1]

    sampling_rate = 0.04 # Hz

    df = pd.DataFrame(data = np.vstack([corr, dep, dc, str1, dip1, rake1, ishift]).T,
                       columns=['corr', 'dep', 'dc', 'str1', 'dip1', 'rake1', 'ishift'])
    df = df.sort_values(by='corr', ascending=True).reset_index(drop=True)

    fig=plt.figure(figsize=figsize)
    ax=fig.add_axes([0.1, 0.1, 0.3, 0.69])

    cmap = cm.roma_r
    for j in range(len(dep)):
        focmec = [df['str1'][j],df['dip1'][j],df['rake1'][j]]
        x, y = (df['corr'][j], df['dep'][j])
        ax.scatter(x, y, zorder=0, s=0)
        tmp = beachball.plot_beachball_mpl(focmec, ax, size=10, position=(df['corr'][j], df['dep'][j]),
                                 beachball_type='dc', edgecolor='k', color_t=cmap(df['dc'][j]/100),
                                 color_p='w', linewidth=0.75)

    ax.axhline(df.loc[df['corr'].idxmax(), 'dep'], color='C7', linestyle='--', zorder=0, lw=0.75)

    timeshift = df.loc[df['corr'].idxmax(), 'ishift'] * sampling_rate

    with open(datarootdir+'R-CMT/LowFreq_SingleSubEvent_LineDepths/inv1.dat') as f:
        data = f.readlines()[51]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]

    mw = 2.0 / 3.0 * (np.log10(moment) - 9.1)
    axp = ax.get_position()
    fig.text(axp.x0+0.01, axp.y0+0.01, 'Mw '+str('{:.1f}'.format(mw))+' | OT+'+str('{:.1f}'.format(timeshift))+' s',
             va='bottom', fontsize=8)

    ax.set_xlim(df.loc[df['corr'].idxmin(), 'corr']-0.01, df.loc[df['corr'].idxmax(), 'corr']+0.01)
    ax.set_ylim(-5, 150)
    ax.invert_yaxis()
    ax.set_xlabel('Correlation')
    ax.set_ylabel('Depth (km)')

    axp = ax.get_position()
    cax=fig.add_axes([axp.x1+0.005, axp.y0, 0.01, axp.height*0.5])
    norm=mpl.colors.Normalize(vmin=0, vmax=100)
    cb=mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, label='DC (%)')


    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################

    axpbase = ax.get_position()

    datadir=datarootdir+'R-CMT/LowFreq_SingleSubEvent_LineDepths'
    df = pd.read_table(datadir+'/allstat.dat',
                     header=None, names=['code', 'used', 'NS', 'EW', 'Z', 'f1', 'f2', 'f3', 'f4'], delim_whitespace=True)

    numsta = 0
    axpwidth = 0.3
    axpheight= 0.06
    axpy0 = axpbase.y1-0.06
    axpx0 = axpbase.x1+0.13
    normamp = 1e-2*0.3

    obsnormamp = []
    synnormamp = []

    usedstationindex = []
    for j in range(len(df)):
        if df['used'][j] == 1: usedstationindex.append(j)


    for j in usedstationindex:

        obsdata = np.loadtxt(datadir+'/waveform_fits/'+df['code'][j]+'fil.dat')
        syndata = np.loadtxt(datadir+'/waveform_fits/'+df['code'][j]+'syn.dat')

        ax1 = fig.add_axes([axpx0, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        ax2 = fig.add_axes([axpx0+axpwidth+0.01, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        ax3 = fig.add_axes([axpx0+(axpwidth+0.01)*2, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        if numsta == 0:
            ax1.text(175, 1.6, 'NS', ha='center', va='bottom')
            ax2.text(175, 1.6, 'EW', ha='center', va='bottom')
            ax3.text(175, 1.6, 'Z', ha='center', va='bottom')

        ax1.set_yticklabels([])
        ax2.set_yticklabels([])
        ax3.set_yticklabels([])

        if j < len(usedstationindex)-1:
            ax1.set_xticklabels([])
            ax2.set_xticklabels([])
            ax3.set_xticklabels([])

        ax1.set_xlim(-15, 350)
        ax2.set_xlim(-15, 350)
        ax3.set_xlim(-15, 350)

        ax1.set_ylim(-1.5, 1.5)
        ax2.set_ylim(-1.5, 1.5)
        ax3.set_ylim(-1.5, 1.5)

        #normamp = max(obsdata[:,1]) * 3
        obsnormamp = np.max([np.max(obsdata[:,1]), np.max(obsdata[:,2]), np.max(obsdata[:,3])])
        synnormamp = np.max([np.max(syndata[:,1]), np.max(syndata[:,2]), np.max(syndata[:,3])])

        if df['NS'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        ax1.plot(obsdata[:,0], obsdata[:,1]/obsnormamp, color=color, lw=1)
        ax1.plot(syndata[:,0], syndata[:,1]/synnormamp, color='r', lw=1, alpha=0.85)

        if df['EW'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        #normamp = max(obsdata[:,2]) * 3
        ax2.plot(obsdata[:,0], obsdata[:,2]/obsnormamp, color=color, lw=1)
        ax2.plot(syndata[:,0], syndata[:,2]/synnormamp, color='r', lw=1, alpha=0.85)

        if df['Z'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        #normamp = max(obsdata[:,3]) * 3
        ax3.plot(obsdata[:,0], obsdata[:,3]/obsnormamp, color=color, lw=1)
        ax3.plot(syndata[:,0], syndata[:,3]/synnormamp, color='r', lw=1, alpha=0.85)

        axp = ax1.get_position()
        fig.text(axp.x0+0.005, axp.y1-0.005, df['code'][j], va='top', ha='left', fontsize=8,
                path_effects=[path_effects.Stroke(linewidth=1 , foreground='w', alpha=1), path_effects.Normal()])

        #print(df['code'][j])

        numsta += 1

    ax1.set_xlabel('Time (s)')

    ax3.plot([], [], color='k', lw=1, label='Obs.')
    ax3.plot([], [], color='C7', lw=1, label='Obs. (unused)')
    ax3.plot([], [], color='r', lw=1, label='Syn.', alpha=0.85)
    ax3.legend(fontsize=6, loc=(0, -1.2), ncol=3, borderaxespad=0., handlelength=1)
    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################

    data = np.loadtxt(datarootdir+'work/NZstation_info.txt', delimiter='|', usecols=(2,3), skiprows=1)
    lon, lat = data[:,1], data[:,0]
    code = np.loadtxt(datarootdir+'work/NZstation_info.txt', delimiter='|', usecols=(1), skiprows=1, dtype=str)

    df_station = pd.DataFrame(data = np.vstack([lon, lat, code]).T,
                       columns=['lon', 'lat', 'code'])

    datadir=datarootdir+'R-CMT/LowFreq_SingleSubEvent_LineDepths'
    df = pd.read_table(datadir+'/allstat.dat',
                     header=None, names=['code', 'used', 'NS', 'EW', 'Z', 'f1', 'f2', 'f3', 'f4'], delim_whitespace=True)

    for n in range(len(lon)):
        if lon[n] < 0:
            lon[n] = lon[n] + 360

    lonmargin, latmargin = 0.1, 0.1
    tickintx,tickinty, tickformat = 5, 3, 0

    lonmin=np.min(lon)-lonmargin; lonmax=np.max(lon)+lonmargin; latmin=np.min(lat)-latmargin; latmax=np.max(lat)+latmargin
    elon, elat = 179.8000000000, -37.3700000000
    lonmin, lonmax, latmin, latmax = elon-5, elon+1, elat-5, elat+2
    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
              rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')

    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

    axpwidth = 0.15
    mapheight=axpwidth/aspect

    axpxloc, axpyloc = axpbase.x1-axpwidth+0.085, axpbase.y1-mapheight

    ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.drawcoastlines(color='k', linewidth=0.5, zorder=10)

    for j in range(len(df)):
        if df['used'][j] == 1:
            if df['code'][j] == 'RAO':
                lon, lat = 177.929000, -29.245000
            else:
                lat = float(df_station.loc[(df_station['code'] == df['code'][j])]['lat'].values[0])
                lon = float(df_station.loc[(df_station['code'] == df['code'][j])]['lon'].values[0])
            x, y = m(lon, lat)
            ax.scatter(x, y, marker='^', facecolor='ivory', edgecolor='k', lw=0.75, s=50, zorder=10)
            #ax.text(x, y, df['code'][j])

    x, y = m(elon, elat)
    ax.scatter(x, y, marker='x', lw=2, s=70, color='k')

    ax.set_zorder(11)
    ax2 = utils.mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
    ax2.set_zorder(10)

    plt.savefig(figname, bbox_inches="tight", pad_inches=0.1, dpi=150)
    plt.show()

    
    
def rcmtHF(figname):
    # awk 'NR >= 4 && NR <= 28 {print $0}' < inv1.dat > inv1_extracted_sub1.txt
    # awk 'NR >= 61 && NR <= 85 {print $0}' < inv1.dat > inv1_extracted_sub2.txt
    # awk 'NR >= 118 && NR <= 142 {print $0}' < inv1.dat > inv1_extracted_sub3.txt
    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths/inv1.dat') as f:
        data = f.readlines()[41]
    moment0 = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]

    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths/inv1.dat') as f:
        data = f.readlines()[98]
    moment1 = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]

    moments = [moment0, moment1]

    data = np.loadtxt(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths/sources.gmt', skiprows=1)
    dep = data[:,2]

    fig=plt.figure(figsize=figsize)

    sampling_rate = 0.04 # Hz

    for i in range(2):
        data = np.loadtxt(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths/inv1_extracted_sub'+str(i+1)+'.txt')
        corr, dc, str1, dip1, rake1, ishift = data[:,2], data[:,4], data[:,5], data[:,6], data[:,7],data[:,1]

        df = pd.DataFrame(data = np.vstack([corr, dep, dc, str1, dip1, rake1, ishift]).T,
                           columns=['corr', 'dep', 'dc', 'str1', 'dip1', 'rake1', 'ishift'])
        df = df.sort_values(by='corr', ascending=True).reset_index(drop=True)

        ax=fig.add_axes([0.1+i*(0.2+0.01), 0.1, 0.2, 0.83])

        cmap = cm.roma_r
        for j in range(len(dep)):
            focmec = [df['str1'][j],df['dip1'][j],df['rake1'][j]]
            x, y = (df['corr'][j], df['dep'][j])
            ax.scatter(x, y, zorder=0, s=0)
            tmp = beachball.plot_beachball_mpl(focmec, ax, size=10, position=(df['corr'][j], df['dep'][j]),
                                     beachball_type='dc', edgecolor='k', color_t=cmap(df['dc'][j]/100),
                                     color_p='w', linewidth=0.75)

        ax.axhline(df.loc[df['corr'].idxmax(), 'dep'], color='C7', linestyle='--', zorder=0, lw=0.75)
        timeshift = df.loc[df['corr'].idxmax(), 'ishift'] * sampling_rate
        mw = 2.0 / 3.0 * (np.log10(moments[i]) - 9.1)
        axp = ax.get_position()
        fig.text(axp.x0+0.01, axp.y0+0.01, 'Sub-event '+str(i+1)+'\nMw '+str('{:.1f}'.format(mw))+' | OT+'+str('{:.1f}'.format(timeshift))+' s',
                 va='bottom', fontsize=8)

        ax.set_xlim(np.min(corr)-0.05, np.max(corr)+0.05)
        ax.set_ylim(0, 113)
        ax.invert_yaxis()
        if i >= 1:
            ax.set_yticklabels([])
        if i == 0:
            ax.set_ylabel('Depth (km)')
            ax.set_xlabel('Correlation')

    axp = ax.get_position()
    cax=fig.add_axes([axp.x1+0.005, axp.y0, 0.01, axp.height*0.5])
    norm=mpl.colors.Normalize(vmin=0, vmax=100)
    cb=mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, label='DC (%)')
    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################

    axpbase = ax.get_position()


    datadir=datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths'
    df = pd.read_table(datadir+'/allstat.dat',
                     header=None, names=['code', 'used', 'NS', 'EW', 'Z', 'f1', 'f2', 'f3', 'f4'], delim_whitespace=True)

    numsta = 0
    axpwidth = 0.3
    axpheight= 0.06
    axpy0 = axpbase.y1-0.06
    axpx0 = axpbase.x1+0.22

    obsnormamp = []
    synnormamp = []

    usedstationindex = []
    for j in range(len(df)):
        if df['used'][j] == 1: usedstationindex.append(j)

    for j in usedstationindex:
        obsdata = np.loadtxt(datadir+'/waveform_fits/'+df['code'][j]+'fil (1).dat')
        syndata = np.loadtxt(datadir+'/waveform_fits/'+df['code'][j]+'syn (1).dat')

        ax1 = fig.add_axes([axpx0, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        ax2 = fig.add_axes([axpx0+axpwidth+0.01, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        ax3 = fig.add_axes([axpx0+(axpwidth+0.01)*2, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        if numsta == 0:
            ax1.text(225, 1.6, 'NS', ha='center', va='bottom')
            ax2.text(225, 1.6, 'EW', ha='center', va='bottom')
            ax3.text(225, 1.6, 'Z', ha='center', va='bottom')

        ax1.set_yticklabels([])
        ax2.set_yticklabels([])
        ax3.set_yticklabels([])

        if j < len(usedstationindex)-1:
            ax1.set_xticklabels([])
            ax2.set_xticklabels([])
            ax3.set_xticklabels([])

        ax1.set_xlim(-15, 450)
        ax2.set_xlim(-15, 450)
        ax3.set_xlim(-15, 450)

        ax1.set_ylim(-1.5, 1.5)
        ax2.set_ylim(-1.5, 1.5)
        ax3.set_ylim(-1.5, 1.5)

        #normamp = max(obsdata[:,1]) * 3
        obsnormamp = np.max([np.max(obsdata[:,1]), np.max(obsdata[:,2]), np.max(obsdata[:,3])])
        synnormamp = np.max([np.max(syndata[:,1]), np.max(syndata[:,2]), np.max(syndata[:,3])])

        if df['NS'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        ax1.plot(obsdata[:,0], obsdata[:,1]/obsnormamp, color=color, lw=1)
        ax1.plot(syndata[:,0], syndata[:,1]/synnormamp, color='r', lw=1, alpha=0.85)

        if df['EW'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        #normamp = max(obsdata[:,2]) * 3
        ax2.plot(obsdata[:,0], obsdata[:,2]/obsnormamp, color=color, lw=1)
        ax2.plot(syndata[:,0], syndata[:,2]/synnormamp, color='r', lw=1, alpha=0.85)

        if df['Z'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        #normamp = max(obsdata[:,3]) * 3
        ax3.plot(obsdata[:,0], obsdata[:,3]/obsnormamp, color=color, lw=1)
        ax3.plot(syndata[:,0], syndata[:,3]/synnormamp, color='r', lw=1, alpha=0.85)

        axp = ax1.get_position()
        fig.text(axp.x0+0.005, axp.y1-0.005, df['code'][j], va='top', ha='left', fontsize=8,
                path_effects=[path_effects.Stroke(linewidth=1 , foreground='w', alpha=1), path_effects.Normal()])

        #print(df['code'][j])

        numsta += 1

    ax1.set_xlabel('Time (s)')

    ax3.plot([], [], color='k', lw=1, label='Obs.')
    ax3.plot([], [], color='C7', lw=1, label='Obs. (unused)')
    ax3.plot([], [], color='r', lw=1, label='Syn.', alpha=0.85)
    ax3.legend(fontsize=6, loc=(0, -1.2), ncol=3, borderaxespad=0., handlelength=1)


    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################

    data = np.loadtxt(datarootdir+'work/NZstation_info.txt', delimiter='|', usecols=(2,3), skiprows=1)
    lon, lat = data[:,1], data[:,0]
    code = np.loadtxt(datarootdir+'work/NZstation_info.txt', delimiter='|', usecols=(1), skiprows=1, dtype=str)

    df_station = pd.DataFrame(data = np.vstack([lon, lat, code]).T,
                       columns=['lon', 'lat', 'code'])

    datadir=datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths'
    df = pd.read_table(datadir+'/allstat.dat',
                     header=None, names=['code', 'used', 'NS', 'EW', 'Z', 'f1', 'f2', 'f3', 'f4'], delim_whitespace=True)

    for n in range(len(lon)):
        if lon[n] < 0:
            lon[n] = lon[n] + 360

    lonmargin, latmargin = 0.1, 0.1
    tickintx,tickinty, tickformat = 5, 5, 0

    lonmin=np.min(lon)-lonmargin; lonmax=np.max(lon)+lonmargin; latmin=np.min(lat)-latmargin; latmax=np.max(lat)+latmargin
    elon, elat = 179.8000000000, -37.3700000000
    lonmin, lonmax, latmin, latmax = elon-5, elon+1, elat-5, elat+9
    #print(lonmin, lonmax, latmin, latmax)

    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
              rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')

    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

    axpwidth = 0.1
    mapheight=axpwidth/aspect

    axpxloc, axpyloc = axpbase.x1-axpwidth+0.18, axpbase.y1-mapheight

    ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.drawcoastlines(color='k', linewidth=0.5, zorder=10)

    for j in range(len(df)):
        if df['used'][j] == 1:
            if df['code'][j] == 'RAO':
                lon, lat = 177.929000, -29.245000
            else:
                lat = float(df_station.loc[(df_station['code'] == df['code'][j])]['lat'].values[0])
                lon = float(df_station.loc[(df_station['code'] == df['code'][j])]['lon'].values[0])
            x, y = m(lon, lat)
            ax.scatter(x, y, marker='^', facecolor='ivory', edgecolor='k', lw=0.75, s=50, zorder=10)
            #ax.text(x, y, df['code'][j])


    x, y = m(179.8000000000, -37.3700000000)
    ax.scatter(x, y, marker='x', lw=2, s=70, color='k')


    ax.set_zorder(11)
    ax2 = utils.mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
    ax2.set_zorder(10)

    plt.savefig(figname, bbox_inches="tight", pad_inches=0.1, dpi=150)
    plt.show()

    
def rcmtHF2D(figname):

    datadir=datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane'
    data = np.loadtxt(datadir+'/sources.gmt', usecols=(0,1,2,3))
    lon, lat, dep = data[:,0][:-2], data[:,1][:-2], data[:,3][:-2]
    for n in range(len(lon)):
        if lon[n] < 0: lon[n] += 360


    data = np.loadtxt(datadir+'/sources_strike_dist.csv', skiprows=1, delimiter=',')
    dist = data[:,5]
    # awk 'NR >= 4 && NR <= 93 {print $0}' < inv1.dat > inv1_extracted_sub1.txt
    # awk 'NR >= 126 && NR <= 215 {print $0}' < inv1.dat > inv1_extracted_sub2.txt

    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[106]
    moment0 = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]

    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[228]
    moment1 = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]


    moments = [moment0, moment1]
    sampling_rate = 0.04
    fig=plt.figure(figsize=figsize)

    cmap = cm.roma_r
    cmap_corr = cm.grayC_r
    bestsol_loc = []
    bestsol_mec = []
    for i in range(2):
        ax = fig.add_axes([0+(0.45)*i, 0, 0.3, 0.4])
        if i == 0:
            axp0 = ax.get_position()
        else:
            axp1 = ax.get_position()
        data = np.loadtxt(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1_extracted_sub'+str(i+1)+'.txt')
        ind, corr, dc, str1, dip1, rake1, ishift = data[:,1], data[:,2], data[:,4], data[:,5], data[:,6], data[:,7],data[:,1]

        ax.axhline(dep[np.argmax(corr)], color='k', linestyle='--', zorder=1, lw=0.75)
        ax.axvline(dist[np.argmax(corr)], color='k', linestyle='--', zorder=1, lw=0.75)

        df = pd.DataFrame(data = np.vstack([corr, dep, dc, str1, dip1, rake1, ishift]).T,
                           columns=['corr', 'dep', 'dc', 'str1', 'dip1', 'rake1', 'ishift'])

        for j in range(len(corr)):
            focmec = [df['str1'][j],df['dip1'][j],df['rake1'][j]]
            x, y = (df['corr'][j], df['dep'][j])

            if corr[j] >= max(corr)*0.9:
                alpha=1
            else:
                alpha=0.2

            if corr[j] == max(corr):
                color_t = 'C5'
                bestsol_loc.append([lon[j], lat[j]])
                bestsol_mec.append(focmec)
            else:
                color_t = 'C7'
            ax.scatter(dist[j], dep[j], zorder=0, s=0)
            tmp = beachball.plot_beachball_mpl(focmec, ax, size=10, position=(dist[j], dep[j]),
                                     beachball_type='dc', edgecolor='k', color_t=color_t,
                                     color_p='w', linewidth=0.75, zorder=2, alpha=alpha)

        #sc=ax.scatter(lat, dep, c=corr, vmin=0, vmax=1, cmap=cm.batlow)

        xi=np.linspace(np.min(dist), np.max(dist), 100)
        yi=np.linspace(np.min(dep), np.max(dep), 100)
        X, Y=np.meshgrid(xi, yi)
        zi=griddata((dist, dep), corr, (X, Y),'linear')
        interval=np.linspace(min(corr), max(corr), 21)
        ax.contourf(X, Y, zi, interval, cmap=cmap_corr, zorder=0)
        ax.contour(X, Y, zi, [max(corr)*0.9, 1], colors='k', zorder=0, linewidths=0.75)



        ax.set_xlim(min(dist)-5, max(dist)+5)
        ax.set_ylim(min(dep)-5, max(dep)+5)
        ax.invert_yaxis()
        ax.invert_xaxis()
        ax.set_xlabel('Distance along strike (km)')
        if i == 0:
            ax.set_ylabel('Depth (km)')

        axp = ax.get_position()
        cax=fig.add_axes([axp.x1+0.005, axp.y0, 0.01, axp.height*0.4])
        norm=mpl.colors.Normalize(vmin=min(corr), vmax=max(corr))
        cb=mpl.colorbar.ColorbarBase(cax, cmap=cmap_corr, norm=norm, label='Correlation')


        timeshift = df.loc[df['corr'].idxmax(), 'ishift'] * sampling_rate
        mw = 2.0 / 3.0 * (np.log10(moments[i]) - 9.1)
        axp = ax.get_position()
        fig.text(axp.x0+0.005, axp.y1+0.005, 'Sub-event '+str(i+1)+'\nMw '+str('{:.1f}'.format(mw))+' | OT+'+str('{:.1f}'.format(timeshift))+' s',
                 va='bottom', fontsize=8)


    ########################################################################
    ########################################################################
    ########################################################################

    lonmargin, latmargin = 0.2, 0.2
    tickintx,tickinty, tickformat = 0.5, 0.5, 1

    lonmin=np.min(lon)-lonmargin; lonmax=np.max(lon)+lonmargin; latmin=np.min(lat)-latmargin; latmax=np.max(lat)+latmargin
    #elon, elat = 179.8000000000, -37.3700000000
    #lonmin, lonmax, latmin, latmax = elon-5, elon+1, elat-5, elat+9
    #print(lonmin, lonmax, latmin, latmax)

    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
              rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')

    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

    axp = ax.get_position()
    #mapheight=axpwidth/aspect
    mapheight = axp.height
    axpwidth = mapheight * aspect


    axpxloc, axpyloc = axp.x1+0.2, axp.y0

    ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
    axplocmap = ax.get_position()
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.drawcoastlines(color='k', linewidth=0.5, zorder=10)

    x, y = m(lon[0:10], lat[0:10])
    ax.scatter(x, y, s=100, facecolor='w', edgecolor='k')

    for i in range(2):
        x, y = m(bestsol_loc[i][0], bestsol_loc[i][1])
        tmp = beachball.plot_beachball_mpl(bestsol_mec[i], ax, size=12, position=(x, y),
                                 beachball_type='dc', edgecolor='k', color_t='C5',
                                 color_p='w', linewidth=0.75, zorder=2, alpha=1)
        ax.text(x-0.1, y, 'Sub-event '+str(i+1), ha='right', va='center', fontsize=6)

    src = shapefile.Reader(datarootdir+'work/tectonicplates/PB2002_boundaries.shp')
    for tmp in src.shapeRecords():
        x = [i[0] for i in tmp.shape.points[:]]
        y = [i[1] for i in tmp.shape.points[:]]
        for n in range(len(x)):
            if x[n] < 0:
                x[n] = x[n] + 360
        x, y = m(x, y)
        ax.plot(x, y, color='C7', lw=0.5, linestyle='--')

    x, y = m(179.875, -37.4)
    ax.text(x, y, 'Trench', ha='left', va='center', fontsize=8, color='gray', rotation=65,
            path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])

    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    model_para = utils.load_fort40(datarootdir+'model_'+str(model[13])+'/fort.40')

    x, y = m(model_para.lon[0], model_para.lat[0])
    sc=ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=100)
    sc.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])

    ax.set_zorder(11)
    ax2 = utils.mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
    ax2.set_zorder(10)

    ######################################################################################################

    datadir=datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane'
    df = pd.read_table(datadir+'/allstat.dat',
                     header=None, names=['code', 'used', 'NS', 'EW', 'Z', 'f1', 'f2', 'f3', 'f4'], delim_whitespace=True)

    numsta = 0
    axpwidth = 0.3
    axpheight= 0.06
    axpy0 = axp0.y0-0.2
    axpx0 = axp0.x0

    obsnormamp = []
    synnormamp = []

    usedstationindex = []
    for j in range(len(df)):
        if df['used'][j] == 1: usedstationindex.append(j)

    for j in usedstationindex:
        obsdata = np.loadtxt(datadir+'/waveform_fits/'+df['code'][j]+'fil.dat')
        syndata = np.loadtxt(datadir+'/waveform_fits/'+df['code'][j]+'syn.dat')

        ax1 = fig.add_axes([axpx0, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        ax2 = fig.add_axes([axpx0+axpwidth+0.01, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        ax3 = fig.add_axes([axpx0+(axpwidth+0.01)*2, axpy0-numsta*(axpheight+0.01), axpwidth, axpheight])
        if numsta == 0:
            ax1.text(225, 1.6, 'NS', ha='center', va='bottom')
            ax2.text(225, 1.6, 'EW', ha='center', va='bottom')
            ax3.text(225, 1.6, 'Z', ha='center', va='bottom')
            axpbase = ax3.get_position()

        ax1.set_yticklabels([])
        ax2.set_yticklabels([])
        ax3.set_yticklabels([])

        if j < len(usedstationindex)-1:
            ax1.set_xticklabels([])
            ax2.set_xticklabels([])
            ax3.set_xticklabels([])

        ax1.set_xlim(-15, 450)
        ax2.set_xlim(-15, 450)
        ax3.set_xlim(-15, 450)

        ax1.set_ylim(-1.5, 1.5)
        ax2.set_ylim(-1.5, 1.5)
        ax3.set_ylim(-1.5, 1.5)

        #normamp = max(obsdata[:,1]) * 3
        obsnormamp = np.max([np.max(obsdata[:,1]), np.max(obsdata[:,2]), np.max(obsdata[:,3])])
        synnormamp = np.max([np.max(syndata[:,1]), np.max(syndata[:,2]), np.max(syndata[:,3])])

        if df['NS'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        ax1.plot(obsdata[:,0], obsdata[:,1]/obsnormamp, color=color, lw=1)
        ax1.plot(syndata[:,0], syndata[:,1]/synnormamp, color='r', lw=1, alpha=0.85)

        if df['EW'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        #normamp = max(obsdata[:,2]) * 3
        ax2.plot(obsdata[:,0], obsdata[:,2]/obsnormamp, color=color, lw=1)
        ax2.plot(syndata[:,0], syndata[:,2]/synnormamp, color='r', lw=1, alpha=0.85)

        if df['Z'][j] <= 0 or df['used'][j] == 0:
            color='C7'
        else:
            color='k'
        #normamp = max(obsdata[:,3]) * 3
        ax3.plot(obsdata[:,0], obsdata[:,3]/obsnormamp, color=color, lw=1)
        ax3.plot(syndata[:,0], syndata[:,3]/synnormamp, color='r', lw=1, alpha=0.85)

        axp = ax1.get_position()
        fig.text(axp.x0+0.005, axp.y1-0.005, df['code'][j], va='top', ha='left', fontsize=8,
                path_effects=[path_effects.Stroke(linewidth=1 , foreground='w', alpha=1), path_effects.Normal()])

        #print(df['code'][j])

        numsta += 1

    ax1.set_xlabel('Time (s)')

    ax3.plot([], [], color='k', lw=1, label='Obs.')
    ax3.plot([], [], color='C7', lw=1, label='Obs. (unused)')
    ax3.plot([], [], color='r', lw=1, label='Syn.', alpha=0.85)
    ax3.legend(fontsize=6, loc=(0, -1.2), ncol=3, borderaxespad=0., handlelength=1)


    #######################################################################
    #######################################################################
    #######################################################################
    #######################################################################

    data = np.loadtxt(datarootdir+'work/NZstation_info.txt', delimiter='|', usecols=(2,3), skiprows=1)
    lon, lat = data[:,1], data[:,0]
    code = np.loadtxt(datarootdir+'work/NZstation_info.txt', delimiter='|', usecols=(1), skiprows=1, dtype=str)

    df_station = pd.DataFrame(data = np.vstack([lon, lat, code]).T,
                       columns=['lon', 'lat', 'code'])

    datadir=datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane'
    df = pd.read_table(datadir+'/allstat.dat',
                     header=None, names=['code', 'used', 'NS', 'EW', 'Z', 'f1', 'f2', 'f3', 'f4'], delim_whitespace=True)

    for n in range(len(lon)):
        if lon[n] < 0:
            lon[n] = lon[n] + 360

    lonmargin, latmargin = 0.1, 0.1
    tickintx,tickinty, tickformat = 5, 5, 0

    lonmin=np.min(lon)-lonmargin; lonmax=np.max(lon)+lonmargin; latmin=np.min(lat)-latmargin; latmax=np.max(lat)+latmargin
    elon, elat = 179.8000000000, -37.3700000000
    lonmin, lonmax, latmin, latmax = elon-5, elon+1, elat-5, elat+9
    #print(lonmin, lonmax, latmin, latmax)

    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,\
              rsphere=(6378137.00,6356752.3142),resolution='i',projection='cyl')

    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))

    axpwidth = 0.15
    mapheight=axpwidth/aspect

    axpxloc, axpyloc = axplocmap.x1-axpwidth, axpbase.y1-mapheight

    ax=fig.add_axes([axpxloc, axpyloc, axpwidth, mapheight])
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.fillcontinents(color='C7', zorder=0, alpha=0.1)
    m.drawcoastlines(color='k', linewidth=0.5, zorder=10)

    for j in range(len(df)):
        if df['used'][j] == 1:
            if df['code'][j] == 'RAO':
                lon, lat = 177.929000, -29.245000
            else:
                lat = float(df_station.loc[(df_station['code'] == df['code'][j])]['lat'].values[0])
                lon = float(df_station.loc[(df_station['code'] == df['code'][j])]['lon'].values[0])
            x, y = m(lon, lat)
            ax.scatter(x, y, marker='^', facecolor='ivory', edgecolor='k', lw=0.75, s=50, zorder=10)
            #ax.text(x, y, df['code'][j])


    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    model_para = utils.load_fort40(datarootdir+'model_'+str(model[13])+'/fort.40')

    x, y = m(model_para.lon[0], model_para.lat[0])
    sc=ax.scatter(x, y, s=100, marker='*', facecolor='none', edgecolor='k', alpha=1, lw=1, zorder=100)


    ax.set_zorder(11)
    ax2 = utils.mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
    ax2.set_zorder(10)

    plt.savefig(figname, bbox_inches="tight", pad_inches=0.1, dpi=150)
    plt.show()

    
def mtcompilation(figname):
    # This study FFM
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    for j in np.arange(13, 14, 1):
        model_para = utils.load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')
        focmec = utils.convertUSEtoNED([0.265294, -1.010916,  0.745621, -0.306763, -0.106606, -0.030296])

    moment = model_para.moment[0]
    mw = 2.0 / 3.0 * (np.log10(moment) - 9.1)
    focmec_ffm = focmec
    mw_ffm = mw

    # USGS
    from pyrocko.io import quakeml
    qml = quakeml.QuakeML.load_xml(filename=datarootdir+'work/quakeml_USGS_Wphase.xml')
    tensor = qml.event_parameters.event_list[0].focal_mechanism_list[0].moment_tensor_list[0].tensor

    mt = pmt.MomentTensor.from_values(
        [tensor.mrr.value, tensor.mtt.value, tensor.mpp.value, tensor.mrt.value, tensor.mrp.value, tensor.mtp.value]
    )
    focmec = mt.m6()
    moment = qml.event_parameters.event_list[0].focal_mechanism_list[0].moment_tensor_list[0].scalar_moment.value
    mw = 2.0 / 3.0 * (np.log10(moment) - 9.1)
    focmec = utils.convertUSEtoNED(focmec)
    focmec_usgs = focmec
    mw_usgs = mw

    # GEOFON
    # https://geofon.gfz-potsdam.de/data/alerts/2021/gfz2021ekhv/mt.txt
    import requests
    import re
    URL = 'https://ftp.webdc.eu/old/data/alerts/2021/gfz2021ekhv/mt.txt'
    response = requests.get(URL, verify=False)
    text = response.content.decode('UTF-8')
    for line in text.splitlines():
        if 'Mrr' in line:
            Mrr = float(line[6:11])
            Mtt = float(line[22:30])
        elif 'Mpp' in line:
            Mpp = float(line[6:11])
            Mrt = float(line[22:30])
        elif 'Mrp' in line:
            Mrp = float(line[6:11])
            Mtp = float(line[22:30])
        elif 'MW' in line:
            mw = float(line.replace('MW', ''))

    #Mrr= 1.74
    #Mtt=-7.84
    #Mpp= 6.10
    #Mrt=-1.13
    #Mrp= 2.37
    #Mtp=-1.28
    focmec = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
    #moment = 7.7*10**19
    #mw = 2.0 / 3.0 * (np.log10(moment) - 9.1)
    focmec = utils.convertUSEtoNED(focmec)
    focmec_geofon = focmec
    mw_geofon = mw

    import obspy
    cat_gcmt = obspy.read_events(datarootdir+'work/SPUD_QUAKEML_bundle.xml')
    for event in cat_gcmt.events:
        if event.resource_id == 'smi:service.iris.edu/fdsnws/event/1/query?eventid=11384594':
            mw = event.magnitudes[0].mag
            focmec=[event.focal_mechanisms[0].moment_tensor.tensor.m_rr,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_tt,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_pp,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_rt,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_rp,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_tp]
    #focmec = [0.364, -1.030, 0.664, 0.036, 0.071, -0.213 ]
    #mw = 7.2
    focmec = utils.convertUSEtoNED(focmec)
    focmec_gcmt = focmec
    mw_gcmt = mw

    # Steve's R-CMT (low freq.)
    #a1, a2, a3, a4, a5, a6 = .574851E+19,  .302531E+20,  .389741E+20,  .880749E+20, -.253843E+20,  .000000E+00
    with open(datarootdir+'R-CMT/LowFreq_SingleSubEvent_LineDepths/inv1.dat') as f:
        data = f.readlines()[47]
    a1, a2, a3, a4, a5, a6 = np.array(data.strip().split()).astype(float)
    Mxx=-a4+a6
    Myy=-a5+a6
    Mzz=a4+a5+a6
    Mxy=a1
    Mxz=a2
    Myz=-a3
    focmec = [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
    focmec_rcmt = focmec
    with open(datarootdir+'R-CMT/LowFreq_SingleSubEvent_LineDepths/inv1.dat') as f:
        data = f.readlines()[51]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]
    mw_rcmt = 2.0 / 3.0 * (np.log10(moment) - 9.1)

    # Steve's R-CMT (high freq., multi-point, 2D)
    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[102]
    a1, a2, a3, a4, a5, a6 = np.array(data.strip().split()).astype(float)
    #a1, a2, a3, a4, a5, a6 = .238361E+20,  .545402E+20,  .558320E+20,  .602305E+20, -.363396E+20, .000000E+00
    Mxx=-a4+a6
    Myy=-a5+a6
    Mzz=a4+a5+a6
    Mxy=a1
    Mxz=a2
    Myz=-a3
    focmec = [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
    focmec_rcmt_1 = focmec
    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[106]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]
    mw_rcmt_1 = 2.0 / 3.0 * (np.log10(moment) - 9.1)

    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[224]
    a1, a2, a3, a4, a5, a6 = np.array(data.strip().split()).astype(float)
    #a1, a2, a3, a4, a5, a6 = -.781991E+19,  .897654E+19, -.446557E+19, -.103638E+20, -.927517E+19,  .000000E+00
    Mxx=-a4+a6
    Myy=-a5+a6
    Mzz=a4+a5+a6
    Mxy=a1
    Mxz=a2
    Myz=-a3
    focmec = [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
    focmec_rcmt_2 = focmec
    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[228]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]
    mw_rcmt_2 = 2.0 / 3.0 * (np.log10(moment) - 9.1)


    # Steve's R-CMT (high freq., multi-point, fixed 1D)
    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths/inv1.dat') as f:
        data = f.readlines()[37]
    a1, a2, a3, a4, a5, a6 = np.array(data.strip().split()).astype(float)
    #a1, a2, a3, a4, a5, a6 = .140224E+20,  .233149E+20,  .459788E+20,  .105113E+21, -.376678E+20,  .000000E+00,
    Mxx=-a4+a6
    Myy=-a5+a6
    Mzz=a4+a5+a6
    Mxy=a1
    Mxz=a2
    Myz=-a3
    focmec = [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
    focmec_rcmt_1d_1 = focmec
    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths/inv1.dat') as f:
        data = f.readlines()[41]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]
    mw_rcmt_1d_1 = 2.0 / 3.0 * (np.log10(moment) - 9.1)

    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths/inv1.dat') as f:
        data = f.readlines()[94]
    a1, a2, a3, a4, a5, a6 = np.array(data.strip().split()).astype(float)
    #a1, a2, a3, a4, a5, a6 = -.990366E+19, .671121E+19, -.656904E+19, -.620102E+19, -.941617E+19,  .000000E+00,
    Mxx=-a4+a6
    Myy=-a5+a6
    Mzz=a4+a5+a6
    Mxy=a1
    Mxz=a2
    Myz=-a3
    focmec = [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
    focmec_rcmt_1d_2 = focmec
    with open(datarootdir+'R-CMT/HighFreq_MultipleSubEvents_LineDepths/inv1.dat') as f:
        data = f.readlines()[98]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]
    mw_rcmt_1d_2 = 2.0 / 3.0 * (np.log10(moment) - 9.1)

    # SCARDEC (GEOSCOPE) http://geoscope.ipgp.fr/index.php/en/catalog/earthquake-description?seis=at00qpg5dz
    focmec_scardec = [289, 45, 113]
    mw_scardec = 7.41


    # Regional CMT by GeoNet
    url = 'https://raw.githubusercontent.com/GeoNet/data/main/moment-tensor/GeoNet_CMT_solutions.csv'
    dfGeoNetCMTData = pd.read_csv(url)
    dfGeoNetCMTData.Date = pd.to_datetime(dfGeoNetCMTData.Date, format='%Y%m%d%H%M%S')
    dfGeoNetCMTData.Longitude[dfGeoNetCMTData.Longitude < 0] += 360

    tmp = dfGeoNetCMTData.loc[(dfGeoNetCMTData['PublicID'] == '2021p169083')]
    focmec_geonet = [tmp.Mxx.values[0],tmp.Myy.values[0],tmp.Mzz.values[0],tmp.Mxy.values[0],tmp.Mxz.values[0],tmp.Myz.values[0]]
    mw_geonet = tmp.Mw.values[0]

    # Steve's W-phase
    from pyrocko.model import load_events
    mean = load_events("../materials/WphaseCMT/mean_wphase.yaml")[0]  # Read in event file containing best solution
    focmec_wphase = mean.moment_tensor
    mw_wphase = 2.0 / 3.0 * (np.log10(mean.moment_tensor.moment) - 9.1)

    # list
    focmecs = [focmec_ffm, focmec_wphase, focmec_rcmt, focmec_rcmt_1d_1, focmec_rcmt_1d_2, focmec_rcmt_1, focmec_rcmt_2,
               focmec_gcmt, focmec_usgs, focmec_geofon, focmec_scardec, focmec_geonet]
    mws = [mw_ffm, mw_wphase, mw_rcmt, mw_rcmt_1d_1, mw_rcmt_1d_2, mw_rcmt_1, mw_rcmt_2,
           mw_gcmt, mw_usgs, mw_geofon, mw_scardec, mw_geonet]
    labels = ['This study\nFFM', 'This study\nW-phase','This study\nR-CMT LF',
              'This study\nR-CMT HF1', 'This study\nR-CMT HF2', 'This study\nR-CMT 2D HF1', 'This study\nR-CMT 2D HF2',
              'GCMT', 'USGS', 'GEOFON', 'SCARDEC', 'GeoNet']
    colors = ['C5', 'C0', 'C1', 'C2', 'C2', 'C3', 'C3',
              'C7', 'C7', 'C7', 'C7', 'C7']

    fig=plt.figure(figsize=figsize)
    ax=fig.add_axes([0.1, 0.1, 1.25, 0.15])

    x = -0.1
    for focmec, mw, label, color in zip(focmecs, mws, labels, colors):
        x += 0.3
        beachball.plot_beachball_mpl(focmec, ax, size=3*mw, position=(x, 0.35),
                                     beachball_type='deviatoric', linewidth=1, edgecolor=color, color_t=color)
        beachball.plot_beachball_mpl(focmec, ax, size=3*mw, position=(x, 0.35),
                                     beachball_type='dc', linewidth=1, color_t='none', color_p='none')
        ax.text(x, 0.68, label, ha='center', va='center', fontsize=6)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(0, 3.67)
    plt.savefig(figname, bbox_inches='tight', pad_inches=0.1)
    plt.show()

    
def bathy3D(figname):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca(projection='3d', alpha=1)
    locp = []
    basep = []
    for i in np.arange(0, 101, 1):
        data = np.loadtxt(datarootdir+'work/bathymetry/Xsection_'+str(i+1)+'.txt')
        x, y, z  = data[:,2], 51-np.ones(len(data))*i*0.5, data[:,3]
        xlist = [ x[loc] for loc in np.arange(0, len(x), 2) ]
        ylist = [ y[loc] for loc in np.arange(0, len(x), 2) ]
        zlist = [ z[loc] for loc in np.arange(0, len(x), 2) ]
        locp.append([xlist, ylist, zlist])
        tmpx = xlist[ np.argmin(zlist) ]
        tmpy = ylist[ np.argmin(zlist) ]
        tmpz = zlist[ np.argmin(zlist) ]
        basep.append([tmpx, tmpy, tmpz])
        if i > 0: # for frames
            for j in [len(xlist)-1]: # for background frame
                ax.plot([xlist[j], locp[i-1][0][j]], [ylist[j], locp[i-1][1][j]], [zlist[j], locp[i-1][2][j]],
                        color='k', lw=1, solid_capstyle='round')

            for j in [0]: # for foreground frame
                ax.plot([xlist[j], locp[i-1][0][j]], [ylist[j], locp[i-1][1][j]], [zlist[j], locp[i-1][2][j]],
                        color='k', lw=1, solid_capstyle='round')

    locp = np.array(locp)
    ax.plot_surface(locp[:,0], locp[:,1], locp[:,2], cmap=cm.lapaz, vmin=-6200, vmax=0, alpha=0.75, rcount=201, ccount=201)
    basep = np.array(basep)
    ax.plot(locp[-1][0],locp[-1][1],locp[-1][2], color='k', lw=1, solid_capstyle='round', solid_joinstyle='round') # for background frame
    ax.plot(locp[0][0],locp[0][1],locp[0][2], color='k', lw=1, solid_capstyle='round', solid_joinstyle='round') # for foreground frame
    ax.view_init(azim=123, elev=21)
    ax.set_zlim(-20000, 0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_ylim(-78, 50)
    ax.set_axis_off()
    plt.savefig(figname, transparent=True, pad_inches=0, bbox_inches="tight", dpi=300)
    plt.show()
