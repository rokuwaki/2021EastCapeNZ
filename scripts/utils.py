import numpy as np
from netCDF4 import Dataset
from cmcrameri import cm
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.patheffects as path_effects
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import obspy
from pyrocko.plot import beachball
from pyrocko import moment_tensor as pmt
from pyrocko.model import load_events
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84
import pandas as pd
import shapefile

# some default font setting using Open Sans https://fonts.google.com/specimen/Open+Sans
initfontsize = 10
mpl.rc('axes', labelsize=initfontsize, titlesize=initfontsize)
mpl.rc('xtick', labelsize=initfontsize)
mpl.rc('ytick', labelsize=initfontsize)
mpl.rc('legend', fontsize=initfontsize, edgecolor='none')
mpl.rc('savefig', dpi=600, transparent=False)
mpl.rc('font', size=initfontsize)

mpl.rcParams['font.weight'] = 400
mpl.rcParams['font.family'] = 'Open Sans'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Open Sans'
mpl.rcParams['mathtext.it'] = 'Open Sans:italic'
mpl.rcParams['mathtext.bf'] = 'Open Sans:bold'

# default colour cycle ('C0', 'C1', ..., 'C7') using `Set2`
cmap = plt.get_cmap('Set2', 8)
from cycler import cycler
custom_color_cycle=[]
for i in range(cmap.N):
    rgb = cmap(i)[:3]
    custom_color_cycle.append(str(mpl.colors.rgb2hex(rgb)))
plt.rc('axes', prop_cycle=(cycler(color=custom_color_cycle)))


def drawRCMT(fig, ax, m, datarootdir, axes, model_para):
    dx = model_para.xx[0]
    meclocs, mos, focmecs, labels = loadmecs(datarootdir)
    x, y = m(179.62, -37.23)
    if axes == 0:
        i = 5
        cmtdepth = 12
        cmaprcmt = cm.nuuk
        color_t = mpl.colors.rgb2hex(cmaprcmt(2/5))
    else:
        i = 4
        cmtdepth = 52
        cmaprcmt = cm.nuuk
        color_t = mpl.colors.rgb2hex(cmaprcmt(1/5))
    tmp = beachball.plot_beachball_mpl(focmecs[i], ax, size=dx*3, position=(x, y),
                             beachball_type='dc', edgecolor=color_t, color_t=color_t,
                             color_p='w', linewidth=0.5, alpha=0.5, zorder=100)
    tmp = beachball.plot_beachball_mpl(focmecs[i], ax, size=dx*3, position=(x, y),
                             beachball_type='dc', edgecolor='k', color_t='none',
                             color_p='none', linewidth=0.5, alpha=1, zorder=100)
    ax.text(x, y-0.025, 'R-CMT\n'+r'$z=$'+str(cmtdepth)+' km', ha='center', va='top', fontsize=8, color='k',
            path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])
    
    axp = ax.get_position()
    axinset=fig.add_axes([axp.x0, axp.y1-0.1, 0.1, 0.1])
    axinset.set_alpha(0)
    axinset.set_facecolor('none')
    axinset.set_xticks([])
    axinset.set_yticks([])
    
    


def drawPaxis(fig, ax, ax2, m, axes, datarootdir, modelid, model_para):
    cmap = cm.bilbao
    
    data=np.loadtxt(datarootdir+'model_'+str(modelid)+'/FFM_MT.txt')
    m1,m2,m3,m4,m5,m6 = data[:,4],data[:,5],data[:,6],data[:,7],data[:,8],data[:,9]
        
    data=np.loadtxt(datarootdir+'model_'+str(modelid)+'/FFM_DCall.txt', skiprows=1)
    lon, lat, depth, slip = data[:,1], data[:,2], data[:,3], data[:,4]
    stk0, dip0, rake0 = data[:,5], data[:,6],data[:,7]
    stk1, dip1, rake1 = data[:,8], data[:,9],data[:,10]
    strike, dip, rake = [],[],[]
    for i in range(len(stk0)):
        tmpstr, tmpdip, tmprake = selectplane(model_para.strike[0], model_para.dip[0], stk0[i], dip0[i], rake0[i], stk1[i], dip1[i], rake1[i])
        strike.append(tmpstr)
        dip.append(tmpdip)
        rake.append(tmprake)

    if axes == 0:
        dep0, dep1 = 0, 50
        axp = ax.get_position()
        fig.text(axp.x1-0.01, axp.y0+0.01, 'P-axes: Shallow layers\n(0–50 km depth)', ha='right', va='bottom', fontsize=8,
                path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])
    else:
        dep0, dep1 = 50, 1000
        axp = ax.get_position()
        fig.text(axp.x1-0.01, axp.y0+0.01, 'P-axes: Deep layers\n(>50 km depth)', ha='right', va='bottom', fontsize=8,
                path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=1), path_effects.Normal()])

    for i in range(len(strike)):
        if depth[i] > dep0 and depth[i] <= dep1:
            if lon[i] < 0:
                lon[i] = lon[i] + 360
            x, y = m(lon[i], lat[i])
            focmecs=[strike[i], dip[i], rake[i]]
            color=cmap(slip[i]/max(slip))

            dx = model_para.xx[0]
            focmecs=[m1[i],m2[i],m3[i],m4[i],m5[i],m6[i]]
            focmecs = convertUSEtoNED(focmecs)
            tmp = beachball.plot_beachball_mpl(focmecs, ax, size=dx*3, position=(x, y),
                                     beachball_type='dc', edgecolor='k', color_t=cmap(slip[i]/max(slip)),
                                     color_p='w', linewidth=0.5, alpha=slip[i]/max(slip), zorder=int(slip[i]/max(slip)*100))

            mt = pmt.MomentTensor(
                strike=strike[i],
                dip=dip[i],
                rake=rake[i])
            pnt = get_plunge_azimuth(mt.m6_up_south_east())
            t_azi = pnt[5] # T-axis azimuth
            p_azi = pnt[1] # p-axis azimuth
            for pm in [-1, 1]:
                #tmp = geod.Direct(lat[i], lon[i], t_azi, pm*15*1e3*slip[i]/max(slip))
                tmp = geod.Direct(lat[i], lon[i], p_azi, pm*20*1e3*slip[i]/max(slip))
                if tmp['lon1'] < 0: tmp['lon1'] += 360
                if tmp['lon2'] < 0: tmp['lon2'] += 360
                x0, y0 = m(tmp['lon1'], tmp['lat1'])
                x1, y1 = m(tmp['lon2'], tmp['lat2'])
                ax.plot([x0,x1], [y0,y1], color=cmap(slip[i]/max(slip)), lw=1,
                        zorder=int(slip[i]/max(slip)*100), alpha=slip[i]/max(slip))
                
    if axes == 1:
        ax2.set_yticklabels([])

        axp = ax.get_position()
        cax=fig.add_axes([axp.x1+0.005, axp.y0, 0.01, axp.height/2])
        norm=mpl.colors.Normalize(vmin=0, vmax=max(slip))
        cb=mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, label='Slip (m)',
                                     ticks=np.linspace(0, max(slip), 5), format='%.2f', alpha=1)
    
    


def drawslip3D(datarootdir, model_para, modelid, fig, ax, axlabel, lon, lat, depth, slip):
    
    slipvertslist = getSlipCell3D(datarootdir, model_para, modelid, ax, 1)
    cmap = cm.bilbao
    for i,verts in enumerate(slipvertslist):
        pc = Poly3DCollection(verts)
        pc.set_alpha(0.85)
        pc.set_facecolor(cmap(slip[i]/max(slip)))
        ax.add_collection3d(pc)
    verts = plotModelEdge3D(model_para, ax, 500)
    vertslist = []
    vertslist.append(verts)
    
    tmplon, tmplat = np.array(vertslist[0][0])[:,0], np.array(vertslist[0][0])[:,1]
    ax.plot([max(tmplon), min(tmplon)], [max(tmplat), min(tmplat)], [50, 50], zorder=1000, color='k', lw=0.75, linestyle='--', alpha=0.3)
    for lonlabel, latlabel, deplabel, text in zip([max(tmplon), max(tmplon), max(tmplon)],
                                                  [max(tmplat)+0.017, max(tmplat)+0.017, max(tmplat)+0.017],
                                                  [50, 20, 78.7],
                                                  ['50 km', 'Shallow', 'Deep']):    
        ax.text(lonlabel, latlabel, deplabel, text, ha='right', va='center', fontsize=8, color='k',
                    path_effects=[path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])
    sc=ax.scatter(model_para.lon[0]-0.1, model_para.lat[0]-0.014, model_para.depth[0]-3.4, marker='*', s=100, edgecolor='k', facecolor='w', zorder=1000)
    sc.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])

    axp = axlabel.get_position()
    cax=fig.add_axes([axp.x1-0.15, axp.y0+0.3, 0.01, axp.height*0.4])
    norm=mpl.colors.Normalize(vmin=0, vmax=max(slip))
    cb=mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, label='Slip (m)',
                                 ticks=np.linspace(0, max(slip), 5), format='%.2f', alpha=0.85)
    


def drawbaseframe3D(model_para, fig, ax, azim, elev):
    lonmin=model_para.lon[0]-1*0.85
    lonmax=model_para.lon[0]+1*0.85
    latmin=model_para.lat[0]-1*0.5
    latmax=model_para.lat[0]+0.75*0.5

    ax.view_init(azim=azim, elev=elev)
    ax.set_zlim(0., 110)
    ax.set_xlim(lonmin, lonmax)
    ax.set_ylim(latmin, latmax)
    ax.invert_zaxis()
    ax.xaxis.pane.set_facecolor('none')
    ax.yaxis.pane.set_facecolor('none')
    ax.zaxis.pane.set_facecolor('none')
    ax.grid(False)
    ax.set_axis_off()

    tmpx = np.array([lonmin, lonmax])
    tmpy = np.array([latmax, latmax])
    tmpz = np.array([[0, 0],[115, 115]])
    ax.plot_surface(tmpx, tmpy, tmpz, color='whitesmoke', alpha=0.1)

    tmpx = np.array([lonmax, lonmax])
    tmpy = np.array([latmin, latmax])
    tmpz = np.array([[0, 0],[115, 115]])
    ax.plot_surface(tmpx, tmpy, tmpz, color='whitesmoke', alpha=0.1)

    Xs = np.linspace(lonmin, lonmax, 2)
    Ys = np.linspace(latmin, latmax, 2)
    Xs, Ys = np.meshgrid(Xs, Ys)
    Zs = np.array([[115,115],[115,115]])
    ax.plot_surface(Xs, Ys, Zs, color='whitesmoke', alpha=0.1)
    
    ax = fig.gca(projection='3d', alpha=1)

    ax.plot([lonmin, lonmax], [latmin, latmin], [0, 0], color='k', lw=0.75)
    ax.plot([lonmin, lonmax], [latmax, latmax], [0, 0], zorder=1000, color='k', lw=0.75)
    ax.plot([lonmin, lonmin], [latmin, latmax], [0, 0], zorder=1000, color='k', lw=0.75)
    ax.plot([lonmax, lonmax], [latmin, latmax], [0, 0], zorder=1000, color='k', lw=0.75)

    ax.plot([lonmin, lonmax], [latmin, latmin], [115, 115], color='k', lw=0.75)
    ax.plot([lonmin, lonmax], [latmax, latmax], [115, 115], zorder=1000, color='k', lw=0.75)
    ax.plot([lonmin, lonmin], [latmin, latmax], [115, 115], zorder=1000, color='k', lw=0.75)
    ax.plot([lonmax, lonmax], [latmin, latmax], [115, 115], zorder=1, color='k', lw=0.75)

    ax.plot([lonmin, lonmin], [latmin, latmin], [0, 115], color='k', lw=0.75)
    ax.plot([lonmax, lonmax], [latmin, latmin], [0, 115], color='k', lw=0.75)
    ax.plot([lonmin, lonmin], [latmax, latmax], [0, 115], color='k', lw=0.75)
    ax.plot([lonmax, lonmax], [latmax, latmax], [0, 115], color='k', lw=0.75)
    
    
    xtickint = 0.5
    ytickint = 0.3
    xticklabels = np.arange(lonmin - lonmin%xtickint+xtickint, lonmax - lonmax%xtickint + xtickint, xtickint)
    yticklabels = np.arange(latmin - latmin%ytickint+ytickint, latmax - latmax%ytickint + ytickint, ytickint)[:-1]
    zticklabels = np.arange(0, 120, 20)

    for xtick in xticklabels:
        ax.plot([xtick, xtick], [latmin-0.01, latmin], [115, 115], zorder=1, color='k', lw=0.75)
        ax.text(xtick, latmin-0.03, 115, str('{:.1f}'.format(xtick)), ha='left', va='center')
    for ytick in yticklabels:
        ax.plot([lonmin-0.05, lonmin], [ytick, ytick], [115, 115], zorder=1, color='k', lw=0.75)
        ax.text(lonmin-0.25, ytick, 115, str('{:.1f}'.format(ytick)), ha='center', va='center')
    for ztick in zticklabels:
        ax.plot([lonmin-0.05, lonmin], [latmax, latmax], [ztick, ztick], zorder=1, color='k', lw=0.75)
        ax.text(lonmin-0.15, latmax, ztick, str(ztick), ha='right', va='center')

    ax.view_init(azim=azim, elev=elev)
    ax.set_zlim(0., 110)
    ax.set_xlim(lonmin, lonmax)
    ax.set_ylim(latmin, latmax)
    ax.invert_zaxis()
    ax.xaxis.pane.set_facecolor('none')
    ax.yaxis.pane.set_facecolor('none')
    ax.zaxis.pane.set_facecolor('none')
    ax.grid(False)
    ax.set_axis_off()

    axlabel = fig.add_axes([0.1, 0.1, 0.85, 0.85])
    axlabel.set_facecolor('none')
    axlabel.text(0.4, 0.08, 'Latitude ($\degree$)', ha='center', va='top', rotation=-5)
    axlabel.text(0.8, 0.1, 'Longitude ($\degree$)', ha='left', va='bottom', rotation=65)
    axlabel.text(0.12, 0.4, 'Depth (km)', ha='right', va='center', rotation=92)
    
    x0, y0 = 0.33, 0.18
    theta, scale = 176.5, 0.1
    x1, y1 = x0+np.cos(np.deg2rad(theta))*scale, y0+np.sin(np.deg2rad(theta))*scale
    axlabel.annotate('N', xytext=(x1,y1), xy=(x0, y0), arrowprops=dict(arrowstyle="<-", linewidth=0.75, color='k'))
    axlabel.set_axis_off()

    ax.set(facecolor='none')
    
    return ax, axlabel


def drawlabels(elat, elon, ax, m):
    x, y=m(elon-0.1, elat)
    text=ax.text(x, y, r'$M_{\rm{W}}$ 7.3 2021-03-04', fontsize=8, ha='right', va='center', zorder=100,)
    text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='w', alpha=1), path_effects.Normal()])

    x, y=m(179.7, -37.9)
    text=ax.text(x, y, 'Ruatoria\nIndentation', fontsize=8, ha='left', va='top', zorder=100,)
    text.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])

    x, y=m(179.4, -38.6)
    text=ax.text(x, y, 'Hikurangi margin', fontsize=8, ha='left', va='top', zorder=100,)
    text.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])

    x, y=m(180.35, -37.1)
    text=ax.text(x, y, 'Kermadec trench', fontsize=8, ha='center', va='center', zorder=100, rotation=57)
    text.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])
    
    plon, plat, azi, vel = [180.6, 180.4, 180.2], [-37.8, -38.37, -38.9], [262.61, 262.17, 261.72], [46.40, 45.38, 44.37]
    for i in range(len(plat)):
        g = geod.Direct(plat[i], plon[i], azi[i], vel[i]/2*1e3)
        if g['lon2'] < 0: g['lon2'] += 360
        spoint = m(plon[i], plat[i])
        epoint = m(g['lon2'], g['lat2'])
        an = ax.annotate('', xy=(epoint[0], epoint[1]), xytext=(spoint[0], spoint[1]),
                         arrowprops=dict(arrowstyle="simple,head_length=0.3,head_width=0.3,tail_width=0.1", edgecolor='w', facecolor='k', lw=0.2))

        x, y = m(g['lon2']+0.05, g['lat2']-0.02)
        text=ax.text(x, y, '{:.1f}'.format(vel[i])+' mm/yr', va='top', ha='left', size=6, color='k',
                    path_effects=[path_effects.Stroke(linewidth=1, foreground='w', alpha=0.75), path_effects.Normal()])
    


def drawbeachball(data,datarootdir,lonmin,lonmax,latmin,latmax,ax,m):
    cat_gcmt = obspy.read_events(data)
    for event in cat_gcmt.events:
        mag = event.magnitudes[0].mag
        lon, lat, dep = event.origins[1].longitude, event.origins[1].latitude, event.origins[1].depth
        if lon < 0: lon += 360
        if lon > lonmin+0.01 and lon < lonmax-0.01 and lat > latmin+0.01 and lat < latmax-0.01 and event.origins[0].time <= obspy.UTCDateTime('2021-03-04'):
            focmecs=[event.focal_mechanisms[0].moment_tensor.tensor.m_rr,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_tt,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_pp,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_rt,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_rp,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_tp]
            focmecs = convertUSEtoNED(focmecs)

            if mag >= 7:
                size = 16
            elif mag >= 5 and mag < 7:
                size = 8
            else:
                size = 4

            x, y = m(lon, lat)
            tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                     beachball_type='deviatoric', edgecolor='none', color_t='C7',
                                     color_p='w', linewidth=0.5, alpha=1, zorder=int(100-mag*10))
            tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                     beachball_type='dc', edgecolor='k', color_t='none',
                                     color_p='none', linewidth=0.5, alpha=1, zorder=int(100-mag*10))

            if event.resource_id == 'smi:service.iris.edu/fdsnws/event/1/query?eventid=11384594': # mainshock
                tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                         beachball_type='deviatoric', edgecolor='none', color_t='C5',
                                         color_p='w', linewidth=0.5, alpha=1, zorder=int(100))
                tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                         beachball_type='dc', edgecolor='k', color_t='none',
                                         color_p='none', linewidth=0.5, alpha=1, zorder=int(100))
                
    loclist = ([lonmin+0.18, latmax-0.18], [lonmin+0.18, latmax-0.195], [lonmin+0.18, latmax-0.255])
    sizelist = [16, 8, 4]
    for loc, size in zip(loclist, sizelist):
        x, y = m(loc[0]-0.035, loc[1]+0.035)
        focmecs = [0.265296, -1.010935,  0.745638, -0.306760, -0.106615, -0.030290]
        focmecs = convertUSEtoNED(focmecs)
        tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                 beachball_type='deviatoric', edgecolor='none', color_t='C7',
                                 color_p='w', linewidth=0.5, alpha=1, zorder=int(111))
        tmp = beachball.plot_beachball_mpl(focmecs, ax, size=size, position=(x, y),
                                 beachball_type='dc', edgecolor='k', color_t='none',
                                 color_p='none', linewidth=0.5, alpha=1, zorder=int(111))
    ax.scatter([], [], label=r'$M_{\rm{W}} \geq 7$', s=0)
    ax.scatter([], [], label=r'$M_{\rm{W}} \geq 5$', s=0)
    ax.scatter([], [], label=r'$M_{\rm{W}} \geq 3$', s=0)
    ax.legend(loc='upper left', fontsize=6, labelspacing=0.01, handletextpad=1, borderpad=0.5).set_zorder(110)

    meclocs, mos, focmecs, labels = loadmecs(datarootdir)
    meclocs_fix_lon = np.linspace(179.3, 180.25, 4)
    meclocs_fix_lat = np.ones(4) * -36.45
    colors = ['C5', 'C5', 'C5', 'C7']
    for i in range(len(colors)):
        x0, y0 = m(meclocs[i][0], meclocs[i][1])
        ax.scatter(x0, y0, facecolor='k', edgecolor='none', s=5, zorder=1000)
        x1, y1 = m(meclocs_fix_lon[i], meclocs_fix_lat[i])
        ax.plot([x0,x1], [y0,y1], color='k', lw=0.5, zorder=10)
        size=16
        tmp = beachball.plot_beachball_mpl(focmecs[i], ax, size=size, position=(x1, y1),
                                 beachball_type='deviatoric', edgecolor='none', color_t=colors[i],
                                 color_p='w', linewidth=0.5, alpha=1, zorder=int(111))
        tmp = beachball.plot_beachball_mpl(focmecs[i], ax, size=size, position=(x1, y1),
                                 beachball_type='dc', edgecolor='k', color_t='none',
                                 color_p='none', linewidth=0.5, alpha=1, zorder=int(111))
        x, y = m(meclocs_fix_lon[i], meclocs_fix_lat[i]+0.16)
        ax.text(x, y-0.02, labels[i], fontsize=6, ha='center', va='center', zorder=100,
                     path_effects=[path_effects.Stroke(linewidth=1.5, foreground='w', alpha=1), path_effects.Normal()])
    


def drawbasemap(fig,lonmin,lonmax,latmin,latmax,axpxloc,axpyloc,axpwidth,tickintx,tickinty,tickformat,coastflag=0,resolution='i',continentalpha=0):
    m=Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax,rsphere=(6378137.00,6356752.3142),resolution=resolution,projection='cyl')
    x, y=m([lonmin, lonmax], [latmin, latmax])
    aspect=abs(max(x)-min(x))/abs(max(y)-min(y))
    axpheight=axpwidth/aspect
    ax=fig.add_axes([axpxloc, axpyloc, axpwidth, axpheight])
    m.fillcontinents(color='C7', zorder=0, alpha=continentalpha)
    if coastflag == 1:
        m.drawcoastlines(color='k', linewidth=0.5, zorder=0)
    ax2 = mapTicksBasemap(fig,m,ax,tickintx,tickinty,lonmin,lonmax,latmin,latmax,tickformat)
    return m, ax, ax2


# mpl-ish tick style for Basemap
def mapTicksBasemap(fig,m,ax,xtickint,ytickint,lonmin,lonmax,latmin,latmax, dpoint):
    axp=ax.get_position()
    ax2 = fig.add_axes([axp.x0, axp.y0, axp.width, axp.height])

    #xticklabels = np.arange(int(lonmin - lonmin%xtickint), int(lonmax - lonmax%xtickint) + xtickint, xtickint)
    #yticklabels = np.arange(int(latmin - latmin%ytickint), int(latmax - latmax%ytickint) + ytickint, ytickint)

    xticklabels = np.arange(lonmin - lonmin%xtickint, lonmax - lonmax%xtickint + xtickint, xtickint)
    yticklabels = np.arange(latmin - latmin%ytickint, latmax - latmax%ytickint + ytickint, ytickint)

    tmp = [m(i, latmin) for i in xticklabels]
    xticks = [tmp[i][0] for i in range(len(xticklabels))]
    ax2.set_xticks(xticks)

    tmp = [m(lonmin, i) for i in yticklabels]
    yticks = [tmp[i][1] for i in range(len(yticklabels))]
    ax2.set_yticks(yticks)

    #print(xticklabels)
    #tmp = [ str(xticklabels[i])+'$\degree$' for i in range(len(xticklabels)) ]
    if dpoint == 0:
        tmp = [ str('{:.0f}'.format(xticklabels[i]))+'$\degree$E' if xticklabels[i] > 0 and xticklabels[i] <= 180 else \
               str('{:.0f}'.format(abs(xticklabels[i])))+'$\degree$W' if  xticklabels[i] < 0 else \
               str('{:.0f}'.format(abs(360-xticklabels[i])))+'$\degree$W' if  xticklabels[i] > 180 else \
               str('{:.0f}'.format(xticklabels[i]))+'$\degree$' for i in range(len(xticklabels)) ]
        ax2.set_xticklabels(tmp)

        #tmp = [ str(yticklabels[i])+'$\degree$' for i in range(len(yticklabels)) ]
        tmp = [ str('{:.0f}'.format(yticklabels[i]))+'$\degree$N' if yticklabels[i] > 0 else \
               str('{:.0f}'.format(abs(yticklabels[i])))+'$\degree$S' if  yticklabels[i] < 0 else \
               str('{:.0f}'.format(yticklabels[i]))+'$\degree$' for i in range(len(yticklabels)) ]
        ax2.set_yticklabels(tmp)

    elif dpoint == 1:
        tmp = [ str('{:.1f}'.format(xticklabels[i]))+'$\degree$E' if xticklabels[i] > 0 and xticklabels[i] <= 180 else \
               str('{:.1f}'.format(abs(xticklabels[i])))+'$\degree$W' if  xticklabels[i] < 0 else \
               str('{:.1f}'.format(abs(360-xticklabels[i])))+'$\degree$W' if  xticklabels[i] > 180 else \
               str('{:.1f}'.format(xticklabels[i]))+'$\degree$' for i in range(len(xticklabels)) ]
        ax2.set_xticklabels(tmp)

        #tmp = [ str(yticklabels[i])+'$\degree$' for i in range(len(yticklabels)) ]
        tmp = [ str('{:.1f}'.format(yticklabels[i]))+'$\degree$N' if yticklabels[i] > 0 else \
               str('{:.1f}'.format(abs(yticklabels[i])))+'$\degree$S' if  yticklabels[i] < 0 else \
               str('{:.1f}'.format(yticklabels[i]))+'$\degree$' for i in range(len(yticklabels)) ]
        ax2.set_yticklabels(tmp)

    elif dpoint == 2:
        tmp = [ str('{:.2f}'.format(xticklabels[i]))+'$\degree$E' if xticklabels[i] > 0 and xticklabels[i] <= 180 else \
               str('{:.2f}'.format(abs(xticklabels[i])))+'$\degree$W' if  xticklabels[i] < 0 else \
               str('{:.2f}'.format(abs(360-xticklabels[i])))+'$\degree$W' if  xticklabels[i] > 180 else \
               str('{:.2f}'.format(xticklabels[i]))+'$\degree$' for i in range(len(xticklabels)) ]
        ax2.set_xticklabels(tmp)

        #tmp = [ str(yticklabels[i])+'$\degree$' for i in range(len(yticklabels)) ]
        tmp = [ str('{:.2f}'.format(yticklabels[i]))+'$\degree$N' if yticklabels[i] > 0 else \
               str('{:.2f}'.format(abs(yticklabels[i])))+'$\degree$S' if  yticklabels[i] < 0 else \
               str('{:.2f}'.format(yticklabels[i]))+'$\degree$' for i in range(len(yticklabels)) ]
        ax2.set_yticklabels(tmp)

    elif dpoint == 3:
        tmp = [ str('{:.3f}'.format(xticklabels[i]))+'$\degree$E' if xticklabels[i] > 0 and xticklabels[i] <= 180 else \
               str('{:.3f}'.format(abs(xticklabels[i])))+'$\degree$W' if  xticklabels[i] < 0 else \
               str('{:.3f}'.format(abs(360-xticklabels[i])))+'$\degree$W' if  xticklabels[i] > 180 else \
               str('{:.3f}'.format(xticklabels[i]))+'$\degree$' for i in range(len(xticklabels)) ]
        ax2.set_xticklabels(tmp)

        #tmp = [ str(yticklabels[i])+'$\degree$' for i in range(len(yticklabels)) ]
        tmp = [ str('{:.3f}'.format(yticklabels[i]))+'$\degree$N' if yticklabels[i] > 0 else \
               str('{:.3f}'.format(abs(yticklabels[i])))+'$\degree$S' if  yticklabels[i] < 0 else \
               str('{:.3f}'.format(yticklabels[i]))+'$\degree$' for i in range(len(yticklabels)) ]
        ax2.set_yticklabels(tmp)

    xlimmin, ylimmin = m(lonmin, latmin)
    xlimmax, ylimmax = m(lonmax, latmax)
    ax2.set_xlim(xlimmin, xlimmax)
    ax2.set_ylim(ylimmin, ylimmax)
    minZorder=min([_.zorder for _ in ax.get_children()])
    ax2.set_zorder(minZorder-1)
    ax2.set_facecolor('none')
    #ax.set_xticks([]); ax.set_yticks([])

    return ax2


def drawbathy(data, fig, ax, m):
    fh = Dataset(data, mode='r')
    lons = fh.variables['lon'][:]; lats = fh.variables['lat'][:]; tmax = fh.variables['z'][:]
    fh.close()
    vmin, vmax = -6200, 0
    ls = LightSource(azdeg=225, altdeg=75)
    rgb = ls.shade(tmax, cmap=cm.grayC_r, vmin=vmin, vmax=vmax)
    im = ax.imshow(rgb, origin='lower',alpha=1, zorder=-1, interpolation='gaussian')
    x0 = m(np.min(lons), np.min(lats))
    x1 = m(np.max(lons), np.max(lats))
    im.set_extent([x0[0], x1[0], x0[1], x1[1]])
    cf = ax.contourf(lons, lats, tmax, levels=np.arange(vmin, vmax+620, 620), cmap=cm.lapaz, alpha=0.75, antialiased=True)
    ax.contour(lons, lats, tmax, levels=np.arange(vmin, vmax, 620), zorder=10, linewidths=0.3, colors='k')
    ax.contour(lons, lats, tmax, levels=np.arange(0, 0+620, 620), zorder=10, linewidths=0.5, colors='k')
    
    axp = ax.get_position()
    cax=fig.add_axes([axp.x1+0.005, axp.y1-axp.height/2, 0.01, axp.height/2])
    plt.colorbar(cf, cax=cax, label='Elevation (m)')

    
# connvert coordinates of moment tensor from USE to NED
def convertUSEtoNED(focmec):
    '''
    convert basis-moment tensors
        from mrr, mtt, mpp, mrt, mrp, mtp
          to mnn, mee, mdd, mne, mnd, med
    https://pyrocko.org/docs/current/library/reference/moment_tensor.html#pyrocko.moment_tensor.MomentTensor
    usage:
    focmec = convertUSEtoNED(focmec[0],focmec[1],focmec[2],focmec[3],focmec[4],focmec[5])
    '''
    mrr = focmec[0]
    mtt = focmec[1]
    mpp = focmec[2]
    mrt = focmec[3]
    mrp = focmec[4]
    mtp = focmec[5]
    #
    return [mtt, mpp, mrr, -mtp, mrt, -mrp]
    
    
# Load focal mechanism compilation
def loadmecs(datarootdir):
    meclocs = [] # true location
    mos = [] # seismic moment Nm
    focmecs = [] # moment tensors in USE system (Pyrocko)
    labels = [] # label appended to beachball

    # This study FFM
    model=np.loadtxt(datarootdir+'modellist.txt', dtype=int, usecols=0)
    for j in np.arange(13, 14, 1):
        model_para = load_fort40(datarootdir+'model_'+str(model[j])+'/fort.40')

        data=np.loadtxt(datarootdir+'model_'+str(model[j])+'/FFM_DCpreferred.txt', skiprows=1)
        lon, lat, slip, strike, dip, rake = data[:,2], data[:,3],data[:,1], data[:,7], data[:,8],data[:,9]
        meclocs.append([lon[np.argmax(slip)], lat[np.argmax(slip)]])

        mos.append(model_para.moment[0])

        focmec = convertUSEtoNED([0.265296, -1.010935,  0.745638, -0.306760, -0.106615, -0.030290])
        focmecs.append(focmec)
    labels.append('FFM')

    # Steve's R-CMT Single, Low-frequency
    meclocs.append([179.8000000000, -37.3700000000])
    with open('../materials/R-CMT/LowFreq_SingleSubEvent_LineDepths/inv1.dat') as f:
        data = f.readlines()[51]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]
    mos.append(moment)
    #a1, a2, a3, a4, a5, a6 = .108957E+20,.314176E+20,  .410898E+20,  .820662E+20, -.254468E+20,  .000000E+00
    #a1, a2, a3, a4, a5, a6=.574851E+19,  .302531E+20,  .389741E+20,  .880749E+20, -.253843E+20,  .000000E+00
    with open('../materials/R-CMT/LowFreq_SingleSubEvent_LineDepths/inv1.dat') as f:
        data = f.readlines()[47]
    a1, a2, a3, a4, a5, a6 = np.array(data.strip().split()).astype(float)
    Mxx=-a4+a6
    Myy=-a5+a6
    Mzz=a4+a5+a6
    Mxy=a1
    Mxz=a2
    Myz=-a3
    focmecs.append([Mxx, Myy, Mzz, Mxy, Mxz, Myz])
    labels.append('R-CMT')

    # Steve's W-phase
    mean = load_events(datarootdir+"WphaseCMT/mean_wphase.yaml")[0]  # Read in event file containing best solution
    azi = 90-np.rad2deg(np.arctan2(mean.north_shift, mean.east_shift))
    distance = np.sqrt(mean.north_shift**2+mean.east_shift**2)
    tmp = geod.Direct(mean.lat, mean.lon, azi, distance)
    meclocs.append([tmp['lon2'], tmp['lat2']])
    mos.append(mean.moment_tensor.moment)
    focmecs.append(mean.moment_tensor)
    labels.append('W-phase')

    # GCMT https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr=2021&mo=3&day=4&oyr=1976&omo=1&oday=1&jyr=1976&jday=1&ojyr=1976&ojday=1&otype=nd&nday=1&lmw=0&umw=10&lms=0&ums=10&lmb=0&umb=10&llat=-90&ulat=90&llon=-180&ulon=180&lhd=0&uhd=1000&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=0
    meclocs.append([179.87, -37.36])
    import obspy
    cat_gcmt = obspy.read_events('../materials/work/SPUD_QUAKEML_bundle.xml')
    for event in cat_gcmt.events:
        if event.resource_id == 'smi:service.iris.edu/fdsnws/event/1/query?eventid=11384594':
            moment = event.focal_mechanisms[0].moment_tensor.scalar_moment
            #mw = event.magnitudes[0].mag
            focmec=[event.focal_mechanisms[0].moment_tensor.tensor.m_rr,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_tt,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_pp,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_rt,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_rp,
                     event.focal_mechanisms[0].moment_tensor.tensor.m_tp]


    mos.append(moment)
    #focmec = [0.364, -1.030, 0.664, 0.036, 0.071, -0.213 ]
    focmec = convertUSEtoNED(focmec)
    focmecs.append(focmec)
    labels.append('GCMT')

    # Steve's R-CMT, Multi-point, high-frequency, along trench 1    52
    meclocs.append([179.9461002352,   -37.0848574759])
    with open('../materials/R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[106]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]
    mos.append(moment)
    with open('../materials/R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[102]
    a1, a2, a3, a4, a5, a6 = np.array(data.strip().split()).astype(float)
    #a1, a2, a3, a4, a5, a6=.238361E+20,  .545402E+20,  .558320E+20,  .602305E+20, -.363396E+20, .000000E+00
    Mxx=-a4+a6
    Myy=-a5+a6
    Mzz=a4+a5+a6
    Mxy=a1
    Mxz=a2
    Myz=-a3
    focmecs.append([Mxx, Myy, Mzz, Mxy, Mxz, Myz])
    labels.append('R-CMT HF1')

    # Steve's R-CMT, Multi-point, high-frequency, along trench 2    14
    meclocs.append([179.8308942630,   -37.3389833860])
    with open('../materials/R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[228]
    moment = np.array(data.replace('moment (Nm):', '').strip().split()).astype(float)[0]
    mos.append(moment)
    with open('../materials/R-CMT/HighFreq_MultipleSubEvents_TrenchParallelPlane/inv1.dat') as f:
        data = f.readlines()[224]
    a1, a2, a3, a4, a5, a6 = np.array(data.strip().split()).astype(float)
    #a1, a2, a3, a4, a5, a6=-.781991E+19,  .897654E+19, -.446557E+19, -.103638E+20, -.927517E+19,  .000000E+00
    Mxx=-a4+a6
    Myy=-a5+a6
    Mzz=a4+a5+a6
    Mxy=a1
    Mxz=a2
    Myz=-a3
    focmecs.append([Mxx, Myy, Mzz, Mxy, Mxz, Myz])
    labels.append('R-CMT HF2')


    return meclocs, mos, focmecs, labels


# load model parameters from `fort.40` and store the values as pandas data frame
def load_fort40(infile):
    col_names = [ 'c{0:02d}'.format(i) for i in range(10)]
    df = pd.read_table(infile, names=col_names, header=None, delimiter='\t')
    tmp1 = [ float(x) for x in df['c00'][1].replace(' ', ',').split(',') if x ]
    tmp3 = [ float(x) for x in df['c00'][3].replace(' ', ',').split(',') if x ]
    tmp5 = [ float(x) for x in df['c00'][5].replace(' ', ',').split(',') if x ]
    tmp7 = [ float(x) for x in df['c00'][7].replace(' ', ',').split(',') if x ]
    df = pd.DataFrame({ 'moment'   : tmp1[0],
                        'mw'       : tmp1[1],
                        'rigidity' : tmp1[2],
                        'lat'      : tmp1[3],
                        'lon'      : tmp1[4],
                        'depth'    : tmp1[5],
                        'vr'       : tmp1[6],
                        'nsurface' : tmp1[7],

                        'strike'   : tmp3[0],
                        'dip'      : tmp3[1],

                        'xx'       : tmp5[0],
                        'yy'       : tmp5[1],
                        'mn'       : int(tmp5[2]),
                        'nn'       : int(tmp5[3]),
                        'm0'       : int(tmp5[4]),
                        'n0'       : int(tmp5[5]),
                        'tr'       : tmp5[6],
                        'jtn'      : int(tmp5[7]),
                        'icmn'     : int(tmp5[8]),

                        'variance' : tmp7[0],

                      }, index=[0])
    return df


def getSlipCell3D(datarootdir, model_para, modelid, ax, zorder=1):
    ###
    cmap = cm.bilbao
    elon, elat = model_para.lon[0], model_para.lat[0]
    nx, ny, dx, dy, x0, y0 = model_para.mn[0], model_para.nn[0], model_para.xx[0], model_para.yy[0], model_para.m0[0], model_para.n0[0]
    model_dip = model_para.dip[0]
    model_stk = model_para.strike[0]

    data=np.loadtxt(datarootdir+'model_'+str(modelid)+'/FFM_DCpreferred.txt', skiprows=1)
    lon, lat, slip, strike, dip, rake = data[:,2], data[:,3],data[:,1], data[:,7], data[:,8],data[:,9]
    x, y, depth = data[:,10], data[:,11], data[:,4]


    depmin=model_para.depth[0] - np.sin( np.deg2rad(model_para.dip[0]) )*((ny-y0)*dy+dy/2)
    depmax=model_para.depth[0] + np.sin( np.deg2rad(model_para.dip[0]) )*((y0-1)*dy+dy/2)
    #print(ny, y0, dy, depmin, model_para.depth[0])

    if elon < 0: elon += 360

    vertslist = []
    # top: center to right
    for i in range(len(depth)):
    #for i in [0]:
        shiftk = (dy/2) * np.cos(np.deg2rad(model_dip))
        tmp0 = geod.Direct(lat[i], lon[i], model_stk-90, shiftk*1e3)
        if tmp0['lon2'] < 0: tmp0['lon2'] = tmp0['lon2'] + 360
        if tmp0['lon1'] < 0: tmp0['lon1'] += 360
        tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk, (dx/2)*1e3)
        if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
        if tmp1['lon1'] < 0: tmp1['lon1'] += 360
        x1, y1 = (tmp1['lon1'], tmp1['lat1'])
        x2, y2 = (tmp1['lon2'], tmp1['lat2'])

        RTx = x2
        RTy = y2
        RTz = depth[i] - (dy/2) * np.sin(np.deg2rad(model_dip))


        # top: center to left
        tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk, (-dx/2)*1e3)
        if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
        if tmp1['lon1'] < 0: tmp1['lon1'] += 360
        x1, y1 = (tmp1['lon1'], tmp1['lat1'])
        x2, y2 = (tmp1['lon2'], tmp1['lat2'])

        LTx = x2
        LTy = y2
        LTz = depth[i] - (dy/2) * np.sin(np.deg2rad(model_dip))


        # bottom: center to right
        shiftk = (dy/2) * np.cos(np.deg2rad(model_dip))
        tmp0 = geod.Direct(lat[i], lon[i], model_stk+90, shiftk*1e3)
        if tmp0['lon2'] < 0: tmp0['lon2'] = tmp0['lon2'] + 360
        if tmp0['lon1'] < 0: tmp0['lon1'] += 360
        tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk, (dx/2)*1e3)
        if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
        if tmp1['lon1'] < 0: tmp1['lon1'] += 360
        x1, y1 = (tmp1['lon1'], tmp1['lat1'])
        x2, y2 = (tmp1['lon2'], tmp1['lat2'])

        RBx = x2
        RBy = y2
        RBz = depth[i] + (dy/2) * np.sin(np.deg2rad(model_dip))


        # bottom: center to left
        tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk, (-dx/2)*1e3)
        if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
        if tmp1['lon1'] < 0: tmp1['lon1'] += 360
        x1, y1 = (tmp1['lon1'], tmp1['lat1'])
        x2, y2 = (tmp1['lon2'], tmp1['lat2'])

        LBx = x2
        LBy = y2
        LBz = depth[i] + (dy/2) * np.sin(np.deg2rad(model_dip))

        polyx = [LBx, RBx, RTx, LTx, LBx]
        polyy = [LBy, RBy, RTy, LTy, LBy]
        polyz = [LBz, RBz, RTz, LTz, LBz]
        verts = [list(zip(polyx, polyy, polyz))]

        vertslist.append(verts)

    return vertslist


def plotModelEdge3D(model_para, ax, zorder=1):
    elon, elat = model_para.lon[0], model_para.lat[0]
    nx, ny, dx, dy, x0, y0 = model_para.mn[0], model_para.nn[0], model_para.xx[0], model_para.yy[0], model_para.m0[0], model_para.n0[0]
    model_dip = model_para.dip[0]
    model_stk = model_para.strike[0]

    depmin=model_para.depth[0] - np.sin( np.deg2rad(model_para.dip[0]) )*((ny-y0)*dy+dy/2)
    depmax=model_para.depth[0] + np.sin( np.deg2rad(model_para.dip[0]) )*((y0-1)*dy+dy/2)

    if elon < 0: elon += 360

    # top: center to right
    shiftk = ((ny-y0)*dy+dy/2) * np.cos(np.deg2rad(model_dip))
    tmp0 = geod.Direct(elat, elon, model_stk-90, shiftk*1e3)
    if tmp0['lon2'] < 0: tmp0['lon2'] = tmp0['lon2'] + 360
    if tmp0['lon1'] < 0: tmp0['lon1'] += 360
    tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk, ((nx-x0)*dx+dx/2)*1e3)
    if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
    if tmp1['lon1'] < 0: tmp1['lon1'] += 360
    x1, y1 = (tmp1['lon1'], tmp1['lat1'])
    x2, y2 = (tmp1['lon2'], tmp1['lat2'])
    if model_dip > 0:
        ax.plot([x1, x2], [y1, y2], [depmin, depmin], lw=1, color='k', solid_capstyle='projecting', zorder=zorder+1)
    else:
        ax.plot([x1, x2], [y1, y2], [depmin, depmin], lw=1, color='C7', solid_capstyle='projecting', zorder=zorder+1)


    RTx = x2
    RTy = y2
    RTz = depmin


    # top: center to left
    tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk, ((-x0)*dx+dx/2)*1e3)
    if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
    if tmp1['lon1'] < 0: tmp1['lon1'] += 360
    x1, y1 = (tmp1['lon1'], tmp1['lat1'])
    x2, y2 = (tmp1['lon2'], tmp1['lat2'])
    if model_dip > 0:
        ax.plot([x1, x2], [y1, y2], [depmin, depmin], lw=1, color='k', solid_capstyle='projecting', zorder=zorder+1)
    else:
        ax.plot([x1, x2], [y1, y2], [depmin, depmin], lw=1, color='C7', solid_capstyle='projecting', zorder=zorder+1)

    LTx = x2
    LTy = y2
    LTz = depmin


    # bottom: center to right
    shiftk = ((y0-1)*dy+dy/2) * np.cos(np.deg2rad(model_dip))
    tmp0 = geod.Direct(elat, elon, model_stk+90, shiftk*1e3)
    if tmp0['lon2'] < 0: tmp0['lon2'] = tmp0['lon2'] + 360
    if tmp0['lon1'] < 0: tmp0['lon1'] += 360
    tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk, ((nx-x0)*dx+dx/2)*1e3)
    if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
    if tmp1['lon1'] < 0: tmp1['lon1'] += 360
    x1, y1 = (tmp1['lon1'], tmp1['lat1'])
    x2, y2 = (tmp1['lon2'], tmp1['lat2'])
    ax.plot([x1, x2], [y1, y2], [depmax, depmax], lw=1, color='C7', solid_capstyle='butt', zorder=zorder)

    RBx = x2
    RBy = y2
    RBz = depmax


    # bottom: center to left
    tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk, ((-x0)*dx+dx/2)*1e3)
    if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
    if tmp1['lon1'] < 0: tmp1['lon1'] += 360
    x1, y1 = (tmp1['lon1'], tmp1['lat1'])
    x2, y2 = (tmp1['lon2'], tmp1['lat2'])
    ax.plot([x1, x2], [y1, y2], [depmax, depmax], lw=1, color='C7', solid_capstyle='butt', zorder=zorder)

    LBx = x2
    LBy = y2
    LBz = depmax


    # right: hypo depth to bottom
    shiftl = (nx-x0)*dx+dx/2
    shiftk = ((y0-1)*dy+dy/2) * np.cos(np.deg2rad(model_dip))
    tmp0 = geod.Direct(elat, elon, model_stk, shiftl*1e3)
    if tmp0['lon2'] < 0: tmp0['lon2'] = tmp0['lon2'] + 360
    if tmp0['lon1'] < 0: tmp0['lon1'] += 360
    tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk+90, shiftk*1e3)
    if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
    if tmp1['lon1'] < 0: tmp1['lon1'] += 360
    x1, y1 = (tmp1['lon1'], tmp1['lat1'])
    x2, y2 = (tmp1['lon2'], tmp1['lat2'])
    ax.plot([x1, x2], [y1, y2], [model_para.depth[0], depmax], lw=1, color='C7', solid_capstyle='projecting', zorder=zorder)


    # right: hypo depth to top
    shiftk = ((ny-y0)*dy+dy/2) * np.cos(np.deg2rad(model_dip))
    tmp0 = geod.Direct(elat, elon, model_stk, shiftl*1e3)
    if tmp0['lon2'] < 0: tmp0['lon2'] = tmp0['lon2'] + 360
    if tmp0['lon1'] < 0: tmp0['lon1'] += 360
    tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk-90, shiftk*1e3)
    if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
    if tmp1['lon1'] < 0: tmp1['lon1'] += 360
    x1, y1 = (tmp1['lon1'], tmp1['lat1'])
    x2, y2 = (tmp1['lon2'], tmp1['lat2'])
    ax.plot([x1, x2], [y1, y2], [model_para.depth[0], depmin], lw=1, color='C7', solid_capstyle='projecting', zorder=zorder)

    # left: hypo depth to bottom
    shiftl = (x0-1)*dx+dx/2
    shiftk = ((y0-1)*dy+dy/2) * np.cos(np.deg2rad(model_dip))
    tmp0 = geod.Direct(elat, elon, model_stk-180, shiftl*1e3)
    if tmp0['lon2'] < 0: tmp0['lon2'] = tmp0['lon2'] + 360
    if tmp0['lon1'] < 0: tmp0['lon1'] += 360
    tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk+90, shiftk*1e3)
    if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
    if tmp1['lon1'] < 0: tmp1['lon1'] += 360
    x1, y1 = (tmp1['lon1'], tmp1['lat1'])
    x2, y2 = (tmp1['lon2'], tmp1['lat2'])
    ax.plot([x1, x2], [y1, y2], [model_para.depth[0], depmax], lw=1, color='C7', solid_capstyle='projecting', zorder=zorder)


    # right: hypo depth to top
    shiftk = ((ny-y0)*dy+dy/2) * np.cos(np.deg2rad(model_dip))
    tmp0 = geod.Direct(elat, elon, model_stk-180, shiftl*1e3)
    if tmp0['lon2'] < 0: tmp0['lon2'] = tmp0['lon2'] + 360
    if tmp0['lon1'] < 0: tmp0['lon1'] += 360
    tmp1 = geod.Direct(tmp0['lat2'], tmp0['lon2'], model_stk-90, shiftk*1e3)
    if tmp1['lon2'] < 0: tmp1['lon2'] = tmp1['lon2'] + 360
    if tmp1['lon1'] < 0: tmp1['lon1'] += 360
    x1, y1 = (tmp1['lon1'], tmp1['lat1'])
    x2, y2 = (tmp1['lon2'], tmp1['lat2'])
    ax.plot([x1, x2], [y1, y2], [model_para.depth[0], depmin], lw=1, color='C7', solid_capstyle='projecting', zorder=zorder)

    polyx = [LBx, RBx, RTx, LTx, LBx]
    polyy = [LBy, RBy, RTy, LTy, LBy]
    polyz = [LBz, RBz, RTz, LTz, LBz]
    verts = [list(zip(polyx, polyy, polyz))]


    return verts

    #ax.scatter(LBx, LBy, LBz)
    #ax.scatter(RBx, RBy, RBz)
    #ax.scatter(RTx, RTy, RTz)
    #ax.scatter(LTx, LTy, LTz)

    
# select preferred fault plane based on the model plane geometry
def selectplane(modelstk, modeldip, stk0, dip0, rake0, stk1, dip1, rake1):
    vecmodelplane = faultnormalvec(modelstk, modeldip)
    vecplane0 = faultnormalvec(stk0, dip0)
    vecplane1 = faultnormalvec(stk1, dip1)
    tmp0 = np.inner(vecmodelplane, vecplane0)
    tmp1 = np.inner(vecmodelplane, vecplane1)
    if abs(tmp0) > abs(tmp1):
        stk_s = stk0
        dip_s = dip0
        rake_s = rake0
    elif abs(tmp0) < abs(tmp1):
        stk_s = stk1
        dip_s = dip1
        rake_s = rake1
    else:
        stk_s = stk0
        dip_s = dip0
        rake_s = rake0
    return stk_s, dip_s, rake_s

def faultnormalvec(stk, dip):
    nn = -np.sin(np.deg2rad(stk)) * np.sin(np.deg2rad(dip))
    ne =  np.cos(np.deg2rad(stk)) * np.sin(np.deg2rad(dip))
    nd = -np.cos(np.deg2rad(dip))
    return np.array([ne, nn, nd])


def get_plunge_azimuth(focmec):
    '''
    input: 6 components of  mooment tensor in USE convention (e.g., default GCMT, psmeca format)
    output: A list of [p_plunge, p_azimuth, n_plunge, n_azimuth, t_plunge, t_azimuth] angles in degree
    '''
    m1,m2,m3,m4,m5,m6 = focmec[0],focmec[1],focmec[2],focmec[3],focmec[4],focmec[5]
    mt = obspy.imaging.beachball.MomentTensor(m1,m2,m3,m4,m5,m6, 26)

    (d, v) = np.linalg.eigh(mt.mt)
    pl = np.arcsin(-v[0])
    az = np.arctan2(v[2], -v[1])
    for i in range(0, 3):
        if pl[i] <= 0:
            pl[i] = -pl[i]
            az[i] += np.pi
        if az[i] < 0:
            az[i] += 2 * np.pi
        if az[i] > 2 * np.pi:
            az[i] -= 2 * np.pi

    pl = np.rad2deg(pl)
    az = np.rad2deg(az)

    p_plunge = pl[0]
    p_azimuth = az[0]
    n_plunge = pl[1]
    n_azimuth = az[1]
    t_plunge = pl[2]
    t_azimuth = az[2]

    return [p_plunge, p_azimuth, n_plunge, n_azimuth, t_plunge, t_azimuth]


def get_unique_list(seq):
    seen = []
    return [x for x in seq if x not in seen and not seen.append(x)]


def aziequi_wphase(ax, stalist, loc, distancetextsize):

    azidistlist = []
    #print(len(stalist))
    for i in range(len(stalist)):
        tmp = geod.Inverse(loc[1], loc[0], stalist[i][1], stalist[i][2])
        d, a = tmp['a12'], 90-tmp['azi1']

        x, y=(d*np.cos(a*np.pi/180.0), d*np.sin(a*np.pi/180.0))
        sc=ax.scatter(x, y, s=15, marker='^', edgecolor='k', facecolor='ivory', alpha=1, lw=0.75, zorder=10)

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
