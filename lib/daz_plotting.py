#!/usr/bin/env python3

#from daz_lib import *

#for visualisation and kml export
#import hvplot
#import holoviews as hv
#from holoviews import opts
#hv.extension('bokeh')
import simplekml
import datetime as dt
import pandas as pd
import numpy as np
import xarray as xr
import glob, os

# 2024 changing to pygmt from hv
import pygmt

def plot_esds_from_pd(esds):
    """ will quickly plot esds taken by e.g. dazes=daz_lib_licsar.get_daz_frame(frame)
    """
    if 'epochdate' not in esds:
        esds['epochdate'] = esds.apply(lambda x : pd.to_datetime(str(x.epoch)).date(), axis=1)
    (esds.set_index(esds.epochdate).daz*14000).plot()

'''
def plot_vel_esd(frame_esds, frameta, level1 = 'iono_grad', level2 = None, showitrf=True,
                 mindate = dt.date(2014, 11, 1), maxdate = dt.date(2024, 11, 1)):
    # level1: 'tide', 'orig', 'iono_grad', 'iono_f2'
    #global mindate
    #global maxdate
    # ionotype - grad for gradient only (Gomba, 2016 - partial), f2 to use also the approach with F2 heights (Liang 2019)
    xmin = 0
    xmax = (maxdate-mindate).days/365.25
    acq_times = frame_esds['epochdate'].values
    frame_esds = frame_esds.set_index('years_since_beginning')
    frame = frameta['frame'].values[0]
    #level 1 points
    # which col to use for points?
    if level1 == 'tide':
        col_mm = 'daz_mm_notide'
        col_outliers = 'is_outlier_daz_mm_notide'
        col_label = 'daz notide'
        #colour = 'orange'
        colour = 'brown'
    elif level1 == 'orig':
        col_mm = 'daz_mm'
        col_outliers = 'is_outlier'
        col_label = 'daz orig'
        colour = 'red'
    elif level1 == 'iono_grad':
        col_mm = 'daz_mm_notide_noiono_grad'
        if not col_mm in frame_esds:
            col_mm = 'daz_mm_notide_noiono'
        col_outliers = 'is_outlier_'+col_mm
        col_label = 'daz notide, noiono'
        colour = 'green'
    elif level1 == 'iono_grad_OK':
        col_mm = 'daz_mm_notide_noiono_grad_OK'
        col_outliers = 'is_outlier_'+col_mm
        col_label = 'daz notide, noiono'
        colour = 'green'
    elif level1 == 'iono_f2':
        col_mm = 'daz_mm_notide_noiono_F2'
        col_outliers = 'is_outlier_'+col_mm
        col_label = 'daz notide, noiono (+F2)'
        colour = 'green'
    else:
        print("wrong type - try one of those: 'tide', 'orig', 'iono_grad', 'iono_f2'")
        return False
    if level2:
        pointsize = 2
        line_width = 1
    else:
        pointsize = 5
        line_width=5
    #scatter plot - all points (good and bad)
    hv_scatter_good = hv.Scatter(frame_esds[frame_esds[col_outliers] == False][col_mm], label=col_label).opts(size=pointsize, color = colour)
    #hv_scatter_bad = hv.Scatter(frame_esds[frame_esds[col_outliers] == True][col_mm], label=col_label+' (outliers)').opts(size=pointsize-1, marker='x', color = colour) #color = '#5ebaff')
    hv_scatter_bad = hv.Scatter(frame_esds[frame_esds[col_outliers] == True][col_mm], label='').opts(size=pointsize-1, marker='x', color = colour) #color = '#5ebaff')
    #slopes..
    slope = float(frameta['slope_'+col_mm+'_mmyear'])
    intercept = float(frameta['intercept_'+col_mm+'_mmyear'])
    std = float(frameta[col_mm+'_RMSE_selection'])
    #bycol+'_RMSE'
    rmse = float(frameta[col_mm+'_RMSE_full'])
    rmse2 = float(frameta[col_mm+'_RMSE_mmy_full'])
    #
    hv_slope = hv.Slope(slope, intercept).opts(color=colour, line_width=line_width)
    #hv_slope = hv_slope * hv.Curve([-99999, -99998], label = 'vel: {0} mm/y, RMSE: {1} mm'.format(round(slope),round(rmse))).opts(color=colour, line_width=line_width)
    hv_slope = hv_slope * hv.Curve([-99999, -99998], label = 'vel: {0} +-{1} mm/y (95%)'.format(round(slope),round(rmse2))).opts(color=colour, line_width=line_width)
    hv_slopestd1 = hv.Slope(slope, intercept + 2*std).opts(color=colour, line_width=1, line_dash='dashed')
    #hv_slopestd1 = hv_slopestd1.opts(linestyle='dashed', color='red', line_width=1)
    hv_slopestd2 = hv.Slope(slope, intercept - 2*std).opts(color=colour, line_width=1, line_dash='dashed')
    #hv_slopestd2 = hv_slopestd2.opts(linestyle='dashed', color='red', line_width=1)
    hvplot = hv_scatter_bad * hv_scatter_good * hv_slope
    if not level2:
        hvplot = hvplot * hv_slopestd1 * hv_slopestd2
    # level2
    if level2:
        if level2 == 'tide':
            col_mm = 'daz_mm_notide'
            col_outliers = 'is_outlier_daz_mm_notide'
            col_label = 'daz notide'
            colour = 'brown'
        elif level2 == 'orig':
            col_mm = 'daz_mm'
            col_outliers = 'is_outlier'
            col_label = 'daz orig'
            colour = 'red'
        elif level2 == 'iono_grad':
            col_mm = 'daz_mm_notide_noiono_grad'
            if not col_mm in frame_esds:
                col_mm = 'daz_mm_notide_noiono'
            col_outliers = 'is_outlier_'+col_mm
            col_label = 'daz notide, noiono'
            colour = 'green'
        elif level2 == 'iono_grad_OK':
            col_mm = 'daz_mm_notide_noiono_grad_OK'
            col_outliers = 'is_outlier_'+col_mm
            col_label = 'daz notide, noiono'
            colour = 'green'
        elif level2 == 'iono_f2':
            col_mm = 'daz_mm_notide_noiono_F2'
            col_outliers = 'is_outlier_'+col_mm
            col_label = 'daz notide, noiono (+F2)'
            colour = 'green'
        else:
            print("wrong type - try one of those: 'tide', 'orig', 'iono_grad', 'iono_f2'")
            return False
        pointsize = 5
        #if level2 == 'tide':
        #    ........................................
        #scatter plot - all points (good and bad)
        hv_scatter_good = hv.Scatter(frame_esds[frame_esds[col_outliers] == False][col_mm], label=col_label).opts(size=pointsize, color = colour)
        #hv_scatter_bad = hv.Scatter(frame_esds[frame_esds[col_outliers] == True][col_mm], label=col_label+' (outliers)').opts(size=pointsize-1, marker='x', color = colour) #color = '#5ebaff')
        hv_scatter_bad = hv.Scatter(frame_esds[frame_esds[col_outliers] == True][col_mm], label='').opts(size=pointsize-1, marker='x', color = colour) #color = '#5ebaff')
        #slopes..
        slope = float(frameta['slope_'+col_mm+'_mmyear'])
        intercept = float(frameta['intercept_'+col_mm+'_mmyear'])
        std = float(frameta[col_mm+'_RMSE_selection'])
        rmse = float(frameta[col_mm+'_RMSE_full'])
        rmse2 = float(frameta[col_mm+'_RMSE_mmy_full'])
        #
        hv_slope = hv.Slope(slope, intercept).opts(color=colour, line_width=3)
        #hv_slope = hv_slope * hv.Curve([-99999, -99998], label = 'vel: {0} mm/y, RMSE: {1} mm'.format(round(slope),round(rmse))).opts(color=colour, line_width=2)
        hv_slope = hv_slope * hv.Curve([-99999, -99998], label = 'vel: {0} +-{1} mm/y (95%)'.format(round(slope),round(rmse2))).opts(color=colour, line_width=2)
        hv_slopestd1 = hv.Slope(slope, intercept + 2*std).opts(color=colour, line_width=1, line_dash='dashed')
        hv_slopestd2 = hv.Slope(slope, intercept - 2*std).opts(color=colour, line_width=1, line_dash='dashed')
        hvplot = hvplot * hv_scatter_bad * hv_scatter_good * hv_slope
        hvplot = hvplot * hv_slopestd1 * hv_slopestd2
    #
    if showitrf:
        slope_itrf = float(frameta['slope_plates_vel_azi_itrf2014'])
        intercept_itrf = np.median(frame_esds[frame_esds[col_outliers] == False][col_mm])
        hv_slope_itrf = hv.Slope(slope_itrf, intercept_itrf).opts(color='black', line_width=2, line_alpha=0.5) #, line_dash='dashed')
        hv_slope_itrf = hv_slope_itrf * hv.Curve([-99999, -99998], label = 'ITRF: {} mm/y'.format(round(slope_itrf))).opts(color='black', line_width=2, line_alpha=0.5) #line_dash='dashed')
        #perhaps add value!
        #
        #
        hvplot = hvplot * hv_slope_itrf
    #else:
    #additional_label = ', RMSE (no outliers)= '+str(int(std))+' mm'
    #hvplot.label = 'frame '+frame+': final velocity = '+str(round(slope))+' mm/y'+additional_label
    hvplot.label = 'frame '+frame #+': final velocity = '+str(round(slope))+' mm/y'+additional_label
    hvplot = hvplot.opts(legend_position='bottom_left', legend_cols=1)
    xtickdates = pd.date_range(start=mindate, end=maxdate, freq=pd.offsets.MonthBegin(3)).astype(str)
    xtickvals = []
    for xdate in xtickdates:
        xtickvals.append((pd.Timestamp(xdate)- pd.Timestamp(mindate)).days/365.25)
    #xticksnames = list(zip(frame_esds.index, acq_times.astype(str)))
    xticksnames = list(zip(xtickvals, xtickdates))
    #tickdif = int(frameta['count_all']/30)
    #hvplot = hvplot.opts(xticks=xticksnames[::tickdif],xrotation=90)
    hvplot = hvplot.opts(xticks=xticksnames,xrotation=90)
    #xmin = frame_esds.index[0]
    #xmax = frame_esds.index[-1]
    #hvplot = hvplot.opts(fig_size=300, aspect=2, xlim=(xmin, xmax))
    hvplot = hvplot.opts(xlim=(xmin, xmax))
    hvplot = hvplot.opts(xlabel='', ylabel='azimuth shift [mm]', ylim=(-300, 300), show_legend=True)
    hvplot = hvplot.options(width=900, height=350)
    return hvplot
'''




def export_esds2kml(framespd, esds, kmzfile = 'esds.kmz', level1 = 'tide', level2 = 'iono_grad', overwrite = False, clean = True):
    if os.path.exists('plots') or os.path.exists('doc.kml') or os.path.exists(kmzfile):
        if not overwrite:
            print('The kmz files already exist, please recheck')
            return False
        else:
            os.system('rm -r {} doc.kml plots'.format(kmzfile))
    os.mkdir('plots')
    kml = simplekml.Kml()
    print('generating plots')
    # getting min/max date
    try:
        allepochs = esds[esds.is_outlier_daz_mm == False].epochdate.values
        allepochs.sort()
        mindate = allepochs[0]
        maxdate = allepochs[-1]
    except:
        print('Error getting min/max dates. Trying from all epochdates')
        allepochs = esds.epochdate.values
        allepochs.sort()
        mindate = allepochs[0]
        maxdate = allepochs[-1]
    # extract years since the first date (careful, might be minus then...)
    esds['years_since_beginning'] = esds['epochdate'] - mindate
    esds['years_since_beginning'] = esds['years_since_beginning'].apply(lambda x: float(x.days) / 365.25)
    #this will generate plots
    lenframes = len(framespd['frame'])
    i = 0
    for frame in framespd['frame']:
        i=i+1
        print('  Running for {0:6}/{1:6}th frame...'.format(i, lenframes), flush=True, end='\r')
        frameta = framespd[framespd['frame']==frame]
        selected_frame_esds = esds[esds['frame'] == frame].copy()
        #frameplot = plot_vel_esd(selected_frame_esds, frameta, showtec = False)
        try:
            #frameplot = plot_vel_esd(selected_frame_esds, frameta, level2 = level2, level1 = level1, showitrf=True, mindate = mindate, maxdate = maxdate)
            #hv.save(frameplot, 'plots/{}.png'.format(frame), dpi=100, fmt='png')
            frameplot = plot_vel_esd_gmt(selected_frame_esds, frameta, level2=level2, level1=level1, showitrf=True,
                                     mindate=mindate, maxdate=maxdate)
            frameplot.savefig('plots/{}.png'.format(frame), dpi=150)
        except:
            print('error generating plot for frame '+frame)
    #this will generate kmz
    print('generating kml')
    for dirpass in [('A','ascending'),('D','descending')]:
        sname = dirpass[1]
        fol = kml.newfolder(name=sname)
        for i, frameta in framespd[framespd['frame'].str[3] == dirpass[0]].iterrows():
            frame = frameta['frame']
            plotpath = os.path.join('plots', frame+'.png')
            if os.path.exists(plotpath):
                lon = frameta.center_lon
                lat = frameta.center_lat
                point = fol.newpoint(name = frame , coords = [(lon,lat)])
                point.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
                point.style.balloonstyle.text = "<![CDATA[ <table width=900px cellpadding=0 cellspacing=0> <tr><td><img width=900px src='" + plotpath + "' /></td></tr></table>]]>"
    kml.save("doc.kml")
    os.system('7za a temp.zip doc.kml plots >/dev/null 2>/dev/null; mv temp.zip {}'.format(kmzfile))
    if clean:
        if os.path.exists(kmzfile):
            os.system('rm -r doc.kml plots')


def plot_vel_esd_gmt(selected_frame_esds, frameta, mindate, maxdate, level1, level2=None, showitrf=True):
    frame = frameta['frame'].values[0]
    #
    fig = pygmt.Figure()
    #
    fig.basemap(
        projection="X18c/6c",
        region=[mindate, maxdate, -300, 300], #datetime.date(2010, 1, 1), datetime.date(2020, 6, 1), 0, 10],
        #frame=["WSen", "af"],
        #frame=["WSne", "xaf", "yaf+l'daz [mm]'"]
        frame=["a", "+t "+frame, "xafg", "yafg+ldaz [mm]"]
    )
    #
    if level2:
        fig = figpart_var(level1, selected_frame_esds, frameta, fig, additrf = (True & showitrf), plotstd = False)
        fig = figpart_var(level2, selected_frame_esds, frameta, fig, additrf = (False & showitrf), plotstd = True)
    else:
        fig = figpart_var(level1, selected_frame_esds, frameta, fig, additrf = (True & showitrf), plotstd = True)
    #
    #print('legend')
    fig.legend(position="JBL+jBL+o0.1c", box='+gwhite+p1p')
    fig.basemap(frame=True) #["WSen", "af"])
    #fig.show(dpi=120)
    return fig


def figpart_var(level, esdspart, frameta, fig, additrf=False, plotstd=False):
    frame_esds = esdspart.copy()
    #
    if level == 'tide':
        # which column
        col_mm = 'daz_mm_notide'
        col_outliers = 'is_outlier_' + col_mm
        col_label = 'tide corrected'
        # col_color = 'red'
        col_color = 'olivedrab'
        col_size = '0.1c'
    elif level == 'iono_grad':
        # or:
        col_mm = 'daz_mm_notide_noiono_grad'
        if not col_mm in frame_esds:
            col_mm = 'daz_mm_notide_noiono'
        col_outliers = 'is_outlier_' + col_mm
        col_label = 'tide and iono corrected'
        # col_color = 'olivedrab'
        col_color = 'red'
        col_size = '0.15c'
        col_size = '0.2c'
    #
    # slopes..
    slope = float(frameta['slope_' + col_mm + '_mmyear'].values[0])
    intercept = float(frameta['intercept_' + col_mm + '_mmyear'].values[0])
    std = float(frameta[col_mm + '_RMSE_selection'].values[0])
    rmse = float(frameta[col_mm + '_RMSE_full'].values[0])
    rmse2 = float(frameta[col_mm + '_RMSE_mmy_full'].values[0])
    #
    years_since_beginning = frame_esds['years_since_beginning'].values
    frame_esds['model'] = years_since_beginning * slope + intercept
    frame_esds = frame_esds.set_index('years_since_beginning').sort_index()
    frame_esds['years_since_beginning_dup'] = years_since_beginning
    #
    # get it center:
    # sel=frame_esds[frame_esds[col_outliers] == False]
    # central = float(sel[col_mm].mean())
    # centralmodel = float(sel['model'].mean())
    # frame_esds[col_mm] = frame_esds[col_mm] - central
    # frame_esds['model'] = frame_esds['model'] - centralmodel
    #
    for ab in ['A', 'B']:
        frame_esds_ab = frame_esds[frame_esds['S1AorB'] == ab]
        if frame_esds_ab.empty:
            continue
        if ab == 'A':
            symbol = 'c'
        elif ab == 'B':
            symbol = 't'
        # print('first outliers')
        sel = frame_esds_ab[frame_esds_ab[col_outliers] == True]
        x = sel['epochdate'].values
        y = sel[col_mm].values
        #
        # fig.plot(x=x, y=y, style="x"+col_size, pen="0.1p,"+col_color+'3')
        # fig.plot(x=x, y=y, style="p0.08c", pen="0.1p,"+col_color)
        fig.plot(x=x, y=y, style=symbol + "0.08c", pen="thin," + col_color + '3')
        # fig.plot(x=x, y=y, style="x"+col_size, pen="1p,"+col_color+'3')
        #
        # print('main points')
        sel = frame_esds_ab[frame_esds_ab[col_outliers] == False]
        x = sel['epochdate'].values
        y = sel[col_mm].values
        if ab == 'A':
            fig.plot(x=x, y=y, style=symbol + col_size, fill=col_color + '3', pen='0.1p,' + col_color + '4',
                     label=col_label)
        else:
            fig.plot(x=x, y=y, style=symbol + col_size, fill=col_color + '3', pen='0.1p,' + col_color + '4')
    #
    # print('slope on notide')
    sel = frame_esds[frame_esds[col_outliers] == False]
    x = sel['epochdate'].values
    y = sel['model'].values
    title = 'vel: {0} +-{1} mm/y (95%)'.format(round(slope), round(2 * rmse2))
    # title = 'vel: {0} +-{1} mm/y (95%)'.format(round(slope),round(rmse))
    fig.plot(x=[x[0], x[-1]], y=[y[0], y[-1]], pen="3p," + col_color + '4', label=title)
    #
    if plotstd:
        fig.plot(x=[x[0], x[-1]], y=[y[0] - 2 * std, y[-1] - 2 * std], pen="0.2p," + col_color + '4,-')
        fig.plot(x=[x[0], x[-1]], y=[y[0] + 2 * std, y[-1] + 2 * std], pen="0.2p," + col_color + '4,-')
    #
    if additrf:
        slope_itrf = float(frameta['slope_plates_vel_azi_itrf2014'].values[0])
        title = 'ref vel (ITRF): {0} mm/y'.format(round(slope_itrf))
        y1 = y[-1]
        x0yrs = sel['years_since_beginning_dup'].values[0]
        x1yrs = sel['years_since_beginning_dup'].values[-1]
        yrsdiff = x1yrs - x0yrs
        y0 = y1 - yrsdiff * slope_itrf
        fig.plot(x=[x[0], x[-1]], y=[y0, y1], pen="2p,grey", label=title)
    return fig


######### generate gmt plot (arrows)

def plot_decomposed(dec, col = 'noTI', saveit = False):
    median_diff = dec['VEL_E_'+col].median() - dec['ITRF_E'].median()
    
    import pygmt
    x = dec.centroid_lon.values
    y = dec.centroid_lat.values
    cpxEN = dec.ITRF_E.values + 1j*dec.ITRF_N.values
    cpxEN = dec['VEL_E_'+col].values - median_diff + 1j*dec['VEL_N_'+col].values
    
    direction = np.degrees(np.angle(cpxEN))
    length = np.abs(cpxEN) # in mm/year
    
    fig = pygmt.Figure()
    fig.coast(
        region=[-180, 180, -50, 50],
        projection="M0/0/12c",
        #projection="T35/10c",
        frame=True,
        borders=False,
        shorelines="0.01p,black",
        area_thresh=4000,
        land="lightbrown",
        water="lightblue",
    )
    
    fig.plot(
        x=x,
        y=y,
        style="v0.025c+ea+bc",
        #style="v0.05c+ea+bc",
        direction=[direction, length/100],
        pen="0.1p",
        color="red3",
    )
    if saveit:
        fig.savefig(col+'_nomedian.png', dpi = 320)
    fig.show()



# 2022-05-23 - figure for GRL - POD change effect
'''
dist=framespd[framespd['opass']=='A'].podoff5.rename('ascending frames')*-1
dist2=framespd[framespd['opass']=='D'].podoff5.rename('descending frames')*-1
dist3=framespd.podoff5.rename('all frames')*-1




plt.style.use('seaborn-white')
from palettable.colorbrewer.qualitative import Set2_7
colors = Set2_7.mpl_colors

params = {
   'axes.labelsize': 8,
   'font.size': 9,
   'legend.fontsize': 8,
   'figure.dpi' : 200,
   #'xtick.labelsize': 10,
   #'ytick.labelsize': 10,
   #'text.usetex': False,
   'figure.figsize': [3, 4]
   # 'figure.figsize': [6, 4]
   }
plt.rcParams.update(params)


fig, ax = plt.subplots()
#fig = plt.figure(figsize=(8,4), dpi=100)
#ax.set_xlim(-100,100)
bins1 = 12
bins = bins1*2
ax = dist.plot.hist(density=False, legend=True, color=colors[2], histtype='stepfilled',
                ax=ax, bins=bins1, grid=True, alpha=0.5, edgecolor='k', 
                label='ascending frames')
ax = dist2.plot.hist(density=False, legend=True, color=colors[4], histtype='stepfilled',
                ax=ax, bins=bins1, grid=True, alpha=0.5, edgecolor='k',
                label='descending frames')
#dist3.plot.hist(density=True, legend=True,
#                ax=ax, bins=60, grid=True, edgecolor='k', alpha=0.05, label='bagr')
ax = dist3.plot.hist(density=False, legend=True, color=colors[4], histtype='step',
                ax=ax, bins=bins, grid=True, alpha=0.85, edgecolor='red',linewidth=2,
                label='all frames')

               #title='Density histogram of $\Delta(\sigma_{post}-\sigma_{pre})$', 

median_both = dist3.median()
mean_both = dist3.mean()
print('median values are:')
print(dist.median())
print(dist2.median())
print(dist3.median())

print('mean values are:')
print(dist.mean())
print(dist2.mean())
print(dist3.mean())

print('count:')
print(dist3.count())



ax.set_xlim(-120,80)
ax.set_xlim(-100,100)
#ax.set_ylim(0,0.4)
#ax.set_title("Estimates of $\Delta a$ offset due to POD change")
ax.set_title("$\Delta a$ offset due to orbits change")
ax.grid(True)

stderr_median=np.sqrt(dist3.var()/dist3.count())*2*1.253
ax.axvline(median_both, color='k', linewidth=1, alpha = 0.7, linestyle='dashed', label='$\overline{\delta \Delta a}$ ='+str(int(np.round(median_both)))+'$\pm$'+str(np.round(stderr_median,1))+' mm')
#ax.axvline(mean_both, color='k', linewidth=1, alpha = 0.7, linestyle='dashed')
ax.set_ylim(0,90)
min_ylim, max_ylim = ax.get_ylim()
#ax.text(median_both+2, max_ylim*0.939, '$\overline{\delta \Delta a}$',
#        fontsize=11)
#ax.text(median_both-14, max_ylim*0.939, '$\overline{\delta \Delta a}$',
#        fontsize=11)
#ax.axvline(mean_both, color='k', linewidth=1, alpha = 0.7, linestyle='dashed')
#ax.text(mean_both, max_ylim*0.95, ' $\mu_{\Delta,post}=$'+'{:.2f} mm'.format(mean_both),
#        fontsize=9)

#the daz_ARP = -39 mm
ax.axvline(-39, color='b', linewidth=1, alpha = 0.7, label='$\Delta a_{ARP}$ = -39 mm')
#ax.text(-39+2, max_ylim*0.945, '$daz_{ARP}$ = -39 mm', color='b', fontsize=11)
#ax.text(-39+1, max_ylim*0.941, '$\Delta a_{ARP}$', color='b', fontsize=11)
#ax.legend(facecolor='white', framealpha=1, loc='lower left')
legend = ax.legend(frameon = 1, loc='upper right')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_linewidth(0)
ax.set_xlabel('mm')
'''




#20210809 - figures for GRL article:

#figure 
def df_compare_new_orbits(esds, col = 'daz_mm_notide_noiono_grad_OK'):
    std_diffs = []
    #for frame in framespd['frame'].values:
    for frame, selected_frame_esds in esds.groupby('frame'):
        neworb = selected_frame_esds[selected_frame_esds['epochdate'] > pd.Timestamp('20200731')]
        oldorb = selected_frame_esds[selected_frame_esds['epochdate'] > pd.Timestamp('20190701')]
        oldorb = oldorb[oldorb['epochdate'] < pd.Timestamp('20200730')]
        if (not neworb.empty) and (not oldorb.empty):
            if len(neworb) >= 15:
                if not len(oldorb) < len(neworb):
                    oldorb = oldorb.tail(len(neworb))
                    std_old = oldorb[col].std()
                    std_new = neworb[col].std()
                    diff = std_new - std_old
                    if abs(diff) < 500:
                        std_diffs.append(diff)
                #return oldorb, neworb, std_old, std_new
                
            #number=len(neworb)
    std_diffs = np.array(std_diffs)
    # try:
    # pd.DataFrame(std_diffs).hist()
    return std_diffs



def figure_compare(esds):
    diffs2  = df_compare_new_orbits(esds)
    
    dist = pd.DataFrame(diffs2)
    dist = dist[abs(dist[0])<100]
    #workaround for the label.... in 2021...horrible
    dist = dist.rename(columns={0:'KDE'})
    fig, ax = plt.subplots()
    #dist.plot.hist(density=False, legend=False, ax=ax, bins=200)
    dist.plot.hist(density=True, legend=False, ax=ax, bins=40, grid=True, edgecolor='k', alpha=0.65)
    dist.plot.kde(ax=ax, legend=True, title='Density histogram of $\Delta(\sigma_{post}-\sigma_{pre})$', label='bagr', color='red')
    rmse_median = dist.median().values[0]
    plt.axvline(rmse_median, color='k', linestyle='dashed', linewidth=1)
    plt.rcParams["figure.figsize"] = (7,6)
    min_ylim, max_ylim = plt.ylim()
    plt.text(rmse_median*0.6, max_ylim*0.95, '$M_\Delta=${:.2f} mm'.format(rmse_median))
    
    ax.set_xlabel("mm")
    ax.set_xlim(-70, 70)
    fig.savefig('fig_hist.png', dpi=150)

    print(dist.median())



# 2022-06-15 - final fig for AHB plot:
'''
import pygmt
import numpy as np
#dec = a.copy(deep=True)
#dec = bagr.copy(deep=True)
strTI = 'noTI'
#dec = dec2

dec['VEL_E'] = dec['VEL_E_'+strTI]
dec['VEL_N'] = dec['VEL_N_'+strTI]
dec['RMSE_VEL_N'] = dec['RMSE_VEL_N_'+strTI]
dec['RMSE_VEL_E'] = dec['RMSE_VEL_E_'+strTI]

dec = dec[abs(dec['VEL_E']) < 100]
dec = dec[dec['RMSE_VEL_E'] < 35]

median_diff_E = dec['VEL_E'].median() - dec['ITRF_E'].median()
median_diff_N = dec['VEL_N'].median() - dec['ITRF_N'].median()
print('median E: '+str(median_diff_E))
print('median N: '+str(median_diff_N))
# check without median correction
median_diff_E = 0
median_diff_N = 0

df = pd.DataFrame(
    data={
        "x": dec.centroid_lon.values,
        "y": dec.centroid_lat.values,
        "east_velocity": dec.VEL_E.values, # [0, 3, 4, 6, -6, 6],
        "north_velocity": dec.VEL_N.values, # [0, 3, 6, 4, 4, -4],
        "east_sigma": np.zeros_like(dec.centroid_lon.values), # [4, 0, 4, 6, 6, 6],
        "north_sigma": np.zeros_like(dec.centroid_lon.values), # [6, 0, 6, 4, 4, 4],
        "correlation_EN": np.zeros_like(dec.centroid_lon.values), # [0.5, 0.5, 0.5, 0.5, -0.5, -0.5],
        "SITE": dec.index.values, #["0x0", "3x3", "4x6", "6x4", "-6x4", "6x-4"],
    }
)

df_unc = pd.DataFrame(
    data={
        "x": dec.centroid_lon.values,
        "y": dec.centroid_lat.values,
        "east_velocity": dec.VEL_E.values, # [0, 3, 4, 6, -6, 6],
        "north_velocity": dec.VEL_N.values, # [0, 3, 6, 4, 4, -4],
        "east_sigma": 2*dec.RMSE_VEL_E.values, # [4, 0, 4, 6, 6, 6],
        "north_sigma": 2*dec.RMSE_VEL_N.values, # [6, 0, 6, 4, 4, 4],
        "correlation_EN": np.zeros_like(dec.centroid_lon.values), # [0.5, 0.5, 0.5, 0.5, -0.5, -0.5],
        "SITE": dec.index.values, #["0x0", "3x3", "4x6", "6x4", "-6x4", "6x-4"],
    }
)

df_noT = pd.DataFrame(
    data={
        "x": dec.centroid_lon.values,
        "y": dec.centroid_lat.values,
        "east_velocity": dec.VEL_E_noT.values, # [0, 3, 4, 6, -6, 6],
        "north_velocity": dec.VEL_N_noT.values, # [0, 3, 6, 4, 4, -4],
        "east_sigma": np.zeros_like(dec.centroid_lon.values), #2*dec.RMSE_VEL_E_noT.values, # [4, 0, 4, 6, 6, 6],
        "north_sigma": np.zeros_like(dec.centroid_lon.values), #2*dec.RMSE_VEL_N_noT.values, # [6, 0, 6, 4, 4, 4],
        "correlation_EN": np.zeros_like(dec.centroid_lon.values), # [0.5, 0.5, 0.5, 0.5, -0.5, -0.5],
        "SITE": dec.index.values, #["0x0", "3x3", "4x6", "6x4", "-6x4", "6x-4"],
    }
)

df_itrf = pd.DataFrame(
    data={
        "x": dec.centroid_lon.values,
        "y": dec.centroid_lat.values,
        "east_velocity": dec.ITRF_E.values, # [0, 3, 4, 6, -6, 6],
        "north_velocity": dec.ITRF_N.values, # [0, 3, 6, 4, 4, -4],
        "east_sigma": np.zeros_like(dec.centroid_lon.values),
        "north_sigma": np.zeros_like(dec.centroid_lon.values),
        "correlation_EN": np.zeros_like(dec.centroid_lon.values), # [0.5, 0.5, 0.5, 0.5, -0.5, -0.5],
        "SITE": dec.index.values, #["0x0", "3x3", "4x6", "6x4", "-6x4", "6x-4"],
    }
)

df['east_velocity'] = df['east_velocity'] - median_diff_E
df['north_velocity'] = df['north_velocity'] - median_diff_N
x = dec.centroid_lon.values
y = dec.centroid_lat.values
cpxEN = dec.ITRF_E.values + 1j*dec.ITRF_N.values
#cpxEN = dec['VEL_E'].values - median_diff + 1j*dec['VEL_N'].values

direction = np.degrees(np.angle(cpxEN))
length = np.abs(cpxEN) # in mm/year

region=[25, 113, 22, 45]


fig = pygmt.Figure()
pygmt.config(MAP_FRAME_TYPE="plain")
fig.coast(
    region=region,
    #region=[-180, 180, -50, 50],
    projection="M0/0/30c",
    #projection="T35/10c",
    #frame=True,
    #frame=["WSne", "2g2f"],
    frame=["WNse", "5f", "a5f" ],
    #frame='a.5f.25WNse',
    borders=False,
    shorelines="0.01p,black",
    area_thresh=4000,
    land="lightgray",
    water="lightblue1",
)

# plot faults to the background
fig.plot(data=faults, pen="0.1p,darkgray", label="faults")

# first plot uncertainties (background)
fig.velo(
    data=df_unc,
    region=region,
    pen="0.4p,blue",
    uncertaintycolor="whitesmoke",
    #label='uncertainty',
    transparency=60,
    line=True,
    spec="e0.02/0.39/18",
    #projection="x0.8c",
    vector="0.25c+p1p+e+gblue",
)

# then plot arrows with only tides corrected vels
fig.velo(
    data=df_noT,
    region=region,
    pen="0.1p,black,-",
    #uncertaintycolor="whitesmoke",
    line=True,
    #label='SETcorrected',
    spec="e0.02/0.39/18",
    #projection="x0.8c",
    vector="0.1c+e+gblack",
)

# then the ITRF2014 PMM
fig.velo(
    data=df_itrf,
    region=region,
    pen="0.4p,red",
    #label='ITRF2014',
    #uncertaintycolor="lightblue1",
    line=True,
    spec="e0.02/0.39/18",
    #projection="x0.8c",
    vector="0.25c+e+gred",
)

# finally (blue) final velocities
fig.velo(
    data=df,
    region=region,
    pen="0.4p,blue",
    #uncertaintycolor="whitesmoke",
    line=True,
    #label='finalvel',
    spec="e0.02/0.39/18",
    #projection="x0.8c",
    #vector="0.25c+p1.5p+e+gblue",
    vector="0.25c+e+gblue",
)

#fig.plot(
#    x=x,
#    y=y,
#    #style="v0.1025c+ea+bc",
#    #style="v0.05c+ea+bc",
#    direction=[direction, length/100],
#    pen="0.3p",
#    color="blue3",
#)

fig.savefig('AHBa2b_'+strTI+'.png', dpi = 200)
'''



'''
not updated (but working):
import pandas as pd
import os

def get_pod_offset(dazes, years, thresyears = 4, minsamples = 15):
    r2 = np.ones_like(years)
    r3 = np.ones_like(years)
    r2[years >= thresyears] = 0
    r3[years < thresyears] = 0

    # if too small dataset, cancelling
    if np.sum(r3) < minsamples or np.sum(r2) < minsamples:
        return np.nan
    
    # stack them to get A
    A = np.vstack((years,r2,r3)).T
    model = np.linalg.lstsq(A,dazes, rcond=False)[0]
    offset = model[2] - model[1]
    return offset


esds = pd.read_csv('esds_updated_iono_ok_20220324.csv')
esds['epochdate'] = esds.apply(lambda x : pd.to_datetime(str(x.epoch)).date(), axis=1)

firstdatei = pd.Timestamp('2016-07-30')
poddate = pd.Timestamp('2020-07-30')
thresyears = aaa.days / 365.25

esds['years_since_beginning'] = esds['epochdate'] - firstdatei.date()
esds['years_since_beginning'] = esds['years_since_beginning'].apply(lambda x: float(x.days)/365.25)
esds = esds[esds['epochdate'] > firstdatei.date()]
frames = esds['frame'].unique()

offsets = []
for frame in frames:
    dazes = esds[esds['frame'] == frame].daz_mm_notide_noiono_ok.values
    years = esds[esds['frame'] == frame].years_since_beginning.values
    offsets.append(get_pod_offset(dazes, years))

offsets = np.array(offsets)
# filtering - offset would not be larger than... 150 mm?
offsets = offsets[np.abs(offsets)<150]

meanoffset = np.nanmean(offsets)
stderr = np.nanstd(offsets, ddof=3)/np.sqrt(len(offsets))
'''










'''
this is figure for AHB E,N plot:

#col = 'noTI' #noTI or noT for notide+noiono and notide
#median_diff = dec['VEL_E_'+col].median() - dec['ITRF_E'].median()

import numpy as np
#dec = a.copy(deep=True)
dec = bagr.copy(deep=True)
strTI = 'noT'


dec['VEL_E'] = dec['VEL_E_'+strTI]
dec['VEL_N'] = dec['VEL_N_'+strTI]
dec['RMSE_VEL_N'] = dec['RMSE_VEL_N_'+strTI]
dec['RMSE_VEL_E'] = dec['RMSE_VEL_E_'+strTI]

dec = dec[abs(dec['VEL_E']) < 100]
dec = dec[dec['RMSE_VEL_E'] < 35]

median_diff_E = dec['VEL_E'].median() - dec['ITRF_E'].median()
median_diff_N = dec['VEL_N'].median() - dec['ITRF_N'].median()
print('median E: '+str(median_diff_E))
print('median N: '+str(median_diff_N))
# check without median correction
median_diff_E = 0
median_diff_N = 0

df = pd.DataFrame(
    data={
        "x": dec.centroid_lon.values,
        "y": dec.centroid_lat.values,
        "east_velocity": dec.VEL_E.values, # [0, 3, 4, 6, -6, 6],
        "north_velocity": dec.VEL_N.values, # [0, 3, 6, 4, 4, -4],
        "east_sigma": 2*dec.RMSE_VEL_E.values, # [4, 0, 4, 6, 6, 6],
        "north_sigma": 2*dec.RMSE_VEL_N.values, # [6, 0, 6, 4, 4, 4],
        "correlation_EN": np.zeros_like(dec.centroid_lon.values), # [0.5, 0.5, 0.5, 0.5, -0.5, -0.5],
        "SITE": dec.index.values, #["0x0", "3x3", "4x6", "6x4", "-6x4", "6x-4"],
    }
)

df_itrf = pd.DataFrame(
    data={
        "x": dec.centroid_lon.values,
        "y": dec.centroid_lat.values,
        "east_velocity": dec.ITRF_E.values, # [0, 3, 4, 6, -6, 6],
        "north_velocity": dec.ITRF_N.values, # [0, 3, 6, 4, 4, -4],
        "east_sigma": np.zeros_like(dec.centroid_lon.values),
        "north_sigma": np.zeros_like(dec.centroid_lon.values),
        "correlation_EN": np.zeros_like(dec.centroid_lon.values), # [0.5, 0.5, 0.5, 0.5, -0.5, -0.5],
        "SITE": dec.index.values, #["0x0", "3x3", "4x6", "6x4", "-6x4", "6x-4"],
    }
)

df['east_velocity'] = df['east_velocity'] - median_diff_E
df['north_velocity'] = df['north_velocity'] - median_diff_N
x = dec.centroid_lon.values
y = dec.centroid_lat.values
cpxEN = dec.ITRF_E.values + 1j*dec.ITRF_N.values
#cpxEN = dec['VEL_E'].values - median_diff + 1j*dec['VEL_N'].values

direction = np.degrees(np.angle(cpxEN))
length = np.abs(cpxEN) # in mm/year

region=[25, 113, 22, 45]


fig = pygmt.Figure()
pygmt.config(MAP_FRAME_TYPE="plain")
fig.coast(
    region=region,
    #region=[-180, 180, -50, 50],
    projection="M0/0/30c",
    #projection="T35/10c",
    #frame=True,
    #frame=["WSne", "2g2f"],
    frame=["WNse", "5f", "a5f" ],
    #frame='a.5f.25WNse',
    borders=False,
    shorelines="0.01p,black",
    area_thresh=4000,
    land="lightgray",
    water="lightblue1",
)

fig.velo(
    data=df,
    region=region,
    pen="0.2p",
    uncertaintycolor="whitesmoke",
    line=True,
    spec="e0.02/0.39/18",
    #projection="x0.8c",
    vector="0.2c+p1p+e+gblack",
)


fig.velo(
    data=df_itrf,
    region=region,
    pen="0.1p,red",
    #uncertaintycolor="lightblue1",
    line=True,
    spec="e0.02/0.39/18",
    #projection="x0.8c",
    vector="0.1c+p0.5p+e+gred",
)


#fig.plot(
#    x=x,
#    y=y,
#    #style="v0.1025c+ea+bc",
#    #style="v0.05c+ea+bc",
#    direction=[direction, length/100],
#    pen="0.3p",
#    color="blue3",
#)

#fig.savefig('AHBa2b_'+strTI+'.png', dpi = 200)
fig.savefig(os.path.join('AHBa_MAY_'+strTI+'.pdf'), dpi = 320)
fig.show()

'''

def plot_daz_frame_licsar(frame, limit = 8000, newold=True):
    import daz_lib_licsar as dl
    dazes = dl.get_daz_frame(frame)
    msab = dl.fc.get_frame_master_s1ab(frame)
    mdatetime=dl.fc.get_master(frame, asdatetime=True)
    epochdates=dazes['epoch'].tolist()
    ABs = dl.flag_s1b(epochdates,mdatetime,msab,True)
    dazes['AB'] = ABs
    import numpy as np
    if not newold:
        toplotA = dazes[dazes['AB']=='A'].set_index('epoch').daz*14000
        toplotB = dazes[dazes['AB']=='B'].set_index('epoch').daz*14000
        toplotA=toplotA[np.abs(toplotA)<limit]
        toplotB=toplotB[np.abs(toplotB)<limit]
        toplotA.plot(title=frame, ylabel='$u_{az}$ [mm]', marker='o', linestyle='')#-.')
        toplotB.plot(title=frame, ylabel='$u_{az}$ [mm]', marker='o', linestyle='')#-.')
    else:
        azioffs=dl.get_azioffs_old_new_POD(frame)
        dazesep=dazes.set_index('epoch')
        dazesep['pod_diff_mm']=azioffs.set_index('epochdate')['pod_diff_azi_mm']
        dazesep['pod_diff_mm']=dazesep['pod_diff_mm'].fillna(0)
        toplotA = dazesep[dazesep['AB']=='A'].daz*14000
        toplotB = dazesep[dazesep['AB']=='B'].daz*14000
        toplotA=toplotA+dazesep[dazesep['AB']=='A'].pod_diff_mm
        toplotB=toplotB+dazesep[dazesep['AB']=='B'].pod_diff_mm
        toplotA=toplotA[np.abs(toplotA)<limit]
        toplotB=toplotB[np.abs(toplotB)<limit]
        title='azioff (PODdiff-corrected): '+frame
        toplotA.plot(title=title, ylabel='$u_{az}$ [mm]', marker='o', linestyle='')#-.')
        toplotB.plot(title=title, ylabel='$u_{az}$ [mm]', marker='o', linestyle='')#-.')
