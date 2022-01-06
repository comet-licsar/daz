#!/usr/bin/env python3

from daz_lib import *

#for visualisation and kml export
import hvplot
import holoviews as hv
from holoviews import opts
hv.extension('bokeh')
import simplekml


def plot_vel_esd(frame_esds, frameta, level1 = 'iono_grad', level2 = None, showitrf=True):
    # level1: 'tide', 'orig', 'iono_grad', 'iono_f2'
    global mindate
    global maxdate
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
        col_outliers = 'is_outlier_'+col_mm
        col_label = 'daz notide, noiono (grad)'
        colour = 'green'
    elif level1 == 'iono_grad_OK':
        col_mm = 'daz_mm_notide_noiono_grad_OK'
        col_outliers = 'is_outlier_'+col_mm
        col_label = 'daz notide, noiono (grad)'
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
            col_outliers = 'is_outlier_'+col_mm
            col_label = 'daz notide, noiono (grad)'
            colour = 'green'
        elif level2 == 'iono_grad_OK':
            col_mm = 'daz_mm_notide_noiono_grad_OK'
            col_outliers = 'is_outlier_'+col_mm
            col_label = 'daz notide, noiono (grad)'
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





def export_esds2kml(framespd, esds, kmzfile = 'esds.kmz', overwrite = False, clean = True):
    if os.path.exists('plots') or os.path.exists('doc.kml') or os.path.exists(kmzfile):
        if not overwrite:
            print('The kmz files already exist, please recheck')
            return False
        else:
            os.system('rm -r {} doc.kml plots'.format(kmzfile))
    os.mkdir('plots')
    kml = simplekml.Kml()
    print('generating plots')
    #this will generate plots
    for frame in framespd['frame']:
        frameta = framespd[framespd['frame']==frame]
        selected_frame_esds = esds[esds['frame'] == frame].copy()
        #frameplot = plot_vel_esd(selected_frame_esds, frameta, showtec = False)
        try:
            frameplot = plot_vel_esd(selected_frame_esds, frameta, level2 = 'iono_grad', level1 = 'tide', showitrf=True)
            hv.save(frameplot, 'plots/{}.png'.format(frame), dpi=100, fmt='png')
        except:
            print('error generating plot for frame '+frameta['frame'])
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


######### generate gmt plot (arrows)

def plot_decomposed(dec, col = 'noTI', saveit = False):
    median_diff = dec['VEL_E_'+col].median() - dec['ITRF_E'].median()
    
    import pygmt
    x = dec.centroid_lon.values
    y = dec.centroid_lat.values
    cpxEN = dec.ITRF_E.values + 1j*dec.ITRF_N.values
    cpxEN = dec['VEL_E_'+col].values - median_diff + 1j*dec.['VEL_N_'+col].values
    
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




#20210809 - figures for GRL article:

#figure 
def df_compare_new_orbits(esds):
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
                    std_old = oldorb['daz_mm_notide_noiono_grad_OK'].std()
                    std_new = neworb['daz_mm_notide_noiono_grad_OK'].std()
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

