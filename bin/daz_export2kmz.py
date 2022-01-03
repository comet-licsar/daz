###### generate kmz with plots
######### update 2021-05-31: hopefully improved huber regression...:



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



export_esds2kml(framespd, esds, kmzfile = 'esds.final.oct.2021.kmz', overwrite = True, clean = False)


