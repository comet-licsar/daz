#!/usr/bin/env python3

# decompose
##########
import geopandas
import shapely
from daz_lib import *

#function for decomposition (by LS inversion)
def decompose_azi2NE(df, col = 'daz_mm_notide_noiono_grad'):
    velcol = 'slope_'+col+'_mmyear'
    rmscol = col+'_RMSE_mmy_full'
    if not velcol in df.columns:
        #print('workaround for a single column')
        velcol = col
        rmscol = col
    A = []
    d = []
    Qt = []
    for i, row in df.iterrows():
        heading = float(row['heading'])
        At = [np.sin(np.radians(heading)), np.cos(np.radians(heading))]
        Qt.append(1/row[rmscol]**2)
        #b.append(row[rmscol]**2)
        A.append(At)
        d.append(row[velcol])
    A = np.array(A)
    d = np.array(d)
    Q = np.zeros((len(d),len(d)))
    np.fill_diagonal(Q,Qt)
    lstsq = np.linalg.lstsq(A,d, rcond=None)
    # Qm will be variance for V,E and for V,N
    try:
        Qm = np.linalg.inv(A.transpose() @ Q @ A)
    except:
        print('matrix for frames listed below is singular! returning 999999999999 for var')
        print(df['frame'].values)
        Qm = np.zeros([2,2])
        np.fill_diagonal(Qm,999999999)
    Qm = np.abs(Qm.diagonal())
    V_E = lstsq[0][0]
    V_N = lstsq[0][1]
    RMSE_E = np.sqrt(Qm[0])
    RMSE_N = np.sqrt(Qm[1])
    return pd.DataFrame({'V_N': pd.Series(V_N),
                         'V_E': pd.Series(V_E),
                         'RMSE_E': pd.Series(RMSE_E),
                         'RMSE_N': pd.Series(RMSE_N),
                         })


# get ITRF N, E values
def get_itrf_EN(df):
    itrfs_N = []
    itrfs_E = []
    itrfs_rms_N = []
    itrfs_rms_E = []
    iii = 0
    fullcount = len(df)
    for ind, row in df.iterrows():
        iii = iii+1
        print('getting ITRF for {0}/{1} cells'.format(iii, fullcount))
        clon = row['centroid_lon']
        clat = row['centroid_lat']
        # use a median over 'whole' frame:
        Es = []
        Ns = []
        for i in range(round(clon*10-23.4/2),round(clon*10+23.4/2)+1,5):
            lon = i/10
            for j in range(round(clat*10-23.4/2),round(clat*10+23.4/2)+1,5):
                lat = j/10
                try:
                    E, N = get_ITRF_ENU(lat, lon)
                    Es.append(E)
                    Ns.append(N)
                    #itrfs.append(EN2azi(N, E, heading))
                except:
                    print('connection error')
        itrfs_E.append(np.mean(Es))
        itrfs_N.append(np.mean(Ns))
        itrfs_rms_E.append(np.std(Es,ddof=1))
        itrfs_rms_N.append(np.std(Ns,ddof=1))
    df['ITRF_N'] = itrfs_N
    df['ITRF_E'] = itrfs_E
    df['ITRF_RMSE_E'] = itrfs_rms_E
    df['ITRF_RMSE_N'] = itrfs_rms_N
    return df
    




crs = "EPSG:4326"

framespd['opass'] = framespd['frame'].str[3]

gdf = geopandas.GeoDataFrame(framespd, 
            geometry=geopandas.points_from_xy(framespd.center_lon, framespd.center_lat),
            crs=crs)

# establish a grid
# projection of the grid
# total area for the grid
xmin, ymin, xmax, ymax= gdf.total_bounds
cell_size = 2.25  # this is some ~250x250 km

# create the cells in a loop
grid_cells = []
centroid_lon = []
centroid_lat = []
for x0 in np.arange(xmin, xmax+cell_size, cell_size ):
    for y0 in np.arange(ymin, ymax+cell_size, cell_size):
        # bounds
        x1 = x0-cell_size
        y1 = y0+cell_size
        grid_cells.append( shapely.geometry.box(x0, y0, x1, y1)  )
        centroid_lon.append(x0)
        centroid_lat.append(y0)



grid = geopandas.GeoDataFrame(grid_cells, columns=['geometry'], 
                                 crs=crs)

grid['centroid_lon'] = centroid_lon
grid['centroid_lat'] = centroid_lat

# merge framespd and the grid
merged = geopandas.sjoin(gdf, grid, how='left', op='within')

gridgrouped = merged.groupby('index_right')
gridagg = gridgrouped.agg(count=('opass', 'count'),
                        opass=('opass', list),
                        centroid_lon=('centroid_lon', 'mean'),
                        centroid_lat=('centroid_lat', 'mean'))


gridagg = gridagg[gridagg['count'] > 1]
gridagg = gridagg[gridagg['opass'].str.contains('D', regex=False) & gridagg['opass'].str.contains('A', regex=False)]

#now the gridagg contains only A+D cells
# 1. reduce merged and grouped:
for i in merged.index.values:
    if merged.loc[i].index_right not in gridagg.index:
        merged = merged.drop(i)


gridgrouped = merged.groupby('index_right')

# 2. now do the decomposition
decomposed = gridgrouped.apply(decompose_azi2NE, 'daz_mm_notide_noiono_grad')
gridagg['VEL_N_noTI'] = decomposed['V_N'].values
gridagg['VEL_E_noTI'] = decomposed['V_E'].values
gridagg['RMSE_VEL_N_noTI'] = decomposed['RMSE_N'].values
gridagg['RMSE_VEL_E_noTI'] = decomposed['RMSE_E'].values

decomposed = gridgrouped.apply(decompose_azi2NE, 'daz_mm_notide')
gridagg['VEL_N_noT'] = decomposed['V_N'].values
gridagg['VEL_E_noT'] = decomposed['V_E'].values
gridagg['RMSE_VEL_N_noT'] = decomposed['RMSE_N'].values
gridagg['RMSE_VEL_E_noT'] = decomposed['RMSE_E'].values

'''
decomposed = gridgrouped.apply(decompose_azi2NE, 'daz_mm')
gridagg['VEL_N'] = decomposed['V_N'].values
gridagg['VEL_E'] = decomposed['V_E'].values
gridagg['RMSE_VEL_N'] = decomposed['RMSE_N'].values
gridagg['RMSE_VEL_E'] = decomposed['RMSE_E'].values
'''

gridagg = gridagg.dropna()

gridagg = get_itrf_EN(gridagg)

gridagg.to_csv('decomposed_20211012_noroll.csv')

