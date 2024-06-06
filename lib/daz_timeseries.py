#!/usr/bin/env python3

from daz_lib import *


#######################################
# step 5 - estimate velocities (and remove outliers etc.)
###################

def get_rmse(y, y_pred, ddof=2):
    return np.sqrt( np.sum(np.square(y - y_pred)) / (len(y)-ddof) )


### S1AB offset calculation:
def estimate_s1ab(frame_esds, col = 'daz_mm_notide_noiono', rmsiter = 50, printout = True):
    #epochdates = frame_esds['epochdate'].values
    isB = (frame_esds.S1AorB == 'B').values * 1
    if (np.sum(isB) < 20 and len(isB[isB == 0]) < 20) or (np.sum(isB) < 10):
        isB = isB * 0
        # print('cancelling for cAB')
    years = frame_esds.years_since_beginning.values
    dazes = frame_esds[col].values
    A = np.vstack((years,np.ones_like(years),isB)).T
    #res = model_filter(A, dazes, iters=rmsiter,years_to_pod=years_to_pod)
    res = model_filter_v2(A, dazes, iters=rmsiter, target_rmse = 30, printout = printout)
    model=res[0]
    stderr=res[1]
    v = model[0]
    c = model[1]
    c_AB = model[2]
    # now what to return:
    return v,c,stderr,c_AB


def estimate_s1ab_allframes(esds, framespd, col = 'daz_mm_notide_noiono', rmsiter = 50):
    lenframes = len(framespd['frame'])
    framespd['S1AB_offset'] = 0.0
    framespd['slope_daz_rmseiter_mmyear']=0.0
    framespd['intercept_daz_rmseiter_mmyear'] = 0.0
    framespd['stderr_daz_rmseiter_mm'] = 0.0
    i = 0
    for frame in framespd['frame']:
        i=i+1
        print('  Running for {0:6}/{1:6}th frame...'.format(i, lenframes), flush=True, end='\r')
        frameta = framespd[framespd['frame']==frame].copy()
        selected_frame_esds = esds[esds['frame'] == frame].copy()
        v,c,stderr,c_AB = estimate_s1ab(selected_frame_esds, col = col, rmsiter = rmsiter, printout = False)
        frameta['slope_daz_rmseiter_mmyear'] = v
        frameta['S1AB_offset'] = c_AB
        frameta['intercept_daz_rmseiter_mmyear']=c
        frameta['stderr_daz_rmseiter_mm'] = stderr
        framespd.update(frameta)
    return framespd

def correct_s1ab(esds, framespd, cols = ['daz_mm', 'daz_mm_notide', 'daz_mm_notide_noiono'], stderr_thres = 100):
    '''
    Will apply S1AB offset to given columns in esds pd, in case the stderr_daz_rmseiter_mm is lower than stderr_thres
    '''
    if 'S1AB_offset' not in framespd:
        print('ERROR, S1AB_offset not in framespd, cancelling')
        return esds, framespd
    i=0
    lenframes = len(framespd['frame'])
    for frame in framespd['frame']:
        i=i+1
        print('  Running for {0:6}/{1:6}th frame...'.format(i, lenframes), flush=True, end='\r')
        frameta = framespd[framespd['frame'] == frame].copy()
        selected_frame_esds = esds[esds['frame'] == frame].copy()
        selected_frame_esds = selected_frame_esds[selected_frame_esds['S1AorB']=='B']  #[framespd['stderr_daz_rmseiter_mm']>0]
        s1aboff = float(frameta['S1AB_offset'].values[0])
        if float(frameta['stderr_daz_rmseiter_mm'].values[0]) < stderr_thres:
            selected_frame_esds[cols] = selected_frame_esds[cols] - s1aboff
        esds.update(selected_frame_esds)
    return esds, framespd

# reduced version, 2024
def model_filter_v2(A, y, limrms=3, iters=2, target_rmse = 30, full_stderr = False, weighted = False, printout = True):
    '''
    limrms is 'how many RMSEs should be used to remove outliers'
    '''
    model = np.linalg.lstsq(A,y, rcond=False)[0]
    ddof = A.shape[1]
    y_pred = np.sum(A*model,axis=1)
    rmse = get_rmse(y, y_pred, ddof=ddof)
    oldrmse = rmse
    oldlen = len(y)
    for i in range(iters):
        count = len(y)
        # reducing dataset
        sel = np.abs(y_pred - y)<(limrms*rmse)
        y = y[sel]
        if len(y) == count:
            sel = np.abs(y_pred - y) < ((limrms - 1) * rmse)
            y = y[sel]
        if len(y)==count:
            break
        A = A[sel]
        # second iteration (only)
        model = np.linalg.lstsq(A,y, rcond=False)[0]
        y_pred = np.sum(A*model,axis=1)
        rmse = get_rmse(y, y_pred, ddof=ddof)
        #print('vel: '+str(model[0]) +'+-'+str(rmse))
        #print('len: '+str(len(y)))
        if rmse < target_rmse:
            break
    if printout:
        print('RMSE diff: '+str(oldrmse)+' -> '+str(rmse)+' ('+str(oldlen)+'/'+str(len(y))+' samples)')
    if weighted:
        print('not done, as perhaps not needed')
    if full_stderr:
        #properly calculate std errors of model parameters
        # here just assuming sigma of each observation based on expected meas. limits (7 mm processing, 10 mm rmse of satellite position)
        sigma = np.sqrt(7*7 + 100)
        try:
            stderr = get_model_sigma(sigma, A)
        except:
            stderr = np.nan
    else:
        stderr = rmse
    return model, stderr


def model_filter(A, y, limrms=3, iters=2):
    '''
    limrms is 'how many RMSEs should be used to remove outliers'
    '''
    model = np.linalg.lstsq(A,y, rcond=False)[0]
    ddof = A.shape[1]
    y_pred = np.sum(A*model,axis=1)
    rmse = get_rmse(y, y_pred, ddof=ddof)
    for i in range(iters):
        count = len(y)
        # reducing dataset
        sel = np.abs(y_pred - y)<limrms*rmse
        y = y[sel]
        if len(y)==count:
            break
        A = A[sel]
        # second iteration (only)
        model = np.linalg.lstsq(A,y, rcond=False)[0]
        y_pred = np.sum(A*model,axis=1)
        rmse = get_rmse(y, y_pred, ddof=ddof)
    stderr = np.sqrt(rmse**2/len(y))
    return model, stderr


def get_stdvel(rmse, tmatrix):
    """ Gets standard deviation of velocity by error propagation theory.
    
    Args:
        rmse (np.array or float): rmse (or std) of the solution (or of each measurement - better)
        tmatrix (np.array): transposed matrix of samples (time)
    
    Note: written when not fully understanding this..
    rmse should be of each measurement. Also, better to also include reduced chi-squared statistic Chi^2 that is
    r ... residuals between dependent variable and predicted value
    S = np.transpose(r) * W * r
    Chi^2 = S / (n-ddof)
    or we can write that:
    Chi^2 = np.sum((y-y_i)^2/variance) / (n-ddof)
    
    And for further studies: just see https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
    so chi^2 ~ 1 means appropriate fit
    """
    count = tmatrix.shape[0]
    #add column of ones to the tmatrix:
    cones = np.ones(tmatrix.shape)
    A = np.append(tmatrix, cones, axis = 1)
    #
    #now add the rmse to the diagonal of var-covar matrix, already invert it for later
    Qd = np.zeros((count,count)) #,len(d)))
    np.fill_diagonal(Qd,1/rmse**2)
    #
    # do Qm = (G.T Qd^-1 G)^-1  ---- Qm = [var_vel, var_intercept] (if m = [vel, intercept])
    Qm = np.linalg.inv(A.transpose() @ Qd @ A)
    var_vel = Qm.diagonal()[0]
    #var_intercept = Qm.diagonal()[1]
    STD_vel = np.sqrt(var_vel)
    #STD_intercept = np.sqrt(var_intercept)
    return STD_vel


def dates2ordinal(x):
    toord = lambda x: x.toordinal()
    vfunc = np.vectorize(toord)
    return vfunc(x)


def fit_huber(x, y, alpha = 1, epsilon=1.2, outinyears = True):
    '''
    just a Huber regression helper function - in case input x is in date(time),
    you may set either the output will be in years (e.g. mm/year) or if not, in days.
    '''
    typex = str(type(x[0])).split("'")[1].split('.')[0]
    xisdate = False
    if typex == 'datetime' or typex == 'pandas':
        xisdate = True
        #it is datetime, so converting to ordinals
        x = dates2ordinal(x)
        x = x - x[0]
    X = np.array([x]).transpose()
    huber = HuberRegressor(alpha = alpha, epsilon = epsilon)
    huber.fit(X, y)
    slope = huber.coef_[0]
    intercept = huber.intercept_
    y_pred = huber.predict(X)
    rmse = get_rmse(y,y_pred)
    std_vel = get_stdvel(rmse, X)
    if xisdate and outinyears:
        slope = slope*365.25
        std_vel = std_vel*365.25
    return slope, intercept, std_vel, y_pred, huber.outliers_


def df_calculate_slopes(esdsin, framespdin, alpha = 2.5, eps = 1.5, bycol = 'daz_mm_notide', subset = True, roll_assist = False):
    esds = esdsin.copy(deep=True)
    esds = esds.dropna()
    framespd = framespdin.copy(deep=True)
    #must start with True to initialise all 'as outliers'
    esds['is_outlier_'+bycol] = True
    framespd[bycol+'_RMSE_selection'] = 1.0
    framespd[bycol+'_count_selection'] = 1
    framespd[bycol+'_RMSE_full'] = 1.0
    framespd['slope_'+bycol+'_mmyear'] = -999
    framespd['intercept_'+bycol+'_mmyear'] = -999
    #now calculate per frame
    for frame, group in esds.groupby('frame'):
        print(frame)
        frameta = framespd[framespd['frame'] == frame]
        if frameta.empty:
            print('frame data is empty, skipping')
            continue
        # make preselection
        grsel = group.copy(deep=True)
        grsel = grsel[np.isfinite(grsel[bycol])]
        # get the full dataset, as we will filter it further on
        #limiting the dataset here, as data after 2020-07-30 have way different PODs (v. 1.7)
        #grsel = grsel[grsel['epochdate'] < pd.Timestamp('2020-07-30')] 
        # 2021-10-12 - finally corrected, probably, using daz_ARP - so using full dataset now!!!
        if subset:
            #limiting the dataset here, as often data before mid-2016 are influenced by ionosphere (probably)
            grsel = grsel[grsel['epochdate'] > dt.date(2016,3,1)] #pd.Timestamp('2016-03-01')]  #[grsel['epoch'] < 20200601]
        if len(grsel) < 20:
            print('most of data removed - cancelling for this frame')
            continue
        correctback = np.median(grsel['daz_mm'])
        grsel['daz_mm'] = grsel['daz_mm'] - correctback #np.median(grsel['daz_mm'])
        grsel = grsel[grsel['daz_mm'] > -400][grsel['daz_mm'] < 400]
        if len(grsel) < 15:
            print('most of data removed - cancelling for this frame')
            continue
        # center dataset by mean (mainly to avoid need of intercept - itrf estimate)
        x = grsel['years_since_beginning'].values
        y = grsel[bycol].values
        #center around zero as origin
        meany = np.mean(y)
        meanx = np.mean(x)
        X = np.array([x]).transpose()
        if roll_assist:
            #
            # prepare rolling (median 5) mean-centered data:
            roll = pd.DataFrame(y, index=x)
            roll = roll.rolling(5).mean()
            y2 = roll[0]
            y2 = y2[y2.notna()]
            x2 = y2.index
            meanx2 = meanx + np.mean(x2)
            meany2 = meany + np.mean(y2)
            X2 = x2.values.reshape(len(x2),1)
        # prepare first iteration to huber (warm start) using itrf slope
        #print(frameta)
        colwarm = 'slope_plates_vel_azi_gps'
        if not colwarm in frameta:
            colwarm = 'slope_plates_vel_azi_itrf2014'
        itrfslope = frameta[colwarm].values[0]
        itrfintercept = np.mean(x)*itrfslope * -1 #+ np.median(y)
        xitrf = np.arange(-5,5)
        yitrf = xitrf * itrfslope + itrfintercept
        Xitrf = xitrf.reshape((10, 1))
        # do huber regression on the original selection
        huber = HuberRegressor(alpha = alpha, epsilon = eps, warm_start=True)
        huber.fit(Xitrf, yitrf)
        huber.fit(X, y)
        slope = huber.coef_[0]
        intercept = huber.intercept_ #??
        y_pred = huber.predict(X)
        y_pred_all = y_pred
        if roll_assist:
            # do additional huber regression on the rolling dataset
            huber_roll = HuberRegressor(alpha = alpha, epsilon = eps, warm_start=True)
            huber_roll.fit(Xitrf, yitrf)
            huber_roll.fit(X2, y2)
            slope_roll = huber_roll.coef_[0]
            #intercept_roll = (meanx2 * slope_roll + huber_roll.intercept_) + meany2
            intercept_roll = huber_roll.intercept_   # CHECK IT
            # to predict, use all dates that go to the pre-rolling dataset...
            y_pred_roll = huber_roll.predict(X)
            # if huber roll is closer to itrf, use this
            if abs(slope-itrfslope) > abs(slope_roll-itrfslope):
                print('using rolling value '+str(slope_roll)+' instead of '+str(slope))
                print('(itrf shows '+str(itrfslope)+' )')
                slope = slope_roll
                intercept = intercept_roll
                y_pred = y_pred_roll
        # calculate RMSE using all points:
        rmse_full = get_rmse(y,y_pred_all)
        std_vel_full = get_stdvel(rmse_full, X)
        # and use outliers from huber !!!!!! even if using huber_rolling !!!!
        grsel['is_outlier_'+bycol] = huber.outliers_
        x = grsel[grsel['is_outlier_'+bycol] == False]['years_since_beginning'].values
        X = np.array([x]).transpose()
        y = grsel[grsel['is_outlier_'+bycol] == False][bycol].values
        try:
            y_pred = huber.predict(X)
            rmse = get_rmse(y,y_pred)
        except:
            print('error getting RMSE (all points are outliers?), using NaN value')
            rmse = np.nan
        #
        framespd.at[frameta.index.values[0], 'slope_'+bycol+'_mmyear'] = slope
        framespd.at[frameta.index.values[0], 'intercept_'+bycol+'_mmyear'] = intercept + correctback
        # add the rest:
        #keeping the older approach, but perhaps not needed...?
        residuals = y - y_pred
        framespd.at[frameta.index.values[0], bycol+'_RMSE_selection'] = rmse
        framespd.at[frameta.index.values[0], bycol+'_count_selection'] = len(residuals)
        framespd.at[frameta.index.values[0], bycol+'_RMSE_full'] = rmse_full
        # the STD mmy is perhaps the best estimate - kudos to Andy Hooper for this error propagation calculation approach
        framespd.at[frameta.index.values[0], bycol+'_RMSE_mmy_full'] = std_vel_full
        #
        # adding std/rmse for velocity - error propagation - this would be ok, but we do not know exact/real shift - so using based on the model here!
        years = grsel.years_since_beginning.max() - grsel.years_since_beginning.min()
        shift_mm = slope * years
        # std vel is from http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf - 'should' be right, but perhaps it is not?
        rms_vel_full = abs((rmse_full/shift_mm) * slope)
        rms_vel_sel = abs((rmse/shift_mm) * slope)
        framespd.at[frameta.index.values[0], bycol+'_RMSE_mmy_full_error_multi'] = rms_vel_full
        esds.update(grsel['is_outlier_'+bycol])
    return esds, framespd

'''
for col in ['daz_mm_notide', 'daz_mm_notide_noiono_grad']:
     print('estimating velocities of '+col)
     esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = col, subset = False, roll_assist = roll_assist)
'''
