import numpy as np
# import os
import datetime
from scipy import interpolate
import time
from numba import jit, prange
import aacgmv2 #install with: pip install aacgmv2

from Code.subfun_sunAlsoRises_zenith_numba import sunAlsoRises_zenith_numba

#===================================== MODEL CALLER =====================================
def adelphi_modeler(data_AMPERE, hemi, dateRange_curr, dateRange_dayNum_curr,
                    FLG_smoothing_Jr=True,FLG_smoothing_conductance=3,FLG_smoothing_currentClip=True,
                    FLG_smoothing_polarBoost=True,FLG_smoothing_equatorBoost=False,
                    FLG_smoothing_phi_connectMidnight=True,FLG_smoothing_phi_lat=True,
                    FLG_smoothing_phi_long=True,FLG_smooth_phiFix=True,
                    FLG_smoothing_efield=True,FLG_magicks=False,FLG_textOutputs=False):
        
    #========== Constants ==========
    munaught = 4*np.pi*1.e-7; #H/m == N/Amp^2
    dtr = np.pi/180.;
    re = 6300E3; # Earth radius in m

    sitelat = np.asarray([68.36,73.55,77.72,71.58,66.17,71.30,64.87,62.40,58.80,55.27,61.20,64.18]); #lat of 12 sites used to calc AU/AL/AE
    sitelong = np.asarray([18.82,80.57,104.28,129.,190.17,203.25,212.17,245.60,265.9,285.22,314.16,338.3]); #long of 12 sites used to calc AU/AL/AE
    
    #--- Time Calcs (ahead of time) ---
    time_cntStop = 86400//data_AMPERE['data rate']; # number of times to cover
    # time_sec = np.arange(0,86400,data_AMPERE['data rate']);#sec, time in day in seconds
    # time_hr = time_sec/3600; #hr, time in day in hrs
    time_sec = np.unique(np.int32(data_AMPERE['hour'])*3600 + np.int32(data_AMPERE['min'])*60 + np.int32(data_AMPERE['sec']));#sec, time in day in seconds
    time_hr = time_sec/3600; #hr, time in day in hrs
    # data_AMPERE['time'] = time_sec;
    
    start_yr = dateRange_curr[0]; #punpack
    start_mo = dateRange_curr[1];
    start_dy = dateRange_curr[2];
    
    nlat = data_AMPERE['data info']['lat num'];
    nlon = data_AMPERE['data info']['long num'];
    
    time_stepsAvail = data_AMPERE['data info']['num time steps'];
    if( time_cntStop > time_stepsAvail ):
        time_cntStop = time_stepsAvail; #allows for time_cntStop to be independent (if that was ever what you wanted!)
    #END IF
        
    #create zero arrays
    # #    ***************************************************************************
    # #    arrays for calculating potential
    # spb = np.zeros((24,50),dtype=np.float64);
    # shb = np.zeros((24,50),dtype=np.float64);
    # # phib = np.zeros((24,50),dtype=np.float64);
    # # ajrb = np.zeros((24,50),dtype=np.float64);
    # #    ***************************************************************************
    # #    arrays for plotting UT variations
    # uthr = np.zeros((time_stepsAvail),dtype=np.float64);
    # # enflux = np.zeros((time_stepsAvail),dtype=np.float64);
    # # dbmax = np.zeros((time_stepsAvail),dtype=np.float64);
    # # power = np.zeros((time_stepsAvail),dtype=np.float64);
    powerg = np.zeros((time_stepsAvail),dtype=np.float64);
    facamp = np.zeros((time_stepsAvail),dtype=np.float64);
    potdiff = np.zeros((time_stepsAvail),dtype=np.float64);
    jh = np.zeros((time_stepsAvail),dtype=np.float64);
    # # jemax = np.zeros((time_stepsAvail),dtype=np.float64);
    # # jemin = np.zeros((time_stepsAvail),dtype=np.float64);
    # # qje = np.zeros((time_stepsAvail),dtype=np.float64);
    dbh = np.zeros((12,time_stepsAvail),dtype=np.float64);
    # # jrpoker = np.zeros((time_stepsAvail),dtype=np.float64);
    # # ve720 = np.zeros((time_stepsAvail),dtype=np.float64);
    # # vn720 = np.zeros((time_stepsAvail),dtype=np.float64);
    au_sim = np.zeros((time_stepsAvail),dtype=np.float64);
    al_sim = np.zeros((time_stepsAvail),dtype=np.float64);
    aesim = np.zeros((time_stepsAvail),dtype=np.float64);
    # # cpcp = np.zeros((time_stepsAvail),dtype=np.float64);        # for plotting cross polar cap potential
    # # spvsut = np.zeros((time_stepsAvail),dtype=np.float64);    # pedersen conductance at selected location
    # # shvsut = np.zeros((time_stepsAvail),dtype=np.float64);    # hall conductance at selected location
    # # maxped = np.zeros((time_stepsAvail),dtype=np.float64);    # maximum north-south pedersen current
    # # rmsplt = np.zeros((time_stepsAvail),dtype=np.float64);    # for plotting the rms deviations between facin and facout
    # # auroral parameter arrays***************
    # # ef = np.zeros((24,50),dtype=np.float64); # energy flux
    # jouleheat = np.zeros((24,50),dtype=np.float64);      
    # jehall = np.zeros((24,50),dtype=np.float64); #eastward hall current                                                                      
    # jeast = np.zeros((24,50),dtype=np.float64);# eastward current
    # jnorth = np.zeros((24,50),dtype=np.float64); #northward current
    # jrcalc = np.zeros((24,50),dtype=np.float64); #calculated upward current
    # # these arrays are for storing the parameters for printing output files or plotting
    spout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    shout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    efout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    jouleout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    phiout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    enout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    eeout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    jeout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    jehallout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    jnout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    # # jrin = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    jrout = np.zeros((time_stepsAvail,24,50),dtype=np.float64);
    
    #     +++++++++++++++++++++++++++++++++++++++++++++++++
    #    set mlt and mlat range and increment***************************************
    thmax = 50.;
    delth = 1.;
    # delthkm = delth*109.;
    #nth = np.int64(1+(thmax-thmin)/delth);    #number of theta values at which quantities are calculated
    #    arrays that contain the angles in spherical coordinates
    thsn = delth*np.arange(50,dtype=np.float64);
    dph = 15*dtr;
    phmlt = 15*np.arange(24,dtype=np.float64);
    mag_lat = (90.-thmax)+delth*np.arange(50,dtype=np.float64);    #the magnetic latitudes corresponding to the theta array
    th50 = 40.+np.arange(50,dtype=np.float64);
    if (hemi.lower() == 'south') | (hemi.lower() == 'southern') :
        mag_lat = -mag_lat; #FOR SOUTHERN HEMISPHERE
        th50=-40.-np.arange(50,dtype=np.float64); #FOR SOUTHERN HEMISPHERE
    #END IF
    
    #precalc values as needed
    #--- area of chunks at each latitude band ---
    th50 = mag_lat;
    dth0=(1./360.)*2*np.pi; #for area
    dph0=(15./360.)*2*np.pi; #for area
    theta = th50*np.pi/180.; #theta is latitude from south to north in radians
    area = dth0*dph0*re*re*np.cos(theta); #area of cell in m2
    
    #--- AACGMv2 calcs (not numba-friendly) ---
    geo_lat = np.empty((time_sec.size,nlat,nlon));
    geo_long = np.empty((time_sec.size,nlat,nlon));
    
    timeTemp_hr = time_sec//3600; #hr component
    timeTemp_min = np.mod(time_sec,3600)//60; #minute component
    timeTemp_sec = np.int64(np.mod(np.mod(time_sec,3600),60)); #second component
    timeTemp_usec = np.int64(np.round((np.mod(np.mod(time_sec,3600),60) - np.int64(np.mod(np.mod(time_sec,3600),60)))*1E6)); #microsecond component with FP32 basic truncation protection
    mag_mltTemp = np.arange(0,24,nlon/24,dtype=np.float64); #"hrs", precalc MLT variable 0 to 23
    for jk in range(0,time_sec.size):
        dtime = datetime.datetime(start_yr, start_mo, start_dy, hour=timeTemp_hr[jk], minute=timeTemp_min[jk], second=timeTemp_sec[jk], microsecond=timeTemp_usec[jk]); #date time object for aacgmv2    
        for jj in range(0,nlon):
            mag_long = aacgmv2.convert_mlt(mag_mltTemp[jj], dtime, m2a=True); #convert to mag_long [this supports vectorized form but not worth using since gotta loop to make dtime anyway]
            for kk in range(0,nlat):
                [geo_lat[jk,kk,jj], geo_long[jk,kk,jj], _] = aacgmv2.convert_latlon(mag_lat[kk], mag_long, 150., dtime, method_code='A2G'); # see https://aacgmv2.readthedocs.io/en/latest/reference/aacgmv2.html#aacgmv2.wrapper.convert_latlon
            #END FOR jj (longitude grid)
        #END FOR kk (latitude grid)
    #END FOR jk (times)
    geo_long[geo_long < 0] += 360; #wrap around to 0 to 360 instead of -180 to 180
    
    #prep jrin just how it needs to be
    jrin = np.flip(data_AMPERE['Jr'].reshape(time_stepsAvail,24,50).astype(np.float64),axis=2);
    
    ayear = start_yr;
    month = start_mo;
    aday = start_dy;
    
    #================= key loop ===================================
    for ntime in range(0,time_cntStop):
                        
        #----- Get the FAC data -----
        # ajr = jrin[ntime,:,:]; #pull it out
        if( FLG_smoothing_Jr == True ): #smooth if on
            for kk in range(49,0-1,-1): #start of MLAT loop
                jrin[ntime,:,kk] = uniform_filter1d_wrap(jrin[ntime,:,kk],3); #smooth in mlt at each latitude ring
            #END FOR kk
        #END IF
        
        #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        #Calc conducatances
        [spb, shb, efout[ntime,:,:]] = get_model_conductances(np.copy(jrin[ntime,:,:]),ayear,month,aday,time_hr[ntime],geo_lat[ntime,:,:],geo_long[ntime,:,:], \
                                                              FLG_smoothing_conductance=FLG_smoothing_conductance, \
                                                              FLG_smoothing_currentClip=FLG_smoothing_currentClip, \
                                                              FLG_smoothing_polarBoost=FLG_smoothing_polarBoost, \
                                                              FLG_smoothing_equatorBoost=FLG_smoothing_equatorBoost); #latest iteration
        
        #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        #Calc potential
        phib = potentializer_bot2top(jrin[ntime,:,:],spb,shb,thmax,thsn,delth,dtr,dph,re, \
                                     FLG_smoothing_phi_connectMidnight=FLG_smoothing_phi_connectMidnight, \
                                     FLG_smoothing_phi_lat=FLG_smoothing_phi_lat, \
                                     FLG_smoothing_phi_long=FLG_smoothing_phi_long);
        
        #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        #cycle through multipliers of phib to get the best agreement between jrin[ntime,:,:] and facout
        if( FLG_smooth_phiFix == True ):
            mplr = .1+.1*np.arange(31,dtype=np.float64)
            rms11 = np.zeros((31),dtype=np.float64)
            phibt = np.zeros((24,50),dtype=np.float64)     ##phib array to be tested 
            for jm in range(0,30+1):
                mxx = mplr[jm]
                phibt = mxx*phib
                [eeast, enorth, jeast, jnorth, jehall, jouleheat, jrcalc] = \
                    get_auroral_parameters(th50,spb,shb,phibt,FLG_smoothing_efield=FLG_smoothing_efield);
                #find the rmsd between jrin[ntime,:,:] and facout
                ampsum = np.sum(np.abs(jrcalc[:,15:45+1])*area[15:45+1]/1E12); #millions of amps or MA, from 55 to 85
                ampsum_ref = np.sum(np.abs(jrin[ntime,:,15:45+1])*area[15:45+1]/1E12); #millions of amps or MA, from 55 to 85
                rms11[jm] = np.abs(1-ampsum/ampsum_ref); 
            # END FOR     #end of jm loop
            min_subscript = np.where( rms11 == np.min(rms11) )[0][0]; #get an index I guess? guessing what min_subscript does in IDL
            mx0 = mplr[min_subscript]
            phib[:,0:45+1] = mx0*phib[:,0:45+1]; #adjust only below 85 (above is fiesty and has accurate values)
        #END IF
        
        #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        #Calc joule heating and current (for accuracy checking)
        if( FLG_magicks == True ):
            [eeast, enorth, jeast, jnorth, jehall, jouleheat, jrcalc] = \
                get_auroral_parameters(th50,spb,shb,phib,FLG_smoothing_efield=FLG_smoothing_efield);
            print(str(ntime+1)+'/'+str(time_cntStop)+' Error start: '+str(np.sum(np.abs(jrcalc - jrin[ntime,:,:]))));
            tic = time.time();
            
            weighter, _ = np.meshgrid(np.arange(1,51),np.arange(0,24)); #create some weights that favor high latitudes (numba don't like meshgrid)
            spb, shb, phib = optimizer(jrin[ntime,:,:], spb, shb, phib, th50, FLG_smoothing_efield=FLG_smoothing_efield, weighter=weighter);
            
            [eeast, enorth, jeast, jnorth, jehall, jouleheat, jrcalc] = \
                get_auroral_parameters(th50,spb,shb,phib,FLG_smoothing_efield=FLG_smoothing_efield);
            print(str(ntime+1)+'/'+str(time_cntStop)+' Error end: '+str(np.sum(np.abs(jrcalc - jrin[ntime,:,:])))+', took: '+str(time.time()-tic)+' sec');
        else:
            [eeast, enorth, jeast, jnorth, jehall, jouleheat, jrcalc] = \
                get_auroral_parameters(th50,spb,shb,phib,FLG_smoothing_efield=FLG_smoothing_efield);
        #END IF
           
        phimax = np.max(phib); #V, these had ,max_subscript
        phimin = np.min(phib); #V, these had ,min_subscript
        pocapo = phimax-phimin #V
        potdiff[ntime] = pocapo #V        
        
        spout[ntime,:,:] = spb[:,:];
        shout[ntime,:,:] = shb[:,:];
        phiout[ntime,:,:] = phib[:,:];
        eeout[ntime,:,:] = eeast[:,:];
        enout[ntime,:,:] = enorth[:,:];
        jrout[ntime,:,:] = jrcalc[:,:];
        jeout[ntime,:,:] = jeast[:,:];
        jnout[ntime,:,:] = jnorth[:,:];
        jehallout[ntime,:,:] = jehall[:,:];
        jouleout[ntime,:,:] = jouleheat[:,:];
        
    #     #************************************************************************
    #     # calculate the cartesian coordinates of the phib array elements
    #     for k in range(0,49+1):
    #         rr = 90.-mlat[k]
    #         for j in range(0,23+1):
    #             if j <= 23:
    #                 alpha = phmlt[j]*dtr
    #             #END IF
    #             if j == 23:
    #                 alpha = phmlt[0]*dtr
    #             #END IF
    #             rx[j,k] = 50.+rr*np.sin(alpha)
    #             ry[j,k] = 50.-rr*np.cos(alpha)
    #         # END FOR
    #     # END FOR
    #     # diff = np.concatenate(((facout[:,15:35+1]-facin[:,15:35+1])**2)).sum(); #sum it all in one go
    #     # rmst = (diff/1200.)**.5
    #     #    make a contour plot of the convection pattern
    #     # p = -100000+5000*np.arange(39,dtype=np.float64)
    # #    facplt = contour(phib,rx,ry,axis_style = 2,dimensions=[800,800],irregular = 1,     \
    # #        c_value = p,title = ptitl,c_label_show = 0,c_color = 0,color='black',     \
    # #        xmajor = 0,ymajor = 0,margin=[.15,.15,.15,.15])
        
        #************************************************************************
        jh[ntime] = np.sum(jouleout[ntime,:,1:49+1]*area[1:49+1]*1E4)/1E16; #GW, from 1 to 89, (uses area in cm2 b/c JH in erg/cm^2*s)
        powerg[ntime] = np.sum(efout[ntime,:,1:49+1]*area[1:49+1]*1E4)/1E16; #GW, from 1 to 89, (uses area in cm2 b/c JH in erg/cm^2*s)
        facamp[ntime] = np.sum(np.abs(jrout[ntime,:,1:49+1])*area[1:49+1]/1E12); #MA (mega-amps), from 1 to 89
        #************************************************************************
     
        #    calculate the magnetic perturbation at the 12 AE stations
        dtime = datetime.datetime(ayear, month, aday); #date time object for aacgmv2
        for isite in range(0,11+1):
            qlat = sitelat[isite];
            qlong = sitelong[isite];
            #convert to aacgm
            [out_lat, out_lon, _] = aacgmv2.convert_latlon(qlat, qlong, 120., dtime, method_code='G2A'); # see https://aacgmv2.readthedocs.io/en/latest/reference/aacgmv2.html#aacgmv2.wrapper.convert_latlon
            #    convert magnetic longitude to magnetic local time
            # yearin = ayear
            # intime = t0
            #cmlt = calc_mlt(yearin,intime,out_lon) #!!! This dtime needs more time support !!!
            cmlt = aacgmv2.convert_mlt(out_lon,dtime,m2a=False); # see https://aacgmv2.readthedocs.io/en/latest/reference/aacgmv2.html#aacgmv2.wrapper.convert_mlt
            sumdh = 0.
            for nfil in range(0,10+1):
                ifil = nfil-5
                flat = out_lat+np.float64(ifil)-40.
                dfil = 109*np.float64(ifil)
                rfil = np.sqrt(120*120.+dfil*dfil)
                interp_x = np.arange(0, jehall.shape[0], 1); #make indexes
                interp_y = np.arange(0, jehall.shape[1], 1); #make indexes
                # interp_xx, interp_yy = np.meshgrid(interp_x, interp_y); #make a meshgrid of them
                interp_fun = interpolate.interp2d(interp_y, interp_x, jehallout[ntime,:,:], kind='linear'); #IDL does bilinear with a 2D input to interpolate, accordion to comment on https://stackoverflow.com/a/8662355/2403531 interp2d does bilinear as well
                jfil0 = interp_fun(cmlt, flat); # above bits get this done: IDL->interpolate(jehall,cmlt,flat)
                jfil = jfil0/1000.        # current in Amps/meter
                dbfil = munaught*jfil*109./(2*np.pi*rfil)        #multiply by 109 to get current in amps
                dhfil = dbfil*120./rfil     # projection of the total db onto the h component
                sumdh = sumdh+dhfil
            # END FOR    # end of nfil loop
            dbhfac = 1.3
            dbh[isite,ntime] = dbhfac*sumdh*1.e9;
        # END FOR    #end of isite loop
        al_sim[ntime] = np.min(dbh[:,ntime]);
        au_sim[ntime] = np.max(dbh[:,ntime]);
        aesim[ntime] = au_sim[ntime]-al_sim[ntime];
    # END FOR    #end of ntime loop
    
    if( FLG_textOutputs == True ):
        dayin = str(start_yr[0])+str(start_mo[0])+str(start_dy[0]); #calc dayin string (used for text outputs)
        adelphi_textOutputs(dayin,hemi,time_cntStop,time_stepsAvail,time_hr,phmlt,th50,spout,shout, \
                                phiout,eeout,enout,jeout,jnout,jehallout,efout,jouleout, \
                                jrin,jrout,potdiff,powerg,jh,aesim);
    #END IF
    
    #set up output to check out with big code
    data_AMPERE['1D'] = {}; #prep dict
    data_AMPERE['1D']['JH integrated'] = jh; #GW
    data_AMPERE['1D']['EF integrated'] = powerg; #GW
    data_AMPERE['1D']['FAC Amp'] = facamp; #MA
    data_AMPERE['JH'] = np.flip(jouleout,axis=2).ravel(); #erg/(cm^2*s) [flip fixes directional issue since here lat is flipped so it goes lo lat -> hi lat but AMPERE comes in as hi lat -> lo lat]
    data_AMPERE['Phi'] = np.flip(phiout,axis=2).ravel();
    data_AMPERE['Jr_in'] = np.flip(jrin,axis=2).ravel(); #it's been slightly modified so sending it out too (smoothing)
    data_AMPERE['Jr_out'] = np.flip(jrout,axis=2).ravel();
    data_AMPERE['EF'] = np.flip(efout,axis=2).ravel(); #energy flux
    data_AMPERE['elec north'] = np.flip(enout,axis=2).ravel(); #energy flux
    data_AMPERE['elec east'] = np.flip(eeout,axis=2).ravel(); #energy flux
    data_AMPERE['1D']['AL sim'] = al_sim; #nT
    data_AMPERE['1D']['AU sim'] = au_sim; #nT
    data_AMPERE['1D']['AE sim'] = aesim; #nT
    
    cpcp = potdiff/1000.; # kV
    data_AMPERE['1D']['CPCP'] = potdiff/1000; # kV
    
    cpcp = potdiff/1000.; # kV
    data_AMPERE['1D']['CPCP'] = potdiff/1000; # kV
    
    data_AMPERE['1D']['FAC Amp * CPCP'] = cpcp*facamp; #kV*MA = GW
    
    data_AMPERE['1D']['time'] = np.unique(data_AMPERE['time']); #sec since start of current year, helps since larger array time is hard to return back to the 1D data
        
    return data_AMPERE
#END DEF

def adelphi_textOutputs(dayin,hemi,time_cntStop,time_stepsAvail,time_hr,phmlt,th50,spout,shout, \
                        phiout,eeout,enout,jeout,jnout,jehallout,efout,jouleout, \
                        jrin,jrout,potdiff,powerg,jh,aesim):
    #save a file here
    from Code.subfun_textNice import textNice
    #======================================================================
    print('file for 2D output in geomagnetic coordinates');
    pfilout='E:\Big Data\AMPERE\AuroraPHILE 2D Output '+dayin+' '+hemi+' Geomagnetic.txt';
    with open(pfilout, 'w') as fandle:
        for ip in range(0,time_cntStop):
            print(dayin+'\t'+format(time_hr[ip],'.5f'),' UT',file=fandle);
            for jp in range(0,23+1):
                print('MLT = '+format(phmlt[jp]/15.,'.5f')+'  Hours',file=fandle);
                print('MLAT  PED  HALL  PHI  EEAST  ENORTH  JEAST JNORTH  EFLUX   JHEAT   JRIN  JROUT',file=fandle);
                for kp in range(0,49+1):
                    print(format(th50[kp],'.0f')+'  '+format(spout[ip,jp,kp],'.1f')+'  '+format(shout[ip,jp,kp],'.1f')+'  '+format(phiout[ip,jp,kp]/1000.,'.1f')+ \
                        '  '+format(eeout[ip,jp,kp],'.1f')+'  '+format(enout[ip,jp,kp],'.1f')+'  '+format(jeout[ip,jp,kp],'.1f')+'  '+format(jehallout[ip,jp,kp],'.1f')+ \
                        '  '+format(efout[ip,jp,kp],'.1f')+'  '+format(jouleout[ip,jp,kp],'.1f')+ \
                        '  '+format(jrin[ip,jp,kp],'.2f')+'  '+format(jrout[ip,jp,kp],'.2f'),file=fandle);
                #END FOR jp
            #END FOR jp
        #END FOR ip
    #END WITH fandle
    #========================================================================
    print('file for 1D output');
    pfilout='E:\Big Data\AMPERE\AuroraPHILE 1D Output '+dayin+' '+hemi+'.txt'
    with open(pfilout, 'w') as fandle:
        print('  Date       UT     CPCP(kV) EnFlux(GW) JHeat(GW) Model_AE',file=fandle);
        for ip in range(0,time_cntStop):
            print(dayin+'\t'+str(np.round(time_hr[ip],3)).zfill(3)+'\t'+str(np.round(potdiff[ip]/1000.,1)).zfill(1)+'\t'+ \
                  str(np.round(powerg[ip],1)).zfill(1)+'\t'+str(np.round(jh[ip],1)).zfill(1)+'\t'+str(np.round(aesim[ip],1)).zfill(1),file=fandle);
        #END FOR ip
    #END WITH
    
    igeog=0
    if igeog > 0:
        from scipy.io.idl import readsav
        #======================================================================
        geogout = np.zeros((11,time_stepsAvail,24,50),dtype=np.float64);
        print('geog coords is going down');
        clong= np.arange(0,360,15);
        zfile='E:\Big Data\AMPERE\lat-long index_file'; #gotta read this file, it's saved IDL variables, lost in time I suspect [it's not]
        idl_save_data = readsav(zfile);
        zlatx = idl_save_data['zlatx'];
        zlonx = idl_save_data['zlonx'];
        y_colat = np.arange(0,50,1); #puts the data on a grid (who knows how the IDL one worked)
        x_long = np.arange(0,24,1);
        for ntime in range(0,time_cntStop):
            # uttime = (24.:ntime/time_stepsAvail);
            for jc in range(0,23+1):
                for kc in range(0,49+1):
                  latidx = zlatx[kc,jc,ntime];
                  lonidx = zlonx[kc,jc,ntime];
                  geogout[0,ntime,jc,kc] = bilinear(spout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[1,ntime,jc,kc] = bilinear(shout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[2,ntime,jc,kc] = bilinear(phiout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[3,ntime,jc,kc] = bilinear(eeout[ntime,:,:]/1000.,y_colat,x_long,latidx,lonidx);
                  geogout[4,ntime,jc,kc] = bilinear(enout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[5,ntime,jc,kc] = bilinear(jeout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[6,ntime,jc,kc] = bilinear(jnout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[7,ntime,jc,kc] = bilinear(efout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[8,ntime,jc,kc] = bilinear(jouleout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[9,ntime,jc,kc] = bilinear(jrin[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[10,ntime,jc,kc] = bilinear(jrout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                #END FOR kc
            #END FOR jc
            # print(str(ntime));
        #END FOR ntime #end of ntime loop
        #======================================================================
        print('file for 2D output for geographic coordinates');
        pfilout='E:\Big Data\AMPERE\AuroraPHILE 2D Output '+dayin+' '+hemi+' Geographic.txt'
        with open(pfilout, 'w') as fandle:
            for ip in range(0,time_cntStop):
                print(dayin+'\t'+format(time_hr[ip],'.5f')+' UT',file=fandle);
                for jp in range(0,23+1):
                    print('Longitude= ',format(clong[jp],'.5f'),'  Degrees East',file=fandle)
                    print('MLAT  PED  HALL  PHI  EEAST  ENORTH  JEAST JNORTH  EFLUX   JHEAT   JRIN  JROUT',file=fandle)
                    for kp in range(0,49+1):
                        print(format(th50[kp],'.0f')+'  '+format(geogout[0,ip,jp,kp],'.1f')+'  '+format(geogout[1,ip,jp,kp],'.1f')+'  '+format(geogout[2,ip,jp,kp]/1000.,'.1f')+ \
                        '  '+format(geogout[3,ip,jp,kp],'.1f')+'  '+format(geogout[4,ip,jp,kp],'.1f')+'  '+format(geogout[5,ip,jp,kp],'.1f')+ \
                        '  '+format(geogout[6,ip,jp,kp],'.1f')+'  '+format(geogout[7,ip,jp,kp],'.1f')+ \
                        '  '+format(geogout[8,ip,jp,kp],'.1f')+'  '+format(geogout[9,ip,jp,kp],'.2f')+'  '+format(geogout[10,ip,jp,kp],'.2f'),file=fandle)
                    #END FOR kp
                #END FOR jp
            #END FOR ip
        #END WITH
    #END IF igeog
#END DEF

def adelphi_textOutputs(dayin,hemi,time_cntStop,time_stepsAvail,time_hr,phmlt,th50,spout,shout, \
                        phiout,eeout,enout,jeout,jnout,jehallout,efout,jouleout, \
                        jrin,jrout,potdiff,powerg,jh,aesim):
    #save a file here
    # from Code.subfun_textNice import textNice
    #======================================================================
    print('file for 2D output in geomagnetic coordinates');
    pfilout='E:\Big Data\AMPERE\AuroraPHILE 2D Output '+dayin+' '+hemi+' Geomagnetic.txt';
    with open(pfilout, 'w') as fandle:
        for ip in range(0,time_cntStop):
            print(dayin+'\t'+format(time_hr[ip],'.5f'),' UT',file=fandle);
            for jp in range(0,23+1):
                print('MLT = '+format(phmlt[jp]/15.,'.5f')+'  Hours',file=fandle);
                print('MLAT  PED  HALL  PHI  EEAST  ENORTH  JEAST JNORTH  EFLUX   JHEAT   JRIN  JROUT',file=fandle);
                for kp in range(0,49+1):
                    print(format(th50[kp],'.0f')+'  '+format(spout[ip,jp,kp],'.1f')+'  '+format(shout[ip,jp,kp],'.1f')+'  '+format(phiout[ip,jp,kp]/1000.,'.1f')+ \
                        '  '+format(eeout[ip,jp,kp],'.1f')+'  '+format(enout[ip,jp,kp],'.1f')+'  '+format(jeout[ip,jp,kp],'.1f')+'  '+format(jehallout[ip,jp,kp],'.1f')+ \
                        '  '+format(efout[ip,jp,kp],'.1f')+'  '+format(jouleout[ip,jp,kp],'.1f')+ \
                        '  '+format(jrin[ip,jp,kp],'.2f')+'  '+format(jrout[ip,jp,kp],'.2f'),file=fandle);
                #END FOR jp
            #END FOR jp
        #END FOR ip
    #END WITH fandle
    #========================================================================
    print('file for 1D output');
    pfilout='E:\Big Data\AMPERE\AuroraPHILE 1D Output '+dayin+' '+hemi+'.txt'
    with open(pfilout, 'w') as fandle:
        print('  Date       UT     CPCP(kV) EnFlux(GW) JHeat(GW) Model_AE',file=fandle);
        for ip in range(0,time_cntStop):
            print(dayin+'\t'+str(np.round(time_hr[ip],3)).zfill(3)+'\t'+str(np.round(potdiff[ip]/1000.,1)).zfill(1)+'\t'+ \
                  str(np.round(powerg[ip],1)).zfill(1)+'\t'+str(np.round(jh[ip],1)).zfill(1)+'\t'+str(np.round(aesim[ip],1)).zfill(1),file=fandle);
        #END FOR ip
    #END WITH
    
    igeog=0
    if igeog > 0:
        from scipy.io.idl import readsav
        #======================================================================
        geogout = np.zeros((11,time_stepsAvail,24,50),dtype=np.float64);
        print('geog coords is going down');
        clong= np.arange(0,360,15);
        zfile='E:\Big Data\AMPERE\lat-long index_file'; #gotta read this file, it's saved IDL variables, lost in time I suspect [it's not]
        idl_save_data = readsav(zfile);
        zlatx = idl_save_data['zlatx'];
        zlonx = idl_save_data['zlonx'];
        y_colat = np.arange(0,50,1); #puts the data on a grid (who knows how the IDL one worked)
        x_long = np.arange(0,24,1);
        for ntime in range(0,time_cntStop):
            # uttime = (24.:ntime/time_stepsAvail);
            for jc in range(0,23+1):
                for kc in range(0,49+1):
                  latidx = zlatx[kc,jc,ntime];
                  lonidx = zlonx[kc,jc,ntime];
                  geogout[0,ntime,jc,kc] = bilinear(spout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[1,ntime,jc,kc] = bilinear(shout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[2,ntime,jc,kc] = bilinear(phiout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[3,ntime,jc,kc] = bilinear(eeout[ntime,:,:]/1000.,y_colat,x_long,latidx,lonidx);
                  geogout[4,ntime,jc,kc] = bilinear(enout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[5,ntime,jc,kc] = bilinear(jeout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[6,ntime,jc,kc] = bilinear(jnout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[7,ntime,jc,kc] = bilinear(efout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[8,ntime,jc,kc] = bilinear(jouleout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[9,ntime,jc,kc] = bilinear(jrin[ntime,:,:],y_colat,x_long,latidx,lonidx);
                  geogout[10,ntime,jc,kc] = bilinear(jrout[ntime,:,:],y_colat,x_long,latidx,lonidx);
                #END FOR kc
            #END FOR jc
            # print(str(ntime));
        #END FOR ntime #end of ntime loop
        #======================================================================
        print('file for 2D output for geographic coordinates');
        pfilout='E:\Big Data\AMPERE\AuroraPHILE 2D Output '+dayin+' '+hemi+' Geographic.txt'
        with open(pfilout, 'w') as fandle:
            for ip in range(0,time_cntStop):
                print(dayin+'\t'+format(time_hr[ip],'.5f')+' UT',file=fandle);
                for jp in range(0,23+1):
                    print('Longitude= ',format(clong[jp],'.5f'),'  Degrees East',file=fandle)
                    print('MLAT  PED  HALL  PHI  EEAST  ENORTH  JEAST JNORTH  EFLUX   JHEAT   JRIN  JROUT',file=fandle)
                    for kp in range(0,49+1):
                        print(format(th50[kp],'.0f')+'  '+format(geogout[0,ip,jp,kp],'.1f')+'  '+format(geogout[1,ip,jp,kp],'.1f')+'  '+format(geogout[2,ip,jp,kp]/1000.,'.1f')+ \
                        '  '+format(geogout[3,ip,jp,kp],'.1f')+'  '+format(geogout[4,ip,jp,kp],'.1f')+'  '+format(geogout[5,ip,jp,kp],'.1f')+ \
                        '  '+format(geogout[6,ip,jp,kp],'.1f')+'  '+format(geogout[7,ip,jp,kp],'.1f')+ \
                        '  '+format(geogout[8,ip,jp,kp],'.1f')+'  '+format(geogout[9,ip,jp,kp],'.2f')+'  '+format(geogout[10,ip,jp,kp],'.2f'),file=fandle)
                    #END FOR kp
                #END FOR jp
            #END FOR ip
        #END WITH
    #END IF igeog
#END DEF


def bilinear(data, x_pts, y_pts, xInterp_pts, yInterp_pts):
    interper = interpolate.interp2d(x_pts, y_pts, data, kind='linear');
    return interper(xInterp_pts, yInterp_pts);
#END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def uniform_filter1d_wrap(inputz,window_size): #numba-friendly version
    accumulator = np.zeros_like(inputz); #prep
    for jk in range(0+int(1-np.ceil(window_size/2)),window_size+int(1-np.ceil(window_size/2))): #semi-clever numbering lets roll do the work
        accumulator += np.roll(inputz,jk); #acculumate values
    #END FOR jk
    return accumulator/window_size #divide, note to do this super accurate (round-off-safer) record the arrays and then call mean here
#END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def uniform_filter1d_orig(inputz,window_size): #numba-friendly version that matches IDL's default smooth mode (end pts are the same?? the sum is broken?? literally not how filters should work??)
    accumulator = np.zeros_like(inputz); #prep
    for jk in range(0+int(1-np.ceil(window_size/2)),window_size+int(1-np.ceil(window_size/2))): #semi-clever numbering lets roll do the work
        accumulator += np.roll(inputz,jk); #acculumate values
    #END FOR jk
    accumulator = accumulator/window_size; #apply the window
    accumulator[0] = inputz[0]; #end points are unchanged, weird but matching it
    accumulator[-1] = inputz[-1]; #end points are unchanged, weird but matching it
    return accumulator #divide, note to do this super accurate (round-off-safer) record the arrays and then call mean here
#END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def uniform_filter1d_nearest(inputz,window_size): #numba-friendly version
    accumulator = np.zeros_like(inputz);
    for jk in range(0+int(1-np.ceil(window_size/2)),window_size+int(1-np.ceil(window_size/2))):
        midstep = np.roll(inputz,jk); #roll, this code will likely not hold up if window_size > inputz.size
        if( jk < 0 ):
            midstep[jk:] = inputz[-1]; #nearest needs the end values where it rolled over
        else:
            midstep[0:jk] = inputz[0]; #nearest needs the end values where it rolled over
        #END IF
        accumulator += midstep; #acculumate values
    #END FOR jk
    return accumulator/window_size #divide, note to do this super accurate (round-off-safer) record the arrays and then call mean here
#END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def get_model_conductances(ajr,ayear,month,aday,authr,geo_lat,geo_long,FLG_smoothing_conductance=True,FLG_smoothing_currentClip=True,FLG_smoothing_polarBoost=False,FLG_smoothing_equatorBoost=False):
    
    # timeTemp = authr; #get hr variable
    # timeTemp_hr = np.int32(timeTemp); #hr, get the hours
    # timeTemp_min = (timeTemp - timeTemp_hr)*60; #min, get the minutes (still float)
    # if( np.isclose( timeTemp_min, np.round(timeTemp_min) ) ):
    #     timeTemp_min = np.round(timeTemp_min); #round it, deals with the fractional float
    # #END IF
    # timeTemp_min = np.int32(timeTemp_min); #min, get the minutes
    # timeTemp_sec = np.int32( ((timeTemp - timeTemp_hr)*60 - timeTemp_min)*60 ); #sec, get the seconds
    # timeTemp_usec = np.int32(np.round( (((timeTemp - timeTemp_hr)*60 - timeTemp_min)*60 - timeTemp_sec)*1E6 )); #ends with microseconds, so round
    # dtime = datetime.datetime(ayear, month, aday, hour=timeTemp_hr, minute=timeTemp_min, second=timeTemp_sec, microsecond=timeTemp_usec); #date time object for aacgmv2
    
    #================tuning parameters================
    sigmap = 2.
    sigmah = 2.
    facmin = .1
    splatsm = 3    #for latitudinal smoothing
    shlatsm = 3    #for latitudinal smoothing
    #=================================================
    dtr = np.pi/180.
    if( ((np.mod(ayear,4) == 0) & (np.mod(ayear,100) != 0)) | (np.mod(ayear,400) == 0) ):
        firstday = [1,32,61,92,122,153,183,214,245,275,306,336]; #first day for leap years
    else:
        firstday = [1,32,60,91,121,152,182,213,244,274,305,335];    #for calculating julian day
    #END IF
    # phmlt = 15*np.arange(24,dtype=np.float64)
    # phmlt_tru = np.arange(24,dtype=np.float64)
    sp = np.zeros((24,50),dtype=np.float64)
    sh = np.zeros((24,50),dtype=np.float64)
    spb = np.zeros((24,50),dtype=np.float64)
    shb = np.zeros((24,50),dtype=np.float64)
    ef = np.zeros((24,50),dtype=np.float64)
    #    conductance arrays from newer end2end_v2 analysis with sk sigmas
    spm0 = np.zeros((25),dtype=np.float64)
    spm1 = np.zeros((25),dtype=np.float64)
    spp0 = np.zeros((25),dtype=np.float64)
    spp1 = np.zeros((25),dtype=np.float64)
    shm0 = np.zeros((25),dtype=np.float64)
    shm1 = np.zeros((25),dtype=np.float64)
    shp0 = np.zeros((25),dtype=np.float64)
    shp1 = np.zeros((25),dtype=np.float64)
    fluxg = 0.
    #                spm0   spm1    spp0    spp1    shm0   shm1    shp0    shp1
    aa = np.asarray([5.00,  4.20,   7.70,   8.70,   -3.20, 6.80,   -7.30,  14.80])
    bb = np.asarray([-0.80, 1.10,   -1.80,  4.60,   -3.60, -1.50,  5.60,   -10.40])
    cc = np.asarray([60.90, 318.60, 139.00, 327.10, 21.90, 184.90, 100.90, 129.90])
    ang = 15*np.arange(25,dtype=np.float64)
    for iw in range(0,24+1):
        spm0[iw] = aa[0]+bb[0]*np.cos((cc[0]+ang[iw])*dtr)
        spm1[iw] = aa[4]+bb[4]*np.cos((cc[4]+ang[iw])*dtr)
        spp0[iw] = aa[1]+bb[1]*np.cos((cc[1]+ang[iw])*dtr)
        spp1[iw] = aa[5]+bb[5]*np.cos((cc[5]+ang[iw])*dtr)
        shm0[iw] = aa[2]+bb[2]*np.cos((cc[2]+ang[iw])*dtr)
        shm1[iw] = aa[6]+bb[6]*np.cos((cc[6]+ang[iw])*dtr)
        shp0[iw] = aa[3]+bb[3]*np.cos((cc[3]+ang[iw])*dtr)
        shp1[iw] = aa[7]+bb[7]*np.cos((cc[7]+ang[iw])*dtr)
    # END FOR
    # SMOOTHED energy flux fits from new GUVI analysis--WITH THE FOLLOWING CORRECTIONS MADE LATER IN CODE
    # FOR AGREEMENT WITH NUMBERS IN THE PUBLISHED TABLE:
    efp1g = np.asarray([11.60, 11.00,    10.06,    8.25, 6.50, 5.75, 5.00, 3.50, 1.89, 1.68, 1.48, 1.10, 0.82,  \
        0.76, 0.70, 1.70, 2.70, 3.35, 3.99, 5.32, 6.64, 8.80, 11.00,    11.50, 11.6]);
    efp0g = np.asarray([1.70, 1.60, 1.60, 1.65, 1.70, 1.66, 1.63, 1.68, 1.70, 1.67, 1.65, 1.67, 1.75,  \
        1.90, 2.15, 2.20, 2.29, 2.33, 2.38, 2.34, 2.30, 2.15, 2.00, 1.85, 1.7]);
    efm1g = np.asarray([-5.00,    -4.75,    -4.50,    -3.25,    -2.00,    -1.50,    -1.00,    -0.70,    -0.41,    -0.45,    -0.50,  \
        -0.51,    -0.51,    -0.55,    -0.60,    -1.10,    -1.60,    -3.30,    -5.00,    -5.25,    -5.50,    -5.40,    -5.30,    -5.00, -5.]);
    efm0g = np.asarray([3.47, 3.39, 3.32, 3.19, 3.07, 2.88, 2.69, 2.37, 2.06, 1.97, 1.89, 1.91,  \
        1.93, 1.98, 2.14, 2.26, 2.39, 2.85, 3.33, 3.51, 3.70, 3.70, 3.70, 3.58, 3.47]);
    efm0g = efm0g-1.2
    efp0g = efp0g-0.5
    
    if( FLG_smoothing_polarBoost == True ):
        ajr = np.append(np.copy(ajr),np.expand_dims(np.repeat(np.mean(ajr[:,-1]),ajr.shape[0]),axis=1),axis=1); #expand to include 90 via mean of 89 (so that conductance for 89 can be calculated)
        # ajr_arr = np.empty(ajr.shape[0],dtype=np.float64);
        # for j in range(0,23+1):
        #     oppin = np.mod(j+12,24);
        #     ajr_arr[j] = (ajr[oppin,-1]-ajr[j,-1])/2; #90 value is a mean of value below and across
        # #END FOR j
        # ajr = np.append(np.copy(ajr),np.expand_dims(ajr_arr,axis=1),axis=1); #expand to include 90 (so that conductance for 89 can be calculated)
        kk_max = 49;
    else:
        kk_max = 48;
    #END IF
    if( FLG_smoothing_equatorBoost == True ):
        # ajr = np.fliplr(np.append(np.fliplr(np.copy(ajr)),np.expand_dims(np.repeat(np.mean(ajr[:,0]/2),ajr.shape[0]),axis=1),axis=1)); #expand to include 39 via mean of 40 and 0 (40/2 essentially) (so that conductance for 40 can be calculated)
        ajr = np.fliplr(np.append(np.fliplr(np.copy(ajr)),np.expand_dims(np.repeat(0,ajr.shape[0]),axis=1),axis=1)); #expand to include 39 via 0 FAC assumption (so that conductance for 40 can be calculated) [0 is bad for conductivity calcs, causes inf E field]
        kk_min = 0;
    else:
        kk_min = 1;
    #END IF
    
    # calculate conductances using Model07b
    for jj in range(0,23+1):
        for kk in range(kk_min,kk_max+1): #expand to include 49
        # for kk in range(kk_max,kk_min-1,-1): #no diff
            if( FLG_smoothing_currentClip == True ):
                if (np.abs(ajr[jj,kk-1]) < 3) & (np.abs(ajr[jj,kk]) < 3) & (np.abs(ajr[jj,kk+1]) < 3):
                    facm1 = ajr[jj,kk-1]
                    facp1 = ajr[jj,kk+1]
                    fac0 = ajr[jj,kk]
                else:
                    ajr[jj,kk-1] = 0.
                    ajr[jj,kk+1] = 0.
                    ajr[jj,kk] = 0.
                    facm1 = 0.
                    facp1 = 0.
                    fac0 = 0.
                #END IF
            else:
                facm1 = ajr[jj,kk-1]
                facp1 = ajr[jj,kk+1]
                fac0 = ajr[jj,kk]
            #END IF
            # if fac is changing sign, use the average values on either side
            pedm1 = 0.
            pedp1 = 0.
            if fac0*facm1 < 0:
                if facm1 < 0:
                    pedm1 = spm0[jj]+spm1[jj]*facm1
                    halm1 = shm0[jj]+shm1[jj]*facm1
                    fluxm1g = efm0g[jj]+efm1g[jj]*facm1
                # END IF
                if facm1 >= 0:
                    pedm1 = spp0[jj]+spp1[jj]*facm1
                    halm1 = shp0[jj]+shp0[jj]*facm1
                    fluxm1g = efp0g[jj]+efp1g[jj]*facm1
                # END IF
                if facp1 < 0:
                    pedp1 = spm0[jj]+spm1[jj]*facp1
                    halp1 = shm0[jj]+shm1[jj]*facp1
                    fluxp1g = efm0g[jj]+efm1g[jj]*facp1
                # END IF
                if facp1 >= 0:
                    pedp1 = spp0[jj]+spp1[jj]*facp1
                    halp1 = shp0[jj]+shp1[jj]*facp1
                    fluxp1g = efp0g[jj]+efp1g[jj]*facp1
                # END IF
                sigmap = .5*(pedm1+pedp1)
                sigmah = .5*(halm1+halp1)
                fluxg = .5*(fluxm1g+fluxp1g)
                # hoverp = sigmah/sigmap
            else: #if facs are not changing sign
                if fac0 < -facmin:
                    sigmap = spm0[jj]+spm1[jj]*fac0
                    sigmah = shm0[jj]+shm1[jj]*fac0
                    fluxg = efm0g[jj]+efm1g[jj]*fac0
                elif fac0 > facmin:
                    sigmap = spp0[jj]+spp1[jj]*fac0
                    sigmah = shp0[jj]+shp1[jj]*fac0
                    fluxg = efp0g[jj]+efp1g[jj]*fac0
                else:
                    sigmap_a = spp0[jj]+spp1[jj]*fac0
                    sigmah_a = shp0[jj]+shp1[jj]*fac0
                    fluxg_a = efp0g[jj]+efp1g[jj]*fac0
                    sigmap_b = spm0[jj]+spm1[jj]*fac0
                    sigmah_b = shm0[jj]+shm1[jj]*fac0
                    fluxg_b = efm0g[jj]+efm1g[jj]*fac0
                    
                    sigmap = 0.5*(sigmap_a+sigmap_b)
                    sigmah = 0.5*(sigmah_a+sigmah_b)
                    fluxg = 0.5*(fluxg_a+fluxg_b)
                # END IF
            #END IF
            if( FLG_smoothing_currentClip == True ):
                avg3=(facm1+fac0+facp1)/3.
                if (np.abs(np.abs(avg3)-facmin) < 1E-8) | (np.abs(avg3) < facmin):
                    sigmap = 2.
                    sigmah = 2.
                    fluxg = 0.
                # END IF
            #END IF
            ef[jj,kk] = fluxg
            #add in the conductance from solar illumination
            julday = firstday[month-1]+aday-1 #needed -1 to be correct
            # #calculate the geographic coordinates from magnetic coordinates
            # # maglon = phmlt[jj]-15*authr+75.
            # maglon = aacgmv2.convert_mlt(phmlt_tru[jj], dtime, m2a=True);
            # [szlat, szlong, _] = aacgmv2.convert_latlon(th50[kk], maglon, 150., dtime, method_code='A2G'); # see https://aacgmv2.readthedocs.io/en/latest/reference/aacgmv2.html#aacgmv2.wrapper.convert_latlon
            # if szlong < 0:
            #     szlong = 360.+szlong
            # #END IF
            
            # [ thetasun, azimuth] = sunangle(julday, authr, geo_lat[kk,jj], geo_long[kk,jj]) #missing this function
            # [ thetasun, azimuth] = sunangle(julday, authr, geo_lat[kk,jj], geo_long[kk,jj]); #actually did get it, but less accurate
            thetasun = sunAlsoRises_zenith_numba(geo_lat[kk,jj], geo_long[kk,jj], ayear, month, aday, authr, 0, 0); #but I have my own function so using it
            # earthLoc = coord.EarthLocation(lat=geo_lat[kk,jj]*u.deg, lon=szlong*u.deg, height=150.*u.m); #get location
            # sunAltAz = sun.transform_to(coord.AltAz(obstime=timeTime,location=earthLoc)); #get altaz
            # # # azimuth = sunAltAz.az.value;
            # thetasun = 90-sunAltAz.alt.value; #accurate but not using it, too slow
            e1 = ayear+julday/365.-2010.
            f1 = np.pi*e1/8.
            f107 = 80.+70*np.sin(f1)*np.sin(f1)
            #f107 = 70.
            #print('F10.7 set to : ',f107
            if thetasun < 89:
                chi = thetasun*np.pi/180.
                termx=(f107*np.cos(chi))**.5
                spsun=.88*termx
                shsun = 1.5*termx
            else:
                spsun = 0.
                shsun = 0.
            #END IF
            sp[jj,kk] = np.sqrt(sigmap**2+spsun**2)
            sh[jj,kk] = np.sqrt(sigmah**2+shsun**2)
        # END FOR # end of mlat loop
    # END FOR # end of mlt loop
    if( FLG_smoothing_conductance == 1 ):
        # interpolate the arrays to the specified latitude spacing
        for jy in range(0,23+1):
            sp[jy,:] = uniform_filter1d_orig(sp[jy,:],splatsm)
            sh[jy,:]= uniform_filter1d_orig(sh[jy,:],shlatsm)
        # END FOR
    elif( FLG_smoothing_conductance == 2 ):
        for kx in range(0,49+1):
            sp[:,kx] = uniform_filter1d_wrap(sp[:,kx],splatsm)
            sh[:,kx]= uniform_filter1d_wrap(sh[:,kx],shlatsm)
        # END FOR
    elif( FLG_smoothing_conductance == 3 ):
        for kx in range(0,49+1):
            sp[:,kx] = uniform_filter1d_wrap(sp[:,kx],splatsm)
            sh[:,kx]= uniform_filter1d_wrap(sh[:,kx],shlatsm)
        # END FOR
        for jy in range(0,23+1):
            sp[jy,:] = uniform_filter1d_orig(sp[jy,:],splatsm)
            sh[jy,:]= uniform_filter1d_orig(sh[jy,:],shlatsm)
        # END FOR
    #END IF
    spb = sp
    shb = sh
    
    return spb, shb, ef
#END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def potentializer_top2bot(ajr,spb,shb,thmax,thsn,delth,dtr,dph,re,FLG_smoothing_phi_connectMidnight=True,FLG_smoothing_phi_lat=True,FLG_smoothing_phi_long=True):
    phib = np.zeros((24,50),dtype=np.float64);
    #================tuning parameters================SAVE THIS SET
    # facmin=.1
    # splonsm = 3    #for longitudinal smoothing
    # shlonsm = 3    #for longitudinal smoothing
    latlast = 88 #latitude above which interpolation is used# default is 74
    maxlat = 48    #maximum latitude for integration    default(best)=48
    smvallt = 3; #default is 3
    # #================tuning parameters================THIS SET WORKS EVEN BETTER
    # # facmin=.1
    # # splonsm = 5    #for longitudinal smoothing
    # # shlonsm = 5    #for longitudinal smoothing
    # latlast = 70. #latitude above which interpolation is used# default is 75
    # maxlat = 48    #maximum latitude for integration    default(best)=45
    # smvallt = 3
    #=================================================
    # dtr = np.pi/180.
    # thmin = 0.
    # thmax = 50.
    # delth = 1.
    # # delthkm = delth*109.
    # #    arrays that contain the angles in spherical coordinates
    # thsn = delth*np.arange(50,dtype=np.float64)
    # dph = 15*dtr
    # phmlt = 15*np.arange(24,dtype=np.float64)
    # mlat=(90.-thmax)+delth*np.arange(50,dtype=np.float64)    #the magnetic latitudes corresponding to the theta array
    # # th50 = 40.+np.arange(50,dtype=np.float64)
    # re = 6300.e3        # Earth radius in m
    # smphi = np.zeros((24),dtype=np.float64)
    # smth = np.zeros((50),dtype=np.float64)
    # t = np.zeros((9),dtype=np.float64)
    # sumjr = 0. # for integrating facs over the hemisphere
    #start electrostatic potential calculation
    # boundary conditions
    # phib[:,0] = 0. #potential along 40-degree latitude circle is assumed = 0
    # phib[:,1] = 0.
    dthp = -delth*dtr #because the integration goes toward decreasing theta
    
    # print(str(ntime)+' '+format(np.sum(spb),'.5f')+' '+format(np.mean(spb),'.5f')+' '+format(np.sum(shb),'.5f')+' '+format(np.mean(shb),'.5f'))
    lastk = latlast-40
    # for k in range(1,maxlat+1): # start calculation at 41 degrees
    for k in range(maxlat,1-1,-1): # start calculation at 88 degrees
        theta = (thmax-thsn[k]+.5)*dtr #0.5 added to center on half-degree mark
        sinth = (np.sin(theta))
        cotth = np.cos(theta)/sinth
        kprev = k+1;
        # kprev = -1; #essentially sets phi prev to 0
        knext = k-1;
        for j in range(0,23+1):
            # if( k == 1 ):
            #     print('wut')
            # #END IF
            if k <= lastk:
                if j == 0:
                    phibjm1 = phib[23,k]
                    spbjm1 = spb[23,k]
                    shbjm1 = shb[23,k]
                else:
                    phibjm1 = phib[j-1,k]
                    spbjm1 = spb[j-1,k]
                    shbjm1 = shb[j-1,k]
                #END IF
                if j == 23:
                    phibjp1 = phib[0,k]
                    spbjp1 = spb[0,k]
                    shbjp1 = shb[0,k]
                else:
                    phibjp1 = phib[j+1,k]
                    spbjp1 = spb[j+1,k]
                    shbjp1 = shb[j+1,k]
                #END IF
                # print(str(j)+' '+str(k)+' '+format(phibjm1,'.5f')+\
                #       ' '+format(spbjm1,'.5f')+' '+format(shbjm1,'.5f')+\
                #       ' '+format(phibjp1,'.5f')+' '+format(spbjp1,'.5f')+\
                #       ' '+format(shbjp1,'.5f'))
                
                sigmap = spb[j,k]
                # sigmah = shb[j,k]
                dspdph=(spbjp1-spbjm1)/2./dph
                dshdph=(shbjp1-shbjm1)/2./dph
                
                dspdth=((spb[j,k+1]-spb[j,k-1])/2./dthp)
                dshdth=((shb[j,k+1]-shb[j,k-1])/2./dthp)
                dg1 = cotth*sigmap/(2*dthp)
                ng1 = dg1*phib[j,kprev] #k-1==kprev in orig
                dl1 = sigmap/(dthp*dthp)
                nl1 = 2*sigmap*phib[j,k]/(dthp*dthp)
                nl2 = sigmap*phib[j,kprev]/(dthp*dthp)
                dl2 = dspdth/(2*dthp)
                nl3 = dl2*phib[j,kprev]
                nl4 = dshdth*phibjp1/(sinth*2*dph)
                nl5 = dshdth*phibjm1/(sinth*2*dph)
                ds1 = dshdph/(2*dthp)/sinth
                ns1 = sigmap*phibjp1/(sinth*sinth*dph*dph)
                ns2 = 2*sigmap*phib[j,k]/(sinth*sinth*dph*dph)
                ns3 = sigmap*phibjm1/(sinth*sinth*dph*dph)
                ns4 = dspdph*phibjp1/(sinth*sinth*dph*2*dthp)
                ns5 = dspdph*phibjm1/(sinth*sinth*dph*2*dthp)
                ns6 = dshdph*phib[j,kprev]/(2*dthp)
                
                num = -re*re*ajr[j,k]*1.e-6+ng1+nl1-nl2+nl3-nl4+nl5-ns1+ns2-ns3-ns4+ns5-ns6
                denom = dg1+dl1+dl2-ds1
                
                # #END IF
                if denom != 0:
                    phib[j,knext] = num/denom
                #END IF
                if denom == 0:
                    phib[j,knext] = phib[j,k]
                #END IF
            # END IF
            
    #         if( FLG_smoothing_phi_connectMidnight == True ):
    #             phib[23,:] = phib[22,:]+.33*(phib[1,:]-phib[22,:])
    #             phib[0,:] = phib[22,:]+.6667*(phib[1,:]-phib[22,:])
    #         #END IF
    #         if( FLG_smoothing_phi_long == True ):
    #             # smooth in mlat and mlt
    #             if( smvallt < 3 ):
    #                 smphi[0:23+1] = uniform_filter1d_wrap(phib[0:23+1,k+1],smvallt+1);
    #             else:
    #                 smphi[0:23+1] = uniform_filter1d_wrap(phib[0:23+1,k+1],smvallt);
    #             #END IF
    #             phib[0:23+1,k+1] = smphi[0:23+1]
    #         #END IF
    #         if( FLG_smoothing_phi_lat == True ):
    #             for ism in range(0,23+1):
    #                 if k > 3:
    #                     smth[0:k+1+1] = uniform_filter1d_orig(phib[ism,0:k+1+1],2+1);     
    #                     phib[ism,0:k+1] = smth[0:k+1]
    #                 # END IF
    #             # END FOR
    #         #END IF
    #         # print(str(j)+' '+str(k)+' '+format(phib[j,k-1],'.5f')+\
    #         #       ' '+format(phib[j,k],'.5f')+' '+format(phib[j,k+1],'.5f'))
    #     #END FOR j
    # #END FOR k
    # for k in range(1,maxlat+1+1): # start calculation at 41 degrees
    # # for k in range(maxlat,1-1,-1): # start calculation at 88 degrees
    #     theta = (thmax-thsn[k]+.5)*dtr #0.5 added to center on half-degree mark
    #     sinth = (np.sin(theta))
    #     cotth = np.cos(theta)/sinth
    #     for j in range(0,23+1):
            # ****************************************************
            # interpolate the def above latitude latlast
            oppin = j+12
            if oppin >= 24:
                oppin = oppin-24
            #END IF
            kp1 = k+1
            if k > lastk:
                km30 = k-lastk-1
                mkp50 = 2*(50-lastk+1)
                # phib[j,kp1] = phib[j,29]+(km30)*(phib[oppin,29]-phib[j,29])/(2.*(mkp50)) #OG
                phib[j,kp1] = phib[j,lastk-1]+(km30)*(phib[oppin,lastk-1]-phib[j,lastk-1])/(2.*(mkp50)) #adjustment
                # phib[j,k] = phib[j,k-1]+(km30)*(phib[oppin,k-1]-phib[j,k-1])/(2.*(mkp50)) #adjustment
            # END IF
        #END IF
        # END FOR # end of mlt loop
        
        if( FLG_smoothing_phi_connectMidnight == True ):
            # connect the points across midnight
            phib[23,:] = phib[22,:]+.33*(phib[1,:]-phib[22,:])
            phib[0,:] = phib[22,:]+.6667*(phib[1,:]-phib[22,:])
        #END IF
        if( FLG_smoothing_phi_long == True ):
            # smooth in mlat and mlt
            if( smvallt < 3 ):
                smphi = uniform_filter1d_wrap(phib[0:23+1,k+1],smvallt+1);
            else:
                smphi = uniform_filter1d_wrap(phib[0:23+1,k+1],smvallt);
            #END IF
            phib[:,k+1] = smphi;
        #END IF
        if( FLG_smoothing_phi_lat == True ):
            for ism in range(0,23+1):
                # if k > 3:
                smth = uniform_filter1d_orig(phib[ism,0:k+1+1],2+1); #essential for integration control
                phib[ism,0:k+1] = smth[0:k+1];
                # smth = uniform_filter1d_nearest(phib[ism,k-1:k+1+1],2+1); #not enough  
                # phib[ism,k-1:k+1+1] = smth;
                # END IF
            # END FOR
        #END IF
    # END FOR # end of mag lat loop
    return phib
#END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def potentializer_bot2top(ajr,spb,shb,thmax,thsn,delth,dtr,dph,re,FLG_smoothing_phi_connectMidnight=True,FLG_smoothing_phi_lat=True,FLG_smoothing_phi_long=True):
    phib = np.zeros((24,50),dtype=np.float64);
    #================tuning parameters================SAVE THIS SET
    # latlast = 88 #latitude above which interpolation is used# default is 74
    # maxlat = 48    #maximum latitude for integration    default(best)=48
    smvallt = 3; #default is 3
    #=================================================
    dthp = -delth*dtr #because the integration goes toward decreasing theta
    
    if( smvallt == 2 ):
        smvallt = 3; #IDL's code only operates on 3 or more so make 2->3, 1 or 0 won't do anything so leave it to nothing on that I guess
    #END IF
    # lastk = latlast-40
    for k in range(1,48+1,1): # start calculation at 41 degrees
        theta = (thmax-thsn[k]+.5)*dtr #0.5 added to center on half-degree mark
        sinth = (np.sin(theta))
        cotth = np.cos(theta)/sinth
        kprev = k+1;
        knext = k-1;
        for j in range(0,23+1):
            if j == 0:
                phibjm1 = phib[23,k]
                spbjm1 = spb[23,k]
                shbjm1 = shb[23,k]
            else:
                phibjm1 = phib[j-1,k]
                spbjm1 = spb[j-1,k]
                shbjm1 = shb[j-1,k]
            #END IF
            if j == 23:
                phibjp1 = phib[0,k]
                spbjp1 = spb[0,k]
                shbjp1 = shb[0,k]
            else:
                phibjp1 = phib[j+1,k]
                spbjp1 = spb[j+1,k]
                shbjp1 = shb[j+1,k]
            #END IF
            
            sigmap = spb[j,k]
            # sigmah = shb[j,k]
            dspdph=(spbjp1-spbjm1)/2./dph
            dshdph=(shbjp1-shbjm1)/2./dph
            
            dspdth=((spb[j,k+1]-spb[j,k-1])/2./dthp)
            dshdth=((shb[j,k+1]-shb[j,k-1])/2./dthp)
            dg1 = cotth*sigmap/(2*dthp)
            ng1 = dg1*phib[j,kprev] #k-1==kprev in orig
            dl1 = sigmap/(dthp*dthp)
            nl1 = 2*sigmap*phib[j,k]/(dthp*dthp)
            nl2 = sigmap*phib[j,kprev]/(dthp*dthp)
            dl2 = dspdth/(2*dthp)
            nl3 = dl2*phib[j,kprev]
            nl4 = dshdth*phibjp1/(sinth*2*dph)
            nl5 = dshdth*phibjm1/(sinth*2*dph)
            ds1 = dshdph/(2*dthp)/sinth
            ns1 = sigmap*phibjp1/(sinth*sinth*dph*dph)
            ns2 = 2*sigmap*phib[j,k]/(sinth*sinth*dph*dph)
            ns3 = sigmap*phibjm1/(sinth*sinth*dph*dph)
            ns4 = dspdph*phibjp1/(sinth*sinth*dph*2*dthp)
            ns5 = dspdph*phibjm1/(sinth*sinth*dph*2*dthp)
            ns6 = dshdph*phib[j,kprev]/(2*dthp)
            
            num = -re*re*ajr[j,k]*1.e-6+ng1+nl1-nl2+nl3-nl4+nl5-ns1+ns2-ns3-ns4+ns5-ns6
            denom = dg1+dl1+dl2-ds1
            
            # #END IF
            if denom != 0:
                phib[j,knext] = num/denom
            #END IF
            if denom == 0:
                phib[j,knext] = phib[j,k]
            #END IF
        # END FOR # end of mlt loop
        
        if( FLG_smoothing_phi_connectMidnight == True ):
            # connect the points across midnight
            phib[23,:] = phib[22,:]+.33*(phib[1,:]-phib[22,:])
            phib[0,:] = phib[22,:]+.6667*(phib[1,:]-phib[22,:])
        #END IF
        if( FLG_smoothing_phi_long == True ):
            # smooth in mlat and mlt
            smphi = uniform_filter1d_wrap(phib[0:23+1,knext],smvallt);
            phib[:,knext] = smphi;
        #END IF
        if( FLG_smoothing_phi_lat == True ):
            for ism in range(0,23+1):
                smth = uniform_filter1d_orig(phib[ism,(k-1):49+1],smvallt); #essential for integration control
                phib[ism,k:49+1] = smth[1:];
            # END FOR
        #END IF
    # END FOR # end of mag lat loop
    
    return phib
#END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def potentializer_dubs(ajr,spb,shb,thmax,thsn,delth,dtr,dph,re,FLG_smoothing_phi_connectMidnight=True,FLG_smoothing_phi_lat=True,FLG_smoothing_phi_long=True):
    phib = np.zeros((24,50),dtype=np.float64);
    phib_top2bot = np.zeros((24,50),dtype=np.float64);
    phib_bot2top = np.zeros((24,50),dtype=np.float64);
    #================tuning parameters================SAVE THIS SET
    maxlat = 48    #maximum latitude for integration    default(best)=48
    smvallt = 3; #default is 3
    meetpt = 72; #latitude to meet at (think continental railroad)
    # #================tuning parameters================THIS SET WORKS EVEN BETTER
    #=================================================
    dthp = -delth*dtr #because the integration goes toward decreasing theta
    
    # if( smvallt == 2 ):
    #     smvallt = 3; #IDL's code only operates on 3 or more so make 2->3, 1 or 0 won't do anything so leave it to nothing on that I guess
    # #END IF
    meetpt -= 40; #adjust to the colatitude range (but it's backwards but don't worry about that)
    # for k in range(1,maxlat+1): # start calculation at 41 degrees
    for k in range(maxlat,meetpt-smvallt//2-1,-1): # start calculation at 88 degrees
        theta = (thmax-thsn[k]+.5)*dtr #0.5 added to center on half-degree mark
        sinth = (np.sin(theta))
        cotth = np.cos(theta)/sinth
        kprev = k+1;
        # kprev = -1; #essentially sets phi prev to 0
        knext = k-1;
        for j in range(0,23+1):
            if j == 0:
                phibjm1 = phib_top2bot[23,k]
                spbjm1 = spb[23,k]
                shbjm1 = shb[23,k]
            else:
                phibjm1 = phib_top2bot[j-1,k]
                spbjm1 = spb[j-1,k]
                shbjm1 = shb[j-1,k]
            #END IF
            if j == 23:
                phibjp1 = phib_top2bot[0,k]
                spbjp1 = spb[0,k]
                shbjp1 = shb[0,k]
            else:
                phibjp1 = phib_top2bot[j+1,k]
                spbjp1 = spb[j+1,k]
                shbjp1 = shb[j+1,k]
            #END IF
                
            sigmap = spb[j,k]
            # sigmah = shb[j,k]
            dspdph=(spbjp1-spbjm1)/2./dph
            dshdph=(shbjp1-shbjm1)/2./dph
            
            dspdth=((spb[j,k+1]-spb[j,k-1])/2./dthp)
            dshdth=((shb[j,k+1]-shb[j,k-1])/2./dthp)
            dg1 = cotth*sigmap/(2*dthp)
            ng1 = dg1*phib_top2bot[j,kprev] #k-1==kprev in orig
            dl1 = sigmap/(dthp*dthp)
            nl1 = 2*sigmap*phib_top2bot[j,k]/(dthp*dthp)
            nl2 = sigmap*phib_top2bot[j,kprev]/(dthp*dthp)
            dl2 = dspdth/(2*dthp)
            nl3 = dl2*phib_top2bot[j,kprev]
            nl4 = dshdth*phibjp1/(sinth*2*dph)
            nl5 = dshdth*phibjm1/(sinth*2*dph)
            ds1 = dshdph/(2*dthp)/sinth
            ns1 = sigmap*phibjp1/(sinth*sinth*dph*dph)
            ns2 = 2*sigmap*phib_top2bot[j,k]/(sinth*sinth*dph*dph)
            ns3 = sigmap*phibjm1/(sinth*sinth*dph*dph)
            ns4 = dspdph*phibjp1/(sinth*sinth*dph*2*dthp)
            ns5 = dspdph*phibjm1/(sinth*sinth*dph*2*dthp)
            ns6 = dshdph*phib_top2bot[j,kprev]/(2*dthp)
            
            num = -re*re*ajr[j,k]*1.e-6+ng1+nl1-nl2+nl3-nl4+nl5-ns1+ns2-ns3-ns4+ns5-ns6
            denom = dg1+dl1+dl2-ds1
            
            # #END IF
            if denom != 0:
                phib_top2bot[j,knext] = num/denom
            #END IF
            if denom == 0:
                phib_top2bot[j,knext] = phib_top2bot[j,k]
            #END IF
        # END FOR # end of mlt loop
        
        if( FLG_smoothing_phi_connectMidnight == True ):
            # connect the points across midnight
            phib_top2bot[23,:] = phib_top2bot[22,:]+.33*(phib_top2bot[1,:]-phib_top2bot[22,:])
            phib_top2bot[0,:] = phib_top2bot[22,:]+.6667*(phib_top2bot[1,:]-phib_top2bot[22,:])
        #END IF
        if( FLG_smoothing_phi_long == True ):
            # smooth in mlat and mlt
            smphi = uniform_filter1d_wrap(phib_top2bot[0:23+1,knext],smvallt);
            phib_top2bot[:,knext] = smphi;
        #END IF
        if( FLG_smoothing_phi_lat == True ):
            for ism in range(0,23+1):
                smth = uniform_filter1d_orig(phib_top2bot[ism,meetpt-smvallt//2:k+1+1],smvallt); #essential for integration control
                phib_top2bot[ism,meetpt-smvallt//2:k+1] = smth[:-1];
            # END FOR
        #END IF
    # END FOR # end of mag lat loop
    for k in range(1,meetpt+smvallt//2+1,1): # start calculation at 41 degrees
        theta = (thmax-thsn[k]+.5)*dtr #0.5 added to center on half-degree mark
        sinth = (np.sin(theta))
        cotth = np.cos(theta)/sinth
        kprev = k+1;
        knext = k-1;
        for j in range(0,23+1):
            if j == 0:
                phibjm1 = phib_bot2top[23,k]
                spbjm1 = spb[23,k]
                shbjm1 = shb[23,k]
            else:
                phibjm1 = phib_bot2top[j-1,k]
                spbjm1 = spb[j-1,k]
                shbjm1 = shb[j-1,k]
            #END IF
            if j == 23:
                phibjp1 = phib_bot2top[0,k]
                spbjp1 = spb[0,k]
                shbjp1 = shb[0,k]
            else:
                phibjp1 = phib_bot2top[j+1,k]
                spbjp1 = spb[j+1,k]
                shbjp1 = shb[j+1,k]
            #END IF
            
            sigmap = spb[j,k]
            # sigmah = shb[j,k]
            dspdph=(spbjp1-spbjm1)/2./dph
            dshdph=(shbjp1-shbjm1)/2./dph
            
            dspdth=((spb[j,k+1]-spb[j,k-1])/2./dthp)
            dshdth=((shb[j,k+1]-shb[j,k-1])/2./dthp)
            dg1 = cotth*sigmap/(2*dthp)
            ng1 = dg1*phib_bot2top[j,kprev] #k-1==kprev in orig
            dl1 = sigmap/(dthp*dthp)
            nl1 = 2*sigmap*phib_bot2top[j,k]/(dthp*dthp)
            nl2 = sigmap*phib_bot2top[j,kprev]/(dthp*dthp)
            dl2 = dspdth/(2*dthp)
            nl3 = dl2*phib_bot2top[j,kprev]
            nl4 = dshdth*phibjp1/(sinth*2*dph)
            nl5 = dshdth*phibjm1/(sinth*2*dph)
            ds1 = dshdph/(2*dthp)/sinth
            ns1 = sigmap*phibjp1/(sinth*sinth*dph*dph)
            ns2 = 2*sigmap*phib_bot2top[j,k]/(sinth*sinth*dph*dph)
            ns3 = sigmap*phibjm1/(sinth*sinth*dph*dph)
            ns4 = dspdph*phibjp1/(sinth*sinth*dph*2*dthp)
            ns5 = dspdph*phibjm1/(sinth*sinth*dph*2*dthp)
            ns6 = dshdph*phib_bot2top[j,kprev]/(2*dthp)
            
            num = -re*re*ajr[j,k]*1.e-6+ng1+nl1-nl2+nl3-nl4+nl5-ns1+ns2-ns3-ns4+ns5-ns6
            denom = dg1+dl1+dl2-ds1
            
            # #END IF
            if denom != 0:
                phib_bot2top[j,knext] = num/denom
            #END IF
            if denom == 0:
                phib_bot2top[j,knext] = phib_bot2top[j,k]
            #END IF
        # END FOR # end of mlt loop
        
        if( FLG_smoothing_phi_connectMidnight == True ):
            # connect the points across midnight
            phib_bot2top[23,:] = phib_bot2top[22,:]+.33*(phib_bot2top[1,:]-phib_bot2top[22,:])
            phib_bot2top[0,:] = phib_bot2top[22,:]+.6667*(phib_bot2top[1,:]-phib_bot2top[22,:])
        #END IF
        if( FLG_smoothing_phi_long == True ):
            # smooth in mlat and mlt
            smphi = uniform_filter1d_wrap(phib_bot2top[0:23+1,knext],smvallt);
            phib_bot2top[:,knext] = smphi;
        #END IF
        if( FLG_smoothing_phi_lat == True ):
            for ism in range(0,23+1):
                smth = uniform_filter1d_orig(phib_bot2top[ism,(k-1):meetpt+1+smvallt//2],smvallt); #essential for integration control
                phib_bot2top[ism,k:meetpt+1+smvallt//2] = smth[1:];
            # END FOR
        #END IF
    # END FOR # end of mag lat loop
    
    phib[:,0:meetpt] = phib_bot2top[:,0:meetpt];
    phib[:,meetpt:] = phib_top2bot[:,meetpt:];
    for ism in range(0,23+1):
        phib[ism,meetpt-1:meetpt+2] = uniform_filter1d_orig(phib[ism,meetpt-1:meetpt+2],3); #smooth out where the two places meet middle bit
    #END FOR ism
    
    return phib
#END DEF
    

@jit(nopython=True,nogil=True,parallel=False,fastmath=False)
def get_auroral_parameters(th50,spb,shb,phib,FLG_smoothing_efield=True):
    # *****************************************************
    dtr = np.pi/180.
    # thmin = 0.
    # thmax = 50.
    # delth = 1.
    # delthkm = delth*109.
    #    arrays that contain the angles in spherical coordinates
    # thsn = delth*np.arange(50,dtype=np.float64)
    dph = 15*dtr
    # phmlt = 15*np.arange(24,dtype=np.float64)
    # mlat=(90.-thmax)+delth*np.arange(50,dtype=np.float64)    #the magnetic latitudes corresponding to the theta array
    re = 6300.e3        # Earth radius in m
    etheta = np.empty((50),dtype=np.float64);
    ephi = np.empty((50),dtype=np.float64);
    enorth = np.empty((24,50),dtype=np.float64);
    eeast = np.empty((24,50),dtype=np.float64);
    jeast = np.zeros((24,50),dtype=np.float64);
    jnorth = np.zeros((24,50),dtype=np.float64);
    jehall = np.empty((24,50),dtype=np.float64);
    jouleheat = np.empty((24,50),dtype=np.float64);
    jrcalc = np.empty((24,50),dtype=np.float64);
    for jlt in range(0,23+1):
        for jlat in range(0,49+1):
        # for jlat in range(48,1-1,-1):
            dth = 1*dtr
            dthkm = (1./360.)*2*np.pi*re/1000.
            dphkm = (15./360.)*2*np.pi*re*np.cos(th50[jlat]*dtr)/1000.
            # calculate electric fields
            dphidth = phib[jlt,jlat]-phib[jlt,jlat-1]    
            # dphidth = phib[jlt,jlat+1]-phib[jlt,jlat]  #goes with 48->0 above
            etheta[jlat] = (1./re)*dphidth/dth         #from high potential to low potential# positive southward in NH
            if jlt == 0:
                phim1 = phib[23,jlat]
            #END IF
            if jlt > 0:
                phim1 = phib[jlt-1,jlat]
            #END IF
            dphidph = phib[jlt,jlat]-phim1
            ephi[jlat] = (1./re)*dphidph/dph     #from high potential to low potential# positive westward
        # END FOR # end of jlat loop
        if( FLG_smoothing_efield == True ):
            # smooth the latitudial profiles of electric field
            ephi = uniform_filter1d_orig(ephi,3);
            etheta = uniform_filter1d_orig(etheta,3);
        #END IF
        # for jlat in range(0,45+1): #OG
        for jlat in range(0,49+1): #unlimited power
        # for jlat in range(49,0-1,-1): #unlimited power [offset issue]
            dphkm = (15./360.)*2*np.pi*re*np.cos(th50[jlat]*dtr)/1000.
            # calculate joule heating rate and currents
            esq = ephi[jlat]*ephi[jlat]+etheta[jlat]*etheta[jlat]
            enorth[jlt,jlat] = -etheta[jlat]*1000.
            eeast[jlt,jlat] = -ephi[jlat]*1000.
            # e-field in V/m# convert Joules/m2-s to ergs/cm2-s by multiplying by 1.00E+03
            # 1 erg is 1.00E-07 Joule# 1.00E+04 cm2 per m2
            jeast[jlt,jlat] = spb[jlt,jlat]*eeast[jlt,jlat]+shb[jlt,jlat]*enorth[jlt,jlat]
            jehall[jlt,jlat] = shb[jlt,jlat]*enorth[jlt,jlat]
            jnorth[jlt,jlat] = spb[jlt,jlat]*enorth[jlt,jlat]-shb[jlt,jlat]*eeast[jlt,jlat]
            jrcalc[jlt,jlat] = -(jeast[jlt,jlat]-jeast[jlt,jlat-1])/dphkm-(jnorth[jlt,jlat]-jnorth[jlt,jlat-1])/dthkm #A/km2 or microamps/m2
            # jrcalc[jlt,jlat] = -(jeast[jlt,jlat+1]-jeast[jlt,jlat])/dphkm-(jnorth[jlt,jlat+1]-jnorth[jlt,jlat])/dthkm #A/km2 or microamps/m2 [goes with 49->0 dir]
            jouleheat[jlt,jlat] = spb[jlt,jlat]*esq*1.e3    # in ergs/cm2-s
            #***************************************************
        # END FOR # end of jlat loop
        # **************************************************
    # END FOR # end of jlt loop
    
    return eeast, enorth, jeast, jnorth, jehall, jouleheat, jrcalc
#END DEF

@jit(nopython=True,nogil=True,parallel=True,fastmath=False)
def optimizer_vingettes(ajr, th50, spb_adj_iterator, shb_adj_iterator, phib_adj_iterator, weighter, goal, rator_iterator, error_iterator, num_vingettes, FLG_smoothing_efield=True):
    for j in prange(0,num_vingettes+3):
        #--- Recalc jrcalc ---
        [_, _, _, _, _, _, jrcalc] = \
            get_auroral_parameters(th50,spb_adj_iterator[:,:,j],shb_adj_iterator[:,:,j],phib_adj_iterator[:,:,j],FLG_smoothing_efield=FLG_smoothing_efield);
        rator_iterator[:,:,j] = np.log10(np.abs(jrcalc*weighter - ajr)/goal); #basically creates standard deviations to use, so automatically good stuff is zeroed in on I hope
        error_iterator[j] = np.sum(np.abs(jrcalc*weighter - ajr));
    #END FOR j
    return rator_iterator, error_iterator
#END DEF

@jit(nopython=True,nogil=True,parallel=False,fastmath=False) #bad math happens with parallel=True (it parallelizes things it shouldn't even w/ prange)
def optimizer(ajr, spb, shb, phib, th50, FLG_smoothing_efield=True, weighter=None):
    np.random.seed(1138);
    if( weighter != None ):
        ajr = np.copy(ajr)*weighter; #weight ajr
    #END IF
    num_tries = 100;
    num_vingettes = 1000;
    #cast magicks on spb/shb/phib all at once b/c jrcalc is based on spb/shb/phib and while phib is based on spb/shb its calc has compound error that makes recalcin it less useful
    #we have a pretty good guess for all of them - so isn't too hard
    goal = 0.00001; #desired difference between jrout (jrcalc) and jrin (ajr)
    limits_sp = np.array( (0,np.max(spb)*2) ); #only positive valued
    limits_sh = np.array( (0,np.max(shb)*2) );
    limits_phi = np.array( (-np.max(np.abs(phib))*2,np.max(np.abs(phib))*2) ); #- and + sadly (more space to check :( )
    [_, _, _, _, _, _, jrcalc] = \
        get_auroral_parameters(th50,spb,shb,phib,FLG_smoothing_efield=FLG_smoothing_efield);
    num_long = spb.shape[0];
    num_lat = spb.shape[1];
    # spb_adj = np.empty(spb.shape,dtype=np.float64);
    # shb_adj = np.empty(shb.shape,dtype=np.float64);
    # phib_adj = np.empty(phib.shape,dtype=np.float64);
    spb_adj = np.copy(spb);
    shb_adj= np.copy(shb);
    phib_adj= np.copy(phib);
    # resetor = np.zeros(spb.shape,dtype=np.bool_); #reset array
    resetor_spb = np.zeros(spb.shape,dtype=np.int32); #reset counter (ticks when limits are hit)
    resetor_shb = np.zeros(spb.shape,dtype=np.int32);
    resetor_phib = np.zeros(spb.shape,dtype=np.int32);
    resetor_limit = 20 #number of times exceeding the limits before a reset hits (shows value hasn't settled)
    # rator_lock = np.zeros(spb.shape,dtype=np.bool_); #locks values if they're close enough [too restricting?]
    
    #--- arrays for multiple loops ---
    rator_iterator = np.empty((spb.shape[0],spb.shape[1],num_vingettes+3),dtype=np.float64);
    spb_adj_iterator = np.empty((spb.shape[0],spb.shape[1],num_vingettes+3),dtype=np.float64);
    shb_adj_iterator = np.empty((spb.shape[0],spb.shape[1],num_vingettes+3),dtype=np.float64);
    phib_adj_iterator = np.empty((spb.shape[0],spb.shape[1],num_vingettes+3),dtype=np.float64);
    error_iterator = np.empty(num_vingettes+3,dtype=np.float64);
    for jk in range(num_vingettes, num_vingettes+3): #num_vingettes is orig, num_vingettes+1 is best case guess, num_vingettes+2 is amalgamation of best cases for each "pixel"
        spb_adj_iterator[:,:,jk] = np.copy(spb);
        shb_adj_iterator[:,:,jk] = np.copy(shb);
        phib_adj_iterator[:,:,jk] = np.copy(phib);
    #END FOR jk
    # error_ref = np.sum(np.abs(jrcalc - ajr)); #reference error from the default code [not needed, always goes down]
    # print(('Error start: ',np.sum(np.abs(jrcalc - ajr))));
    # tic = time.time()
    for i in range(0,num_tries):
        # #--- Recalc jrcalc ---
        # [_, _, _, _, _, _, jrcalc] = \
        #     get_auroral_parameters(th50,spb_adj,shb_adj,phib_adj,FLG_smoothing_efield=FLG_smoothing_efield);
        #--- estimate rators and establish locked values (met the goal) ---
        rator = np.log10(np.abs(jrcalc*weighter - ajr)/goal); #basically creates standard deviations to use, so automatically good stuff is zeroed in on I hope
        # rator_lock = rator_lock | (rator < 0); #locked values
        # rator_lock = (rator < 0); #locked values
        # rator[rator_lock] = 0; #locked in
        rator = np.abs(rator);
        
        #--- Reset to original values if walked off too far ---
        resetor_spb_logical = resetor_spb > resetor_limit;
        resetor_spb_where = np.where(resetor_spb_logical);
        if( np.any(resetor_spb_logical) ):
            for jj in range(0,resetor_spb_where[0].size): #numba did not like 2D logical indexing
                spb_adj[resetor_spb_where[0][jj],resetor_spb_where[1][jj]] = spb[resetor_spb_where[0][jj],resetor_spb_where[1][jj]];
                resetor_spb[resetor_spb_where[0][jj],resetor_spb_where[1][jj]] = 0; #reset the restor counter
            #END FOR jj
            # spb_adj[resetor_spb_where] = spb[resetor_spb_where]; #return to original values
            # resetor_spb[resetor_spb_where] = 0; #reset the restor counter
        #END IF
        resetor_shb_logical = resetor_shb > resetor_limit;
        resetor_shb_where = np.where(resetor_shb_logical);
        if( np.any(resetor_shb > resetor_limit) ):
            for jj in range(0,resetor_shb_where[0].size): #numba did not like 2D logical indexing
                shb_adj[resetor_shb_where[0][jj],resetor_shb_where[1][jj]] = shb[resetor_shb_where[0][jj],resetor_shb_where[1][jj]];
                resetor_shb[resetor_shb_where[0][jj],resetor_shb_where[1][jj]] = 0; #reset the restor counter
            #END FOR jj
            # shb_adj[resetor_shb > resetor_limit] = shb[resetor_shb > resetor_limit];
            # resetor_shb[resetor_shb > resetor_limit] = 0; #reset the restor counter
        #END IF
        resetor_phib_logical = resetor_phib > resetor_limit;
        resetor_phib_where = np.where(resetor_phib_logical);
        if( np.any(resetor_phib > resetor_limit) ):
            for jj in range(0,resetor_phib_where[0].size): #numba did not like 2D logical indexing
                phib_adj[resetor_phib_where[0][jj],resetor_phib_where[1][jj]] = phib[resetor_phib_where[0][jj],resetor_phib_where[1][jj]];
                resetor_phib[resetor_phib_where[0][jj],resetor_phib_where[1][jj]] = 0; #reset the restor counter
            #END FOR jj
            # phib_adj[resetor_phib > resetor_limit] = phib[resetor_phib > resetor_limit];
            # resetor_phib[resetor_phib > resetor_limit] = 0; #reset the restor counter
        #END IF
        # if( np.any(resetor) ):
        #     spb_adj[resetor] = spb[resetor];
        #     shb_adj[resetor] = shb[resetor];
        #     phib_adj[resetor] = phib[resetor];
        #     resetor[:,:] = False; #reset the resetor
        # #END IF
        #--- Generate values near the mean, more accurate values are generated nearer to the mean ---
        for jlong in range(0,num_long):
            for klat in range(0,num_lat):
                spb_adj_iterator[jlong,klat,0:num_vingettes] = np.random.normal(loc=spb_adj[jlong,klat],scale=rator[jlong,klat],size=num_vingettes);
                shb_adj_iterator[jlong,klat,0:num_vingettes] = np.random.normal(loc=shb_adj[jlong,klat],scale=rator[jlong,klat],size=num_vingettes);
                phib_adj_iterator[jlong,klat,0:num_vingettes] = np.random.normal(loc=phib_adj[jlong,klat],scale=rator[jlong,klat],size=num_vingettes);
            #END FOR klat
        #END FOR jlong
        #--- Enforce limits and keep a counter if values exceed the reset limit ---
        kj = (spb_adj_iterator > limits_sp[1]) | (spb_adj_iterator < limits_sp[0]);
        resetor_spb += kj[:,:,0];
        # spb_adj_iterator[kj] = np.random.uniform(low=limits_sp[0],high=limits_sp[1],size=np.sum(kj)); #numba idd not like
        kj_where = np.where(kj);
        for jj in range(0,kj.sum()):
            spb_adj_iterator[kj_where[0][jj],kj_where[1][jj],kj_where[2][jj]] = np.random.uniform(limits_sp[0],limits_sp[1]);
        #END FOR jj
        kj = (shb_adj_iterator > limits_sh[1]) | (shb_adj_iterator < limits_sh[0]);
        resetor_shb += kj[:,:,0];
        # shb_adj_iterator[kj] = np.random.uniform(low=limits_sh[0],high=limits_sh[1],size=np.sum(kj));
        kj_where = np.where(kj);
        for jj in range(0,kj.sum()):
            shb_adj_iterator[kj_where[0][jj],kj_where[1][jj],kj_where[2][jj]] = np.random.uniform(limits_sp[0],limits_sp[1]);
        #END FOR jj
        kj = (phib_adj_iterator > limits_phi[1]) | (phib_adj_iterator < limits_phi[0]);
        resetor_phib += kj[:,:,0];
        # phib_adj_iterator[kj] = np.random.uniform(low=limits_phi[0],high=limits_phi[1],size=np.sum(kj));
        kj_where = np.where(kj);
        for jj in range(0,kj.sum()):
            phib_adj_iterator[kj_where[0][jj],kj_where[1][jj],kj_where[2][jj]] = np.random.uniform(limits_sp[0],limits_sp[1]);
        #END FOR jj
        #--- Keep locked values locked (redundant with rator[rator_lock] = 0) ---
        # spb_adj[rator_lock] = spb[rator_lock];
        # shb_adj[rator_lock] = shb[rator_lock];
        # phib_adj[rator_lock] = phib[rator_lock]; 
        
        rator_iterator, error_iterator = optimizer_vingettes(ajr, th50, spb_adj_iterator, shb_adj_iterator, phib_adj_iterator, weighter, \
                            goal, rator_iterator, error_iterator, num_vingettes, FLG_smoothing_efield=FLG_smoothing_efield); #call the vingette tester optimizer
        
        #--- figure out the best ---
        jk = np.where(np.min(error_iterator) == error_iterator)[0][0];
        spb_adj[:,:] = np.copy(spb_adj_iterator[:,:,jk]);
        shb_adj[:,:] = np.copy(shb_adj_iterator[:,:,jk]);
        phib_adj[:,:] = np.copy(phib_adj_iterator[:,:,jk]);
        
        # if( jk < num_vingettes ):
        #     print('Guessin\' got new best case!')
        # elif( jk == num_vingettes ):
        #     print('Original is still best. It\'s a reset.')
        # elif( jk == (num_vingettes+1) ):
        #     print('Last best case is still best - no improvement.');
        # elif( jk == (num_vingettes+2) ):
        #     print('Rator amalgomation is new best case!');
        # #END IF
        
        #--- record best case error found (+1 spot) ---
        # jk = np.where(np.min(error_iterator) == error_iterator)[0][0];
        spb_adj_iterator[:,:,num_vingettes+1] = np.copy(spb_adj_iterator[:,:,jk]);
        shb_adj_iterator[:,:,num_vingettes+1] = np.copy(shb_adj_iterator[:,:,jk]);
        phib_adj_iterator[:,:,num_vingettes+1] = np.copy(phib_adj_iterator[:,:,jk]);
        
        #--- create amalgomation iterator (+2 spot) ---
        for jlong in range(0,num_long):
            for klat in range(0,num_lat):
                # rator_iterator_min = np.min(rator_iterator[jlong,klat,:]);
                # if( rator_iterator_min < rator[jlong,klat] ):
                jk = np.where(np.min(rator_iterator[jlong,klat,:]) == rator_iterator[jlong,klat,:])[0][0];
                spb_adj_iterator[jlong,klat,num_vingettes+2] = spb_adj_iterator[jlong,klat,jk];
                shb_adj_iterator[jlong,klat,num_vingettes+2] = shb_adj_iterator[jlong,klat,jk];
                phib_adj_iterator[jlong,klat,num_vingettes+2] = phib_adj_iterator[jlong,klat,jk];
                # else:
                #     resetor[jlong,klat] = True; #reset it
                #END IF
            #END FOR klat
        #END FOR jlong
        
        #--- Recalc jrcalc ---
        [_, _, _, _, _, _, jrcalc] = \
            get_auroral_parameters(th50,spb_adj,shb_adj,phib_adj,FLG_smoothing_efield=FLG_smoothing_efield);
        # print('Error '+str(i)+'/'+str(num_tries)+': '+str(np.sum(np.abs(jrcalc - ajr))));
        # print(np.sum(np.abs(jrcalc - ajr)))
    #END FOR i
    spb = spb_adj;
    shb = shb_adj;
    phib = phib_adj;
    # print('Error end: '+str(np.sum(np.abs(jrcalc - ajr)))+', took: '+str(time.time()-tic)+' sec');
    
    return spb, shb, phib
#END DEF