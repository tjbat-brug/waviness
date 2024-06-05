import numpy as np       
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt 
from datetime import datetime as dt

# Sinuosity index specific
import haversine as hv

# For plotting
from matplotlib import rc


def is_circumglobal(segment):
    '''
    Check if the first entry is below 0.5 and the last entry is above 359
    or if the first entry is above 359 and the last entry is below 0.5
    '''
    first_lon, _ = segment[0]
    last_lon, _ = segment[-1]
    
    return (first_lon < 0.5 and last_lon > 359) or (first_lon > 359 and last_lon < 0.5)

def group_connected_segments(lst1, lst2):
    connected_segments = {}
    used_numbers = set()  # To keep track of which numbers have been assigned to a group
    for num1, num2 in zip(lst1, lst2):
        if num1 in used_numbers and num2 in used_numbers:
            # Both numbers already in groups, do nothing
            continue
        elif num1 in used_numbers:
            # Add num2 to the group containing num1
            group = connected_segments[num1]
            group.add(num2)
            connected_segments[num2] = group
            used_numbers.add(num2)
        elif num2 in used_numbers:
            # Add num1 to the group containing num2
            group = connected_segments[num2]
            group.add(num1)
            connected_segments[num1] = group
            used_numbers.add(num1)
        else:
            # Both numbers not in any group, create new group
            group = {num1, num2}
            connected_segments[num1] = group
            connected_segments[num2] = group
            used_numbers.update(group)
    # Gather unique groups
    unique_groups = {id(group): group for group in connected_segments.values()}
    return [sorted(group) for group in unique_groups.values()]

def lengths_for_SI(segments, resol_thresh, conn_thresh):

    TOTAL_LENGTH = []                               # This list will contain the length of the different segments 
    LENGTH_CIRCUMGLOBAL = []                        # This list will contain the length of the circumglobal isohypse
    LATS_CIRCUMGLOBAL = []                          # This list will contain the latitudes of the circumglobal isohypse
    LENGTH_CUT_OFFS = []                            # This list will contain the length of the cut-off segments
    LATS_CUT_OFFS = []                              # This list will contain the latitudes of the cut-off segments
    NUMBER_CUT_OFFS = 0                             # This list will contain the number of cut-off segments
    TOTAL_LATS = []                                 # This list will contain the latitudes of the different segments 
    
    for i_s, segment in enumerate(segments[1]):     # loop over all segments
        total_length = 0.
        for i_p in range(segment.shape[0]):         # loop over all points in a segment
            ''' 
            The haversine function below calculates the distance on Earth between to different lat/lon points in km. 
            ax1.contour gives lon, lat, but haversine wants lat, lon, therefore the [::-1] is used to reverse the order of the points.
            '''
            length = hv.haversine(segment[i_p][::-1], segment[i_p-1][::-1], normalize=True)     
            ''' 
            If the distance between two points is larger than 150 km, it is assumed that the distance between the last point of the segment 
            and the first point of the segment that crosses the 0/360 meridian
            '''
            if length >= 150:                       
                    continue
            else:
                total_length += length                      # The total length of the segment is calculated
            TOTAL_LATS.append(segment[i_p][1])              # The latitude of the point is added to the list (segment[i_p] = [lon, lat])
        TOTAL_LENGTH.append(total_length)                   # The total length of the segment is added to the list
    
    longest_segment = segments[1][int(np.where(np.max(TOTAL_LENGTH) == TOTAL_LENGTH)[0])] # The longest segment is selected

    '''
    The code below checks which segments are circumglobal and which segments are cut-off segments.
    '''

    circumglobal_longest = False
    length_cutoffs = []
    connected_1 = []
    connected_2 = []

    for i_s, segment in enumerate(segments[1]):
        print('Segment: ', i_s, TOTAL_LENGTH[i_s], np.mean(segment[np.isin(segment[:, 1], TOTAL_LATS)][:, 1]), np.mean(segments[1][i_s][:,0]))
        '''
        Loop over segments
        '''
        if segment is longest_segment:
            ''' 
            If the longest segment is circumglobal circumglobal_longest is set to True
            '''
            if is_circumglobal(longest_segment):
                circumglobal_longest = True
                continue
            else:
                print(i_s, 'Longest segment is strangly not circumglobal') #
                continue

        elif segment[0][0] == 0.0 or segment[-1][0] == 359.296875 or segment[0][0] == 359.296875 or segment[-1][0] == 0.0:
            '''
            If the segment crosses the 0/360 meridian and is not circumglobal, 
            it is checked if it is connected to the longest segment or a cut-off segment.
            '''
            connected = 0
            for i_k, k in enumerate(segments[1]): # loop over all segments
                if k is segment:
                    continue
                elif hv.haversine(k[0][::-1], segment[0][::-1], normalize=True) <= conn_thresh or \
                        hv.haversine(k[-1][::-1], segment[-1][::-1], normalize=True) <= conn_thresh or \
                        hv.haversine(k[0][::-1], segment[-1][::-1], normalize=True) <= conn_thresh or \
                        hv.haversine(k[-1][::-1], segment[0][::-1], normalize=True) <= conn_thresh:
                    # print(hv.haversine(k[0][::-1], l[0][::-1], normalize=True), hv.haversine(k[-1][::-1], l[-1][::-1], normalize=True),\
                    #       hv.haversine(k[0][::-1], l[-1][::-1], normalize=True), hv.haversine(k[-1][::-1], l[0][::-1], normalize=True))
                    connected += 1
                    connected_1.append(i_s)
                    connected_2.append(i_k)
        else:
            if TOTAL_LENGTH[i_s] <= resol_thresh:
                TOTAL_LENGTH[i_s] = 0
                length_cutoffs.append(TOTAL_LENGTH[i_s])
                continue
            else:
                NUMBER_CUT_OFFS += 1
                length_cutoffs.append(TOTAL_LENGTH[i_s])
                LATS_CUT_OFFS.extend(segment[np.isin(segment[:, 1], TOTAL_LATS)][:, 1])

    connected_segments = group_connected_segments(connected_1, connected_2)

    # print('circumglobal longest: ', circumglobal_longest)
    
    for connected_segment in connected_segments:
        length_connected_segment = []
        for i in connected_segment:
            length_connected_segment.append(TOTAL_LENGTH[i])

        if len(connected_segment) == 2 and np.max(length_connected_segment) != np.max(TOTAL_LENGTH) and circumglobal_longest == True:
            NUMBER_CUT_OFFS += 1 
            
            length_cutoffs.append(np.sum(length_connected_segment))
            LATS_CUT_OFFS.extend(segments[1][connected_segment[0]][np.isin(segments[1][connected_segment[0]][:, 1], TOTAL_LATS)][:, 1])
            LATS_CUT_OFFS.extend(segments[1][connected_segment[1]][np.isin(segments[1][connected_segment[1]][:, 1], TOTAL_LATS)][:, 1])
            
        elif len(connected_segment) == 2 and np.max(length_connected_segment) == np.max(TOTAL_LENGTH) and circumglobal_longest == False:
            '''
            We assume that the longest segment is part of circumglobal isohypse 
            but that the circumglobal isohypse crosses the pole and therefore circumglobal = False
            '''
            LENGTH_CIRCUMGLOBAL.append(np.sum(length_connected_segment))
            LATS_CIRCUMGLOBAL.extend(segments[1][connected_segment[0]][np.isin(segments[1][connected_segment[0]][:, 1], TOTAL_LATS)][:, 1])
            LATS_CIRCUMGLOBAL.extend(segments[1][connected_segment[1]][np.isin(segments[1][connected_segment[1]][:, 1], TOTAL_LATS)][:, 1])

        elif len(connected_segment) == 2 and np.max(length_connected_segment) == np.max(TOTAL_LENGTH) and circumglobal_longest == True:
            LENGTH_CIRCUMGLOBAL.append(np.sum(length_connected_segment))
            LATS_CIRCUMGLOBAL.extend(segments[1][connected_segment[0]][np.isin(segments[1][connected_segment[0]][:, 1], TOTAL_LATS)][:, 1])
            LATS_CIRCUMGLOBAL.extend(segments[1][connected_segment[1]][np.isin(segments[1][connected_segment[1]][:, 1], TOTAL_LATS)][:, 1])
        
        elif len(connected_segment) > 2:
            if np.max(length_connected_segment) == np.max(TOTAL_LENGTH):
                LENGTH_CIRCUMGLOBAL.append(np.sum(length_connected_segment))
                for i in connected_segment:
                    LATS_CIRCUMGLOBAL.extend(segments[1][i][np.isin(segments[1][i][:, 1], TOTAL_LATS)][:, 1])
            elif np.max(length_connected_segment) != np.max(TOTAL_LENGTH) and circumglobal_longest:
                NUMBER_CUT_OFFS += 1 
                length_cutoffs.append(np.sum(length_connected_segment))
                for i in connected_segment:
                    LATS_CUT_OFFS.extend(segments[1][i][np.isin(segments[1][i][:, 1], TOTAL_LATS)][:, 1])
            else:
                continue
    
    if circumglobal_longest and len(LENGTH_CIRCUMGLOBAL) == 0:
        LENGTH_CIRCUMGLOBAL.append(TOTAL_LENGTH[int(np.where(np.max(TOTAL_LENGTH) == TOTAL_LENGTH)[0])])
        LATS_CIRCUMGLOBAL.extend(segments[1][int(np.where(np.max(TOTAL_LENGTH) == TOTAL_LENGTH)[0])][np.isin(segments[1][int(np.where(np.max(TOTAL_LENGTH) == TOTAL_LENGTH)[0])][:, 1], TOTAL_LATS)][:, 1])
        
    # Sums of the different lists are calculated
    TOTAL_LENGTH = np.sum(TOTAL_LENGTH)
    if len(LENGTH_CIRCUMGLOBAL) == 0:
        LENGTH_CIRCUMGLOBAL = np.nan
    LENGTH_CIRCUMGLOBAL = np.sum(LENGTH_CIRCUMGLOBAL)
    if len(LATS_CIRCUMGLOBAL) == 0:
        LATS_CIRCUMGLOBAL = np.nan
    if len(length_cutoffs) == 0:
        LENGTH_CUT_OFFS = np.nan
    if len(LATS_CUT_OFFS) == 0:
        LATS_CUT_OFFS = np.nan
    if np.sum(length_cutoffs) > resol_thresh:
        # print(f'Cut-off segments larger than {round(resol_thresh)} km')
        LENGTH_CUT_OFFS = np.sum(length_cutoffs)
    elif np.sum(length_cutoffs) <= resol_thresh:
        # print(f'Cut-off segments smaller than {round(resol_thresh)} km')
        LENGTH_CUT_OFFS = np.nan
        TOTAL_LENGTH = TOTAL_LENGTH - np.sum(length_cutoffs)
    
    plt.show()
    
    return TOTAL_LENGTH, LENGTH_CIRCUMGLOBAL, LENGTH_CUT_OFFS, NUMBER_CUT_OFFS, TOTAL_LATS, LATS_CIRCUMGLOBAL, LATS_CUT_OFFS

### Data parameters ###
experiments = ['CNTRL', 'SST4', 'RTG', 'PA'] 
path = "path_to_data"
plev = '500'

### Thresholds ###
resol_thresh =  2 * np.pi * 78  # 2 * pi * 78 km (= T255, at the equator)
conn_thresh = 3 * 78            # Distance of ~3 gridpoints


for exp in experiments:
    if exp == 'CNTRL':
        expid = 'h1zz'
    elif exp == 'SST4':
        expid = 'h4zz'
    elif exp == 'RTG':
        expid = 'h1aa'
    elif exp == 'PA':
        expid = 'h2aa'
    else:
        print('Experiment not found')
        
    print('Experiment: ' + exp + '\n')
        
    zonmean = xr.open_dataset(path + exp + "/" + exp + "_" + expid + "_SH_timzonmean.nc")

    print('Adjusted SI-method lat range based on individual averages.\n')

    # Take time mean zonal mean to calculate the threshold and latitude range for less arbitrary threshold

    u500_zonmean = zonmean['U'][0,:,:,0].sel(plev=50000)
    v500_zonmean = zonmean['V'][0,:,:,0].sel(plev=50000)
    V500_zonmean = (u500_zonmean**2 + v500_zonmean**2)**.5
    
    print('Latitude ranges are calculated for the different hemisphere at different level' +  '\n')

    # Create or open a text file for output data
    file =  open('Output_SI/' + exp + '_output_yearly.txt', 'w') 

    file1 = open('Output_SI/' + exp + "/" + exp + '_Z' + plev + '_NH_30_70.txt', 'w')
    file2 = open('Output_SI/' + exp + "/" + exp + '_Z' + plev + '_SH_30_70.txt', 'w')
    file3 = open('Output_SI/' + exp + "/" + exp + '_Z' + plev + '_NH_mod.txt', 'w')
    file4 = open('Output_SI/' + exp + "/" + exp + '_Z' + plev + '_SH_mod.txt', 'w')

    file1.write('Date,Z' + plev + '(m),Z' + plev + '_length_(km),SI,Avelat(deg),Mxlat(deg),Mnlat(deg),Merext(deg),Z' + plev + '_circumglobal,SI_circumglobal,Avelat_circumglobal(deg),Maxlat_circumglobal(deg),Minlat_circumglobal(deg),Merext_circumglobal(deg),Z' + plev + '_cutoff,SI_cutoff,Avelat_cutoff(deg),Maxlat_cutoff(deg),Minlat_cutoff(deg),Merext_cutoff(deg),number_cutoff,\n')
    file2.write('Date,Z' + plev + '(m),Z' + plev + '_length_(km),SI,Avelat(deg),Mxlat(deg),Mnlat(deg),Merext(deg),Z' + plev + '_circumglobal,SI_circumglobal,Avelat_circumglobal(deg),Maxlat_circumglobal(deg),Minlat_circumglobal(deg),Merext_circumglobal(deg),Z' + plev + '_cutoff,SI_cutoff,Avelat_cutoff(deg),Maxlat_cutoff(deg),Minlat_cutoff(deg),Merext_cutoff(deg),number_cutoff,\n')
    file3.write('Date,Z' + plev + '(m),Z' + plev + '_length_(km),SI,Avelat(deg),Mxlat(deg),Mnlat(deg),Merext(deg),Z' + plev + '_circumglobal,SI_circumglobal,Avelat_circumglobal(deg),Maxlat_circumglobal(deg),Minlat_circumglobal(deg),Merext_circumglobal(deg),Z' + plev + '_cutoff,SI_cutoff,Avelat_cutoff(deg),Maxlat_cutoff(deg),Minlat_cutoff(deg),Merext_cutoff(deg),number_cutoff,\n')
    file4.write('Date,Z' + plev + '(m),Z' + plev + '_length_(km),SI,Avelat(deg),Mxlat(deg),Mnlat(deg),Merext(deg),Z' + plev + '_circumglobal,SI_circumglobal,Avelat_circumglobal(deg),Maxlat_circumglobal(deg),Minlat_circumglobal(deg),Merext_circumglobal(deg),Z' + plev + '_cutoff,SI_cutoff,Avelat_cutoff(deg),Maxlat_cutoff(deg),Minlat_cutoff(deg),Merext_cutoff(deg),number_cutoff,\n')
    
    ############################################################################################################
    for y in range(2001,2021,1):
        
        ########################
        ### Loading the data ###
        ########################

        data = xr.open_dataset(path + exp + "/" + exp + "_" + expid + "_" + str(y) + ".nc")

        GPH = data['Z'][:,:,:,:].sel(plev=int(plev) * 100) / 9.81 
        time = GPH.time
        lon = GPH.lon
        lat = GPH.lat
        latr = zonmean.lat.values # needed to find lat ranges

        date_times = []

        for j in range(len(time)):  
            x = time.values[j]
            year = int(str(x)[:4])
            month = int(str(x)[4:6])
            day = int(str(x)[6:8])
            if int(int(str(x)[9:12])) == 5:
                hour = 12
            else:
                hour = int(int(str(x)[9:12]) / 100 * 24)
            date_times.append(dt(year=year, month=month, day=day, hour=hour))

        print("Data is loaded for year " + date_times[0].strftime("%Y") + "\n")

        # Reference latitude
        a_earth = 6.371e3
        ref_50  = 2 * np.pi * a_earth * np.cos(np.radians(50))

        print("Parameters are set \n")

        ### Averaging GPH to calculate the mean isohypse value for a certain latitude band for both hemispheres
        GPH_mean_30_70 = np.full([len(time),2], np.nan) # Column one is NH, column two is SH
        GPH_mean_mod = np.full([len(time),2], np.nan)
        
        '''
        Latitude Range for V500 (0.5 x max(V500))
        '''
        
        ind = int(len(latr)/2)
        threshold = V500_zonmean.max() * 0.50

        lat_max_500_NH = latr[:ind][V500_zonmean.where(V500_zonmean > threshold)[:ind].squeeze().notnull()].max()
        lat_min_500_NH = latr[:ind][V500_zonmean.where(V500_zonmean > threshold)[:ind].squeeze().notnull()].min()
        lat_max_500_SH = latr[ind:][V500_zonmean.where(V500_zonmean > threshold)[ind:].squeeze().notnull()].max()
        lat_min_500_SH = latr[ind:][V500_zonmean.where(V500_zonmean > threshold)[ind:].squeeze().notnull()].min()

        print("50% V500, lat ranges")
        print(np.round(threshold.values, 2),
                np.round(lat_max_500_NH, 2),
                np.round(lat_min_500_NH, 2),
                np.round(lat_max_500_SH, 2),
                np.round(lat_min_500_SH, 2),
                )

        for ti in range(len(time)):
            # original SI
            GPH_mean_30_70[ti,0] = GPH[ti,:,:].sel(lat=slice(70,30)).values.mean()  
            GPH_mean_30_70[ti,1] = GPH[ti,:,:].sel(lat=slice(-30,-70)).values.mean()
            
            # modified SI based on 0.5maxV500
            GPH_mean_mod[ti,0] = GPH[ti,:,:].sel(lat=slice(lat_max_500_NH,lat_min_500_NH)).values.mean()
            GPH_mean_mod[ti,1] = GPH[ti,:,:].sel(lat=slice(lat_max_500_SH,lat_min_500_SH)).values.mean()
                
        print("Mean Isohypses are calculated \n")

        ############################################################################################################
        ########################### This block calculates the SI values ############################################
        ############################################################################################################

        ### Creating the data array's ###
        # 1 = total length (incl. cut off lows and highs), 
        # 2 = SI, 
        # 3 = average lat
        # 4 = max lat
        # 5 = min lat
        # 6 = meridional extent
        # 7 = length circumglobal isohypse segment
        # 8 = SI circumglobal isohypse segment
        # 9 = mean lat longest segment
        # 10 = max lat longest segment
        # 11 = min lat longest segment
        # 12 = meridional extent circumglobal isohypse
        # 13 = length Cut-off segments
        # 14 = SI Cut-off segements
        # 15 = mean lat Cut-off segments
        # 16 = max lat Cut-off segments
        # 17 = min lat Cut-off segments
        # 18 = meridional extent Cut-off segments
        # 19 = number of cut-offs

        #################################

        # Original SI metric by Cattiaux et al. (2016)
        exp_NH_30_70 = np.full([len(time),19], np.nan)
        exp_SH_30_70 = np.full([len(time),19], np.nan)
        
        # Adjusted SI metric based on individual averages and V500, normalised by circumference of the Earth at 40 degrees
        exp_NH_mod = np.full([len(time),19], np.nan)
        exp_SH_mod = np.full([len(time),19], np.nan)

        ############################################################################################################
        ############################# Calculating SI for every timestep ############################################
        ############################################################################################################

        for i in range(len(time)):
            print("Calculate for: " + date_times[i].strftime("%Y%m%d%H"))

            ### Making it possilbe to calculate the isohypse length by creating a figure object 
            fig1 = plt.figure()
            ax1 = plt.subplot(111)

            ### Selecting the mean GPH isohypses for each hemisphere ###
            NH_30_70 = np.array([-1e6, GPH_mean_30_70[i][0], 1e6])
            SH_30_70 = np.array([-1e6, GPH_mean_30_70[i][1], 1e6])
            NH_mod = np.array([-1e6, GPH_mean_mod[i][0], 1e6])
            SH_mod = np.array([-1e6, GPH_mean_mod[i][1], 1e6])

            cz_NH_30_70 = ax1.contour(lon, lat.sel(lat=slice(90,0)),  GPH[i,:,:].sel(lat=slice(90,0)),  NH_30_70) 
            cz_SH_30_70 = ax1.contour(lon, lat.sel(lat=slice(0,-90)), GPH[i,:,:].sel(lat=slice(0,-90)), SH_30_70)
            cz_NH_mod = ax1.contour(lon, lat.sel(lat=slice(90,0)),  GPH[i,:,:].sel(lat=slice(90,0)),  NH_mod) 
            cz_SH_mod = ax1.contour(lon, lat.sel(lat=slice(0,-90)), GPH[i,:,:].sel(lat=slice(0,-90)), SH_mod)

            ### 30 to 70 latitude  average, original SI metric ###
            # Be careful with how longitude and latitude comes out. On Snellius it is lon, lat and haversine wants to have lat, lon
            lines_NH = cz_NH_30_70.allsegs[:]  
            lines_SH = cz_SH_30_70.allsegs[:] 

            ### Northern Hemisphere ###
            TOTAL_LENGTH, LENGTH_CIRCUMGLOBAL, LENGTH_CUT_OFFS, NUMBER_CUT_OFFS, TOTAL_LATS, LATS_CIRCUMGLOBAL, LATS_CUT_OFFS = lengths_for_SI(lines_NH, resol_thresh, conn_thresh)

            exp_NH_30_70[i][0] = TOTAL_LENGTH                                      # Total length GPH including cut-off segments
            exp_NH_30_70[i][1] = TOTAL_LENGTH / ref_50                             # SI metric, normalising with circumference Earth at 50 latitude
            exp_NH_30_70[i][2] = np.mean(TOTAL_LATS)                                    # mean lat of all points of the isohypse (this is right, because it takes the mean of all points and not mean of means)
            exp_NH_30_70[i][3] = np.max(TOTAL_LATS)                                     # max lat of all points of the isohypse
            exp_NH_30_70[i][4] = np.min(TOTAL_LATS)                                     # min lat of all points of the isohypse
            exp_NH_30_70[i][5] = np.max(TOTAL_LATS) - np.min(TOTAL_LATS)                # meridional extent of the isohypse
            exp_NH_30_70[i][6] = LENGTH_CIRCUMGLOBAL                                    # length circumglobal segment
            exp_NH_30_70[i][7] = LENGTH_CIRCUMGLOBAL / ref_50                           # SI longest segment
            exp_NH_30_70[i][8] = np.mean(LATS_CIRCUMGLOBAL)                             # mean lat longest segment     
            exp_NH_30_70[i][9] = np.max(LATS_CIRCUMGLOBAL)                              # max lat longest segment
            exp_NH_30_70[i][10] = np.min(LATS_CIRCUMGLOBAL)                             # min lat longest segment
            exp_NH_30_70[i][11] = np.max(LATS_CIRCUMGLOBAL) - np.min(LATS_CIRCUMGLOBAL) # meridional extent longest segment
            exp_NH_30_70[i][12] = LENGTH_CUT_OFFS                                       # length cut-off segments
            exp_NH_30_70[i][13] = LENGTH_CUT_OFFS / ref_50                              # SI cut-off segments
            exp_NH_30_70[i][14] = np.mean(LATS_CUT_OFFS)                                # mean lat cut-off segments
            exp_NH_30_70[i][15] = np.max(LATS_CUT_OFFS)                                 # max lat cut-off segments
            exp_NH_30_70[i][16] = np.min(LATS_CUT_OFFS)                                 # min lat cut-off segments
            exp_NH_30_70[i][17] = np.max(LATS_CUT_OFFS) - np.min(LATS_CUT_OFFS)         # meridional extent cut-off segments
            exp_NH_30_70[i][18] = NUMBER_CUT_OFFS                                       # number of cut-offs

            ### Southern Hemisphere ###
            TOTAL_LENGTH, LENGTH_CIRCUMGLOBAL, LENGTH_CUT_OFFS, NUMBER_CUT_OFFS, TOTAL_LATS, LATS_CIRCUMGLOBAL, LATS_CUT_OFFS = lengths_for_SI(lines_SH, resol_thresh, conn_thresh)

            exp_SH_30_70[i][0] = TOTAL_LENGTH                                          # Total length GPH including cut-off segments
            exp_SH_30_70[i][1] = TOTAL_LENGTH / ref_50                                 # SI metric, normalising with circumference Earth at 50 latitude
            exp_SH_30_70[i][2] = np.mean(TOTAL_LATS)                                        # mean lat of all points of the isohypse (this is right, because it takes the mean of all points and not mean of means)
            exp_SH_30_70[i][3] = np.min(TOTAL_LATS)                                         # max lat of all points of the isohypse (note: Southern Hemisphere)
            exp_SH_30_70[i][4] = np.max(TOTAL_LATS)                                         # min lat of all points of the isohypse (note: Southern Hemisphere)
            exp_SH_30_70[i][5] = np.abs(np.min(TOTAL_LATS)) - np.abs(np.min(TOTAL_LATS))    # meridional extent of the isohypse
            exp_SH_30_70[i][6] = LENGTH_CIRCUMGLOBAL                                        # length circumglobal segment
            exp_SH_30_70[i][7] = LENGTH_CIRCUMGLOBAL / ref_50                               # SI longest segment
            exp_SH_30_70[i][8] = np.mean(LATS_CIRCUMGLOBAL)                                 # mean lat longest segment
            exp_SH_30_70[i][9] = np.min(LATS_CIRCUMGLOBAL)                                  # max lat longest segment
            exp_SH_30_70[i][10] = np.max(LATS_CIRCUMGLOBAL)                                 # min lat longest segment
            exp_SH_30_70[i][11] = np.abs(np.max(LATS_CIRCUMGLOBAL)) - np.abs(np.min(LATS_CIRCUMGLOBAL))     # meridional extent longest segment
            exp_SH_30_70[i][12] = LENGTH_CUT_OFFS                                           # length cut-off segments
            exp_SH_30_70[i][13] = LENGTH_CUT_OFFS / ref_50                                  # SI cut-off segments
            exp_SH_30_70[i][14] = np.mean(LATS_CUT_OFFS)                                    # mean lat cut-off segments
            exp_SH_30_70[i][15] = np.min(LATS_CUT_OFFS)                                     # max lat cut-off segments
            exp_SH_30_70[i][16] = np.max(LATS_CUT_OFFS)                                     # min lat cut-off segments
            exp_SH_30_70[i][17] = np.abs(np.min(LATS_CUT_OFFS)) - np.abs(np.max(LATS_CUT_OFFS))         # meridional extent cut-off segments
            exp_SH_30_70[i][18] = NUMBER_CUT_OFFS                                           # number of cut-offs

            ### Modified, normalize over mean latitude of the isohypse
            lines_NH = cz_NH_mod.allsegs[:]  
            lines_SH = cz_SH_mod.allsegs[:] 

            ### Northern Hemisphere ###
            TOTAL_LENGTH, LENGTH_CIRCUMGLOBAL, LENGTH_CUT_OFFS, NUMBER_CUT_OFFS, TOTAL_LATS, LATS_CIRCUMGLOBAL, LATS_CUT_OFFS = lengths_for_SI(lines_NH, resol_thresh, conn_thresh)

            exp_NH_mod[i][0] = TOTAL_LENGTH
            exp_NH_mod[i][1] = TOTAL_LENGTH / (2 * np.pi * a_earth * np.cos(np.radians(np.mean(TOTAL_LATS)))) # normalise with mean lat
            exp_NH_mod[i][2] = np.mean(TOTAL_LATS)
            exp_NH_mod[i][3] = np.max(TOTAL_LATS)
            exp_NH_mod[i][4] = np.min(TOTAL_LATS)
            exp_NH_mod[i][5] = np.max(TOTAL_LATS) - np.min(TOTAL_LATS)
            exp_NH_mod[i][6] = LENGTH_CIRCUMGLOBAL                                  # length circumglobal segment
            exp_NH_mod[i][7] = LENGTH_CIRCUMGLOBAL / (2 * np.pi * a_earth * np.cos(np.radians(np.mean(LATS_CIRCUMGLOBAL)))) # normalise with mean lat
            exp_NH_mod[i][8] = np.mean(LATS_CIRCUMGLOBAL)                           # mean lat longest segment
            exp_NH_mod[i][9] = np.max(LATS_CIRCUMGLOBAL)                            # max lat longest segment
            exp_NH_mod[i][10] = np.min(LATS_CIRCUMGLOBAL)                           # min lat longest segment
            exp_NH_mod[i][11] = np.max(LATS_CIRCUMGLOBAL) - np.min(LATS_CIRCUMGLOBAL) # meridional extent longest segment
            exp_NH_mod[i][12] = LENGTH_CUT_OFFS                                     # length cut-off segments
            exp_NH_mod[i][13] = LENGTH_CUT_OFFS / (2 * np.pi * a_earth * np.cos(np.radians(np.mean(LATS_CUT_OFFS)))) # normalise with mean lat
            exp_NH_mod[i][14] = np.mean(LATS_CUT_OFFS)                              # mean lat cut-off segments
            exp_NH_mod[i][15] = np.max(LATS_CUT_OFFS)                               # max lat cut-off segments
            exp_NH_mod[i][16] = np.min(LATS_CUT_OFFS)                               # min lat cut-off segments
            exp_NH_mod[i][17] = np.max(LATS_CUT_OFFS) - np.min(LATS_CUT_OFFS)       # meridional extent cut-off segments
            exp_NH_mod[i][18] = NUMBER_CUT_OFFS                                     # number of cut-offs

            ### Southern Hemisphere ###
            TOTAL_LENGTH, LENGTH_CIRCUMGLOBAL, LENGTH_CUT_OFFS, NUMBER_CUT_OFFS, TOTAL_LATS, LATS_CIRCUMGLOBAL, LATS_CUT_OFFS = lengths_for_SI(lines_SH, resol_thresh, conn_thresh)

            exp_SH_mod[i][0] = TOTAL_LENGTH
            exp_SH_mod[i][1] = TOTAL_LENGTH / (2 * np.pi * a_earth * np.cos(np.radians(np.mean(TOTAL_LATS)))) # normalise with mean lat
            exp_SH_mod[i][2] = np.mean(TOTAL_LATS)
            exp_SH_mod[i][3] = np.min(TOTAL_LATS)
            exp_SH_mod[i][4] = np.max(TOTAL_LATS)
            exp_SH_mod[i][5] = np.abs(np.min(TOTAL_LATS)) - np.abs(np.max(TOTAL_LATS))
            exp_SH_mod[i][6] = LENGTH_CIRCUMGLOBAL                                  # length circumglobal segment
            exp_SH_mod[i][7] = LENGTH_CIRCUMGLOBAL / (2 * np.pi * a_earth * np.cos(np.radians(np.mean(LATS_CIRCUMGLOBAL)))) # normalise with mean lat
            exp_SH_mod[i][8] = np.mean(LATS_CIRCUMGLOBAL)                           # mean lat longest segment
            exp_SH_mod[i][9] = np.min(LATS_CIRCUMGLOBAL)                            # max lat longest segment
            exp_SH_mod[i][10] = np.max(LATS_CIRCUMGLOBAL)                           # min lat longest segment
            exp_SH_mod[i][11] = np.abs(np.min(LATS_CIRCUMGLOBAL)) - np.abs(np.max(LATS_CIRCUMGLOBAL)) # meridional extent longest segment
            exp_SH_mod[i][12] = LENGTH_CUT_OFFS                                     # length cut-off segments
            exp_SH_mod[i][13] = LENGTH_CUT_OFFS / (2 * np.pi * a_earth * np.cos(np.radians(np.mean(LATS_CUT_OFFS)))) # normalise with mean lat
            exp_SH_mod[i][14] = np.mean(LATS_CUT_OFFS)                              # mean lat cut-off segments
            exp_SH_mod[i][15] = np.min(LATS_CUT_OFFS)                               # max lat cut-off segments
            exp_SH_mod[i][16] = np.max(LATS_CUT_OFFS)                               # min lat cut-off segments
            exp_SH_mod[i][17] = np.abs(np.min(LATS_CUT_OFFS)) - np.abs(np.max(LATS_CUT_OFFS))       # meridional extent cut-off segments
            exp_SH_mod[i][18] = NUMBER_CUT_OFFS                                     # number of cut-offs
            
            ## Writing the output
            file1.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,\n' % 
            (date_times[i].strftime("%Y%m%d%H"), GPH_mean_30_70[i,0], exp_NH_30_70[i,0], exp_NH_30_70[i,1], exp_NH_30_70[i,2], 
            exp_NH_30_70[i,3], exp_NH_30_70[i,4], exp_NH_30_70[i,5], exp_NH_30_70[i,6], exp_NH_30_70[i,7], exp_NH_30_70[i,8], 
            exp_NH_30_70[i,9], exp_NH_30_70[i,10], exp_NH_30_70[i,11], exp_NH_30_70[i,12], exp_NH_30_70[i,13], exp_NH_30_70[i,14],
            exp_NH_30_70[i,15], exp_NH_30_70[i,16], exp_NH_30_70[i,17], exp_NH_30_70[i,18])
            )

            file2.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,\n' % 
            (date_times[i].strftime("%Y%m%d%H"), GPH_mean_30_70[i,1], exp_SH_30_70[i,0], exp_SH_30_70[i,1], exp_SH_30_70[i,2],
            exp_SH_30_70[i,3], exp_SH_30_70[i,4], exp_SH_30_70[i,5], exp_SH_30_70[i,6], exp_SH_30_70[i,7], exp_SH_30_70[i,8],
            exp_SH_30_70[i,9], exp_SH_30_70[i,10], exp_SH_30_70[i,11], exp_SH_30_70[i,12], exp_SH_30_70[i,13], exp_SH_30_70[i,14],
            exp_SH_30_70[i,15], exp_SH_30_70[i,16], exp_SH_30_70[i,17], exp_SH_30_70[i,18])
            )

            file3.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,\n' % 
            (date_times[i].strftime("%Y%m%d%H"), GPH_mean_mod[i,0], exp_NH_mod[i,0], exp_NH_mod[i,1], exp_NH_mod[i,2],
            exp_NH_mod[i,3], exp_NH_mod[i,4], exp_NH_mod[i,5], exp_NH_mod[i,6], exp_NH_mod[i,7], exp_NH_mod[i,8],
            exp_NH_mod[i,9], exp_NH_mod[i,10], exp_NH_mod[i,11], exp_NH_mod[i,12], exp_NH_mod[i,13], exp_NH_mod[i,14],
            exp_NH_mod[i,15], exp_NH_mod[i,16], exp_NH_mod[i,17], exp_NH_mod[i,18])
            )

            file4.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,\n' % 
            (date_times[i].strftime("%Y%m%d%H"), GPH_mean_mod[i,1], exp_SH_mod[i,0], exp_SH_mod[i,1], exp_SH_mod[i,2],
            exp_SH_mod[i,3], exp_SH_mod[i,4], exp_SH_mod[i,5], exp_SH_mod[i,6], exp_SH_mod[i,7], exp_SH_mod[i,8],   
            exp_SH_mod[i,9], exp_SH_mod[i,10], exp_SH_mod[i,11], exp_SH_mod[i,12], exp_SH_mod[i,13], exp_SH_mod[i,14],
            exp_SH_mod[i,15], exp_SH_mod[i,16], exp_SH_mod[i,17], exp_SH_mod[i,18])
            )

            plt.close()
        del data, GPH, 
    
    file.close()
    file1.close()
    file2.close()
    file3.close()
    file4.close()
    
# close all windowd 
plt.close('all')
