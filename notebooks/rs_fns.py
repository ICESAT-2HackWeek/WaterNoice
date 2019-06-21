import getpass
import h5py
import numpy as np
import pyproj
from astropy.time import Time
import datetime as dt
import rasterio as rio
from rasterio.plot import show
import gdal
import os


'''#----------------------------------
            General functions
'''#----------------------------------
def write_file(header, data, name, path):
#write tab-delimited data file
#header should be a list of strings that is names of the tab-delimited columns
#name is output file name
#path is where you want to write the file
#data is an array that contains all of the variables - all variables must be the same length

    f_out = open(path+name,'w')
    for i in range(0,100,len(header)):
        f_out.write(header[i]+'\t')

    f_out.write('\n')

#     for i in range(len(header[0])):
#         for j in range(len(header)):
#             f_out.write(str(data[j][i])+'\t')
#         f_out.write('\n')

    f_out.close()


def list_files_local(path):
    """ Get file list form local folder. """
    from glob import glob
    return glob(path)


def gps2dyr(time):
    """ Converte GPS time to decimal years. """
    return Time(time, format='gps').decimalyear

def transform_coord(proj1, proj2, x, y):
    """
    Transform coordinates from proj1 to proj2 (EPSG num).

    Example EPSG projs:
        Geodetic (lon/lat): 4326
        Polar Stereo AnIS (x/y): 3031
        Polar Stereo GrIS (x/y): 3413
    """
    # Set full EPSG projection strings
    proj1 = pyproj.Proj("+init=EPSG:"+str(proj1))
    proj2 = pyproj.Proj("+init=EPSG:"+str(proj2))
    return pyproj.transform(proj1, proj2, x, y)  # convert

def read_h5(fname, vnames=[]):
    """ Simple HDF5 reader. """
    with h5py.File(fname, 'r') as f:
        return [f[v][:] for v in vnames]

'''#----------------------------------
            ICESat-2
'''#----------------------------------
   
    # File name format: ATL06_[yyyymmdd][hhmmss]_[RGTccss]_[vvv_rr].h5

#NOTE: Need to simplify this function
def time_from_fname(fname):
    """ IS2 fname -> datatime object. """
    t = fname.split('_')[1]
    y, m , d, h, mn, s = t[:4], t[4:6], t[6:8], t[8:10], t[10:12], t[12:14]
    time = dt.datetime(int(y), int(m), int(d), int(h), int(mn), int(s))
    return time


def segment_from_fname(fname):
    """ IS2 fname -> segment number. """
    s = fname.split('_')[2]
    return int(s[-2:])


def select_files(files, segments=[10,11,12], t1=(2019,1,1), t2=(2019,2,1)):
    t1 = dt.datetime(*t1)
    t2 = dt.datetime(*t2)
    files_out = []
    for f in files:
        fname = os.path.basename(f)
        time = time_from_fname(fname)
        segment = segment_from_fname(fname)
        if t1 <= time <= t2 and segment in segments:
            files_out.append(f)
    return files_out

def track_type(time, lat, tmax=1):
    """
    Separate tracks into ascending and descending.
    
    Defines tracks as segments with time breaks > tmax,
    and tests whether lat increases or decreases w/time.
    """
    tracks = np.zeros(lat.shape)  # generate track segment
    tracks[0:np.argmax(np.abs(lat))] = 1  # set values for segment
    i_asc = np.zeros(tracks.shape, dtype=bool)  # output index array

    # Loop trough individual secments
    for track in np.unique(tracks):
    
        i_track, = np.where(track == tracks)  # get all pts from seg
    
        if len(i_track) < 2: continue
    
        # Test if lat increases (asc) or decreases (des) w/time
        i_min = time[i_track].argmin()
        i_max = time[i_track].argmax()
        lat_diff = lat[i_track][i_max] - lat[i_track][i_min]
    
        # Determine track type
        if lat_diff > 0:  i_asc[i_track] = True
    
    return i_asc, np.invert(i_asc)  # index vectors

def read_atl03(fname, bbox=None):
    print(fname)
    """ 
    Read 1 ATL06 file and output 6 reduced files. 
    
    Extract variables of interest and separate the ATL06 file 
    into each beam (ground track) and ascending/descending orbits.
    """
    f = h5py.File(fname,'r')
    # Each beam is a group
    group = ['/gt1l', '/gt1r', '/gt2l', '/gt2r', '/gt3l', '/gt3r']
    # Loop trough beams
    for k,g in enumerate(group):
    
        
        # 1) Read in data for a single beam #
        
    
        # Load variables into memory (more can be added!)
        with h5py.File(fname, 'r') as fi:
            if g+'/heights' in fi.keys():
                lat = fi[g+'/heights/lat_ph'][:]
                lon = fi[g+'/heights/lon_ph'][:]
                h = fi[g+'/heights/h_ph'][:]
                conf = fi[g+'/heights/signal_conf_ph']
                bkg = fi[g+'/bckgrd_atlas/bckgrd_counts'][:]
                #dist = fi[g+'/heights/dist_ph_along'][:]
                #seg = fi[g+'/geolocation/segment_dist_x'][:]

                land_ice_class = conf[:,3]
          
        
        # 2) Filter data #
        
                mask = (land_ice_class == 4) & (np.abs(h) < 10e3)
                lat,lon,h = lat[mask],lon[mask],h[mask]
        
        
        # 3) Save selected data #
        
        
        # Define output file name
                ofile = fname.replace('.h5', '_'+g[1:]+'.h5')
                
        # Save variables
                with h5py.File(ofile, 'w') as f:
                    f['lon'] = lon
                    f['lat'] = lat
                    f['h_elv'] = h
                    f['bkg_ct'] = bkg
                    #f['x_atc'] = x_atc
            
                    print('out ->', ofile)


def read_atl06(fname, bbox=None):
    """ 
    Read 1 ATL06 file and output 6 reduced files. 
    
    Extract variables of interest and separate the ATL06 file 
    into each beam (ground track) and ascending/descending orbits.
    """

    # Each beam is a group
    group = ['/gt1l', '/gt1r', '/gt2l', '/gt2r', '/gt3l', '/gt3r']

    # Loop trough beams
    for k,g in enumerate(group):
    
        # 1) Read in data for a single beam #
    
        # Load variables into memory (more can be added!)
        with h5py.File(fname, 'r') as fi:
            lat = fi[g+'/land_ice_segments/latitude'][:]
            lon = fi[g+'/land_ice_segments/longitude'][:]
            h_li = fi[g+'/land_ice_segments/h_li'][:]
            s_li = fi[g+'/land_ice_segments/h_li_sigma'][:]
            t_dt = fi[g+'/land_ice_segments/delta_time'][:]
            q_flag = fi[g+'/land_ice_segments/atl06_quality_summary'][:]
            s_fg = fi[g+'/land_ice_segments/fit_statistics/signal_selection_source'][:]
            snr = fi[g+'/land_ice_segments/fit_statistics/snr_significance'][:]
            h_rb = fi[g+'/land_ice_segments/fit_statistics/h_robust_sprd'][:]
            dac = fi[g+'/land_ice_segments/geophysical/dac'][:]
            f_sn = fi[g+'/land_ice_segments/geophysical/bsnow_conf'][:]
            dh_fit_dx = fi[g+'/land_ice_segments/fit_statistics/dh_fit_dx'][:]
            tide_earth = fi[g+'/land_ice_segments/geophysical/tide_earth'][:]
            tide_load = fi[g+'/land_ice_segments/geophysical/tide_load'][:]
            tide_ocean = fi[g+'/land_ice_segments/geophysical/tide_ocean'][:]
            tide_pole = fi[g+'/land_ice_segments/geophysical/tide_pole'][:]
            t_ref = fi['/ancillary_data/atlas_sdp_gps_epoch'][:]
            rgt = fi['/orbit_info/rgt'][:] * np.ones(len(lat))
            orb = np.full_like(h_li, k)

        # 2) Filter data according region and quality #
        
        # Select a region of interest
        if bbox:
            lonmin, lonmax, latmin, latmax = bbox
            bbox_mask = (lon >= lonmin) & (lon <= lonmax) & \
                        (lat >= latmin) & (lat <= latmax)
        else:
            bbox_mask = np.ones_like(lat, dtype=bool)  # get all
            
        # Only keep good data, and data inside bbox
        mask = (q_flag == 0) & (np.abs(h_li) < 10e3) & (bbox_mask == 1)
        
        # Update variables
        lat, lon, h_li, s_li, t_dt, h_rb, s_fg, snr, q_flag, f_sn, \
            tide_earth, tide_load, tide_ocean, tide_pole, dac, rgt, orb = \
                lat[mask], lon[mask], h_li[mask], s_li[mask], t_dt[mask], \
                h_rb[mask], s_fg[mask], snr[mask], q_flag[mask], f_sn[mask], \
                tide_earth[mask], tide_load[mask], tide_ocean[mask], \
                tide_pole[mask], dac[mask], rgt[mask], orb[mask]

        # Test for no data
        if len(h_li) == 0: continue

        # 3) Convert time and separate tracks #
        
        # Time in GPS seconds (secs sinde 1980...)
        t_gps = t_ref + t_dt

        # Time in decimal years
        t_year = gps2dyr(t_gps)

        # Determine orbit type
        i_asc, i_des = track_type(t_year, lat)
        
        # 4) Save selected data #
        
        # Define output file name
        ofile = fname.replace('.h5', '_'+g[1:]+'.h5')
                
        # Save variables
        with h5py.File(ofile, 'w') as f:
            f['orbit'] = orb
            f['lon'] = lon
            f['lat'] = lat
            f['h_elv'] = h_li
            f['t_year'] = t_year
            f['t_sec'] = t_gps
            f['s_elv'] = s_li
            f['h_rb'] = h_rb
            f['s_fg'] = s_fg
            f['snr'] = snr
            f['q_flg'] = q_flag
            f['f_sn'] = f_sn
            f['tide_load'] = tide_load
            f['tide_ocean'] = tide_ocean
            f['tide_pole'] = tide_pole
            f['tide_earth'] = tide_earth
            f['dac'] = dac
            f['rgt'] = rgt
            f['trk_type'] = i_asc

            print('out ->', ofile)
            
            
'''#----------------------------------
                Landsat-8
'''#----------------------------------
            
def normalize(array): #Normalize bands into 0.0 - 1.0 scale
    array_min, array_max = array.min(), array.max()
    return ((array - array_min)/(array_max - array_min))

def stack_landsat(files): #input a list of file names in order: rgb
    red = np.squeeze(rio.open(files[0]).read())
    green = np.squeeze(rio.open(files[1]).read())
    blue = np.squeeze(rio.open(files[2]).read())
    
    blue = normalize(blue)
    red = normalize(red)
    green = normalize(green)
    
    rgb = np.dstack((red, green, blue))
    rgb = rgb.astype(float)
    
    #get stereo polargraphic coords
    Raster = gdal.Open(files[0])
    width = Raster.RasterXSize
    height = Raster.RasterYSize
    gt = Raster.GetGeoTransform()
    array = Raster.ReadAsArray()

    # Pixel numbers
    x = np.arange(0, width)
    y = np.arange(0, height)

    # Grid Cell Coordinates of upper left corner in EPSG:3031 UTM. 
    X = gt[0] + x * gt[1] 
    Y = gt[3] + y * gt[5]
    
    return X,Y,rgb