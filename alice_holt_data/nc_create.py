import numpy as np
import netCDF4 as nC
import ah_nc_utils as utils
import matplotlib.mlab as ml
import time as tt
import pandas as pd
import datetime as dt


def create_netcdf_dataset(path):
    """Creates netcdf dataset for half hourly flux measurements at Alice Holt
    """
    dataset = nC.Dataset(path+'ah_data_half_hourly_compress.nc', 'w', format='NETCDF4_CLASSIC')
    time = dataset.createDimension('time', 298032)
    lat = dataset.createDimension('lat', 1)
    lon = dataset.createDimension('lon', 1)
    times = dataset.createVariable('time', np.float64, ('time',), zlib=True, complevel=9)
    latitudes = dataset.createVariable('latitude', np.float32, ('lat',), zlib=True, least_significant_digit=4
                                       , complevel=9)
    latitudes[0] = 51.153526
    longitudes = dataset.createVariable('longitude', np.float32, ('lon',), zlib=True, least_significant_digit=4
                                        , complevel=9)
    longitudes[0] = -0.85835201

    dataset.description = 'Alice Holt Straits Inclosure half hourly data'
    dataset.history = 'Created ' + tt.ctime(tt.time())
    dataset.source = 'Ewan Pinnington, University of Reading. email: ewan.pinnington@gmail.com'

    latitudes.units = 'degrees north'
    longitudes.units = 'degrees east'
    times.units = 'minutes since 1970-01-01 00:00:00.0'
    times.calendar = 'gregorian'
    is_day = dataset.createVariable('is_day', 'i1', ('time', 'lat', 'lon'), zlib=True, least_significant_digit=1
                                      , complevel=9)
    is_day.description = 'displays 1 or 0 if corresponding time is day (1) or night (0)'

    air_temp = dataset.createVariable('air_temp', 'f4', ('time', 'lat', 'lon'), zlib=True)
    air_temp.units = 'degC'
    air_temp.description = 'air temperature at 27m'
    air_temp.standard_name = 'air_temperature'
    soil_temp = dataset.createVariable('soil_temp', 'f4', ('time', 'lat', 'lon'), zlib=True)
    soil_temp.units = 'degC'
    soil_temp.description = 'soil temperature at 3cm'
    rg = dataset.createVariable('rg', 'f4', ('time', 'lat', 'lon'), zlib=True)
    rg.units = 'W m-2'
    rg.standard_name = 'surface_downwelling_shortwave_flux_in_air'
    co2_flux = dataset.createVariable('co2_flux', 'f4', ('time', 'lat', 'lon'), zlib=True)
    co2_flux.units = 'umol m-2 s-1'
    co2_flux.standard_name = 'surface_upward_mole_flux_of_carbon_dioxide'
    co2_flux.description = 'unprocessed Alice Holt flux tower record'
    qc_co2_flux = dataset.createVariable('qc_co2_flux', 'i1', ('time', 'lat', 'lon'), zlib=True,
                                         least_significant_digit=1, complevel=9)
    qc_co2_flux.units = 'none'
    qc_co2_flux.description ='quality control flag for half hourly co2 flux observations, 0 - good, 2 - bad'
    u_star = dataset.createVariable('u_star', 'f4', ('time', 'lat', 'lon'), zlib=True)
    u_star.units = 'm s-1'
    wind_dir = dataset.createVariable('wind_dir', 'f4', ('time', 'lat', 'lon'), zlib=True)
    wind_dir.units = 'degree'
    wind_dir.standard_name = 'wind_from_direction'
    wind_dir.description = 'wind direction degrees from north'
    foot_print = dataset.createVariable('foot_print', 'f4', ('time', 'lat', 'lon'), zlib=True)
    foot_print.units = 'm'
    foot_print.description = 'distance from tower where 90% of co2 flux is measured'

    start_date = dt.datetime(1999, 1, 1, 0, 0)
    end_date = dt.datetime(2015, 12, 31, 23, 30)
    time_lst = utils.create_date_list(start_date, end_date, del_t='half_hour')
    times[:] = nC.date2num(time_lst, times.units, times.calendar)
    return dataset


def create_netcdf_dataset_daily(path):
    """Creates netcdf dataset for half hourly flux measurements at Alice Holt
    """
    dataset = nC.Dataset(path+'ah_data_daily.nc', 'w', format='NETCDF4_CLASSIC')
    time = dataset.createDimension('time', 6209)
    lat = dataset.createDimension('lat', 1)
    lon = dataset.createDimension('lon', 1)
    times = dataset.createVariable('time', np.float64, ('time',), zlib=True)
    latitudes = dataset.createVariable('latitude', np.float32, ('lat',), zlib=True, least_significant_digit=4)
    latitudes[0] = 51.153526
    longitudes = dataset.createVariable('longitude', np.float32, ('lon',), zlib=True, least_significant_digit=4)
    longitudes[0] = -0.85835201

    dataset.description = 'Alice Holt Straits Inclosure daily averaged data'
    dataset.history = 'Created ' + tt.ctime(tt.time())
    dataset.source = 'Ewan Pinnington, University of Reading. email: ewan.pinnington@gmail.com'

    latitudes.units = 'degrees north'
    longitudes.units = 'degrees east'
    times.units = 'minutes since 1970-01-01 00:00:00.0'
    times.calendar = 'gregorian'

    # Daily global radiation
    rg = dataset.createVariable('rg', 'f4', ('time', 'lat', 'lon'), zlib=True)
    rg.units = 'MJ m-2 day-1'
    rg.standard_name = 'surface_downwelling_shortwave_flux_in_air'

    # Daily temperatures
    daily_max_temp = dataset.createVariable('daily_max_temp', 'f4', ('time', 'lat', 'lon'), zlib=True)
    daily_max_temp.units = 'degC'
    daily_max_temp.description = 'maximum daily air temperature at 27m'
    daily_max_temp.standard_name = 'air_temperature'
    daily_min_temp = dataset.createVariable('daily_min_temp', 'f4', ('time', 'lat', 'lon'), zlib=True)
    daily_min_temp.units = 'degC'
    daily_min_temp.description = 'minimum daily air temperature at 27m'
    daily_min_temp.standard_name = 'air_temperature'
    daily_mean_temp = dataset.createVariable('daily_mean_temp', 'f4', ('time', 'lat', 'lon'), zlib=True)
    daily_mean_temp.units = 'degC'
    daily_mean_temp.description = 'mean daily air temperature at 27m'
    daily_mean_temp.standard_name = 'air_temperature'
    mean_temp_day = dataset.createVariable('mean_temp_day', 'f4', ('time', 'lat', 'lon'), zlib=True)
    mean_temp_day.units = 'degC'
    mean_temp_day.description = 'mean daytime air temperature at 27m'
    mean_temp_day.standard_name = 'air_temperature'
    mean_temp_night = dataset.createVariable('mean_temp_night', 'f4', ('time', 'lat', 'lon'), zlib=True)
    mean_temp_night.units = 'degC'
    mean_temp_night.description = 'mean nighttime air temperature at 27m'
    mean_temp_night.standard_name = 'air_temperature'

    daily_mean_soil_temp = dataset.createVariable('daily_mean_soil_temp', 'f4', ('time', 'lat', 'lon'), zlib=True)
    daily_mean_soil_temp.units = 'degC'
    daily_mean_soil_temp.description = 'mean daily soil temperature at 3cm'
    mean_soil_temp_day = dataset.createVariable('mean_soil_temp_day', 'f4', ('time', 'lat', 'lon'), zlib=True)
    mean_soil_temp_day.units = 'degC'
    mean_soil_temp_day.description = 'mean daytime soil temperature at 3cm'
    mean_soil_temp_night = dataset.createVariable('mean_soil_temp_night', 'f4', ('time', 'lat', 'lon'), zlib=True)
    mean_soil_temp_night.units = 'degC'
    mean_soil_temp_night.description = 'mean nighttime soil temperature at 3cm'

    # Day/night lengths
    day_length = dataset.createVariable('day_length', 'f4', ('time', 'lat', 'lon'), zlib=True,
                                        least_significant_digit=2)
    day_length.units = 'hours'
    day_length.description = 'day length in hours'
    night_length = dataset.createVariable('night_length', 'f4', ('time', 'lat', 'lon'), zlib=True,
                                          least_significant_digit=2)
    night_length.units = 'hours'
    night_length.description = 'night length in hours'
    doy = dataset.createVariable('doy', 'f4', ('time', 'lat', 'lon'), zlib=True,
                                 least_significant_digit=1)
    doy.description = 'day of year'

    # net ecosystem exchange for whole site
    nee = dataset.createVariable('nee', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee.units = 'g C m-2 day-1'
    nee.standard_name = 'surface_upward_mole_flux_of_carbon_dioxide'
    nee.description = 'processed total daily net ecosystem exchange from Alice Holt flux site'
    nee_std = dataset.createVariable('nee_std', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee_std.units = 'g C m-2 day-1'
    nee_std.description = 'standard deviation for total daily net ecosystem exchange'
    nee_day = dataset.createVariable('nee_day', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee_day.units = 'g C m-2 day-1'
    nee_day.standard_name = 'surface_upward_mole_flux_of_carbon_dioxide'
    nee_day.description = 'processed total daytime net ecosystem exchange from Alice Holt flux site'
    nee_day_std = dataset.createVariable('nee_day_std', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee_day_std.units = 'g C m-2 day-1'
    nee_day_std.description = 'standard deviation for total daytime net ecosystem exchange'
    nee_night = dataset.createVariable('nee_night', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee_night.units = 'g C m-2 day-1'
    nee_night.standard_name = 'surface_upward_mole_flux_of_carbon_dioxide'
    nee_night.description = 'processed total nighttime net ecosystem exchange from Alice Holt flux site'
    nee_night_std = dataset.createVariable('nee_night_std', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee_night_std.units = 'g C m-2 day-1'
    nee_night_std.description = 'standard deviation for total nighttime net ecosystem exchange'
    nee_origin = dataset.createVariable('nee_origin', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee_origin.description = 'origin of flux measurement (E, W, undetermined)'
    nee_day_origin = dataset.createVariable('nee_day_origin', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee_day_origin.description = 'origin of flux measurement (E, W, undetermined)'
    nee_night_origin = dataset.createVariable('nee_night_origin', 'f4', ('time', 'lat', 'lon'), zlib=True)
    nee_night_origin.description = 'origin of flux measurement (E, W, undetermined)'

    # leaf area index
    lai = dataset.createVariable('lai', 'f4', ('time', 'lat', 'lon'), zlib=True)
    lai.standard_name = 'leaf_area_index'
    lai.description = 'average lai for whole site'
    lai_east = dataset.createVariable('lai_east', 'f4', ('time', 'lat', 'lon'), zlib=True)
    lai_east.standard_name = 'leaf_area_index'
    lai_east.description = 'average lai for site east of flux tower'
    lai_west = dataset.createVariable('lai_west', 'f4', ('time', 'lat', 'lon'), zlib=True)
    lai_west.standard_name = 'leaf_area_index'
    lai_west.description = 'average lai for site west of flux tower'
    lma = dataset.createVariable('lma', 'f4', ('time', 'lat', 'lon'), zlib=True)
    lma.description = 'leaf mass area for Alice Holt'
    lma.units = 'g C m-2'

    # woody biomass
    c_woo = dataset.createVariable('c_woo', 'f4', ('time', 'lat', 'lon'), zlib=True)
    c_woo.standard_name = 'wood_carbon_content'
    c_woo.units = 'g C m-2'
    c_woo.description = 'average wood carbon stock for whole site'
    c_woo_east = dataset.createVariable('c_woo_east', 'f4', ('time', 'lat', 'lon'), zlib=True)
    c_woo_east.standard_name = 'wood_carbon_content'
    c_woo_east.description = 'average wood carbon stock for site east of flux tower'
    c_woo_west = dataset.createVariable('c_woo_west', 'f4', ('time', 'lat', 'lon'), zlib=True)
    c_woo_west.standard_name = 'wood_carbon_content'
    c_woo_west.description = 'average wood carbon stock for site west of flux tower'

    # fine root biomass
    c_roo = dataset.createVariable('c_roo', 'f4', ('time', 'lat', 'lon'), zlib=True)
    c_roo.units = 'g C m-2'
    c_roo.description = 'average fine root carbon stock for whole site'
    c_roo_east = dataset.createVariable('c_roo_east', 'f4', ('time', 'lat', 'lon'), zlib=True)
    c_roo_east.description = 'average fine root carbon stock for site east of flux tower'
    c_roo_west = dataset.createVariable('c_roo_west', 'f4', ('time', 'lat', 'lon'), zlib=True)
    c_roo_west.description = 'average fine root carbon stock for site west of flux tower'

    # leaf on/leaf off
    d_on = dataset.createVariable('d_onset', 'f4', ('time', 'lat', 'lon'), zlib=True)
    d_on.units = 'day of year'
    d_on.description = 'day of labile release to leaves'
    d_off = dataset.createVariable('d_fall', 'f4', ('time', 'lat', 'lon'), zlib=True)
    d_off.units = 'day of year'
    d_off.description = 'day of leaf fall'
    cronset = dataset.createVariable('cronset', 'f4', ('time', 'lat', 'lon'), zlib=True)
    cronset.units = 'days'
    cronset.description = 'length of labile release to leaves in days'

    start_date = dt.datetime(1999, 1, 1)
    end_date = dt.datetime(2015, 12, 31)
    time_lst = utils.create_date_list(start_date, end_date, del_t='day')
    times[:] = nC.date2num(time_lst, times.units, times.calendar)
    return dataset


def open_netcdf(filename):
    """Opens a netCDF file
    """
    return nC.Dataset(filename, 'a')


def open_xls_sheet(filename, sheet_name):
    return pd.read_excel(filename, sheet_name)


def merge_csv_files(direct):
    f_out = open("ah_data_half_hourly.csv", "a")
    # first file:
    for line in open(direct+"flux_met_ah_1999.csv"):
        f_out.write(line)
    # now the rest:
    for num in range(2000, 2016):
        f = open(direct+"flux_met_ah_"+str(num)+".csv")
        f.next()  # skip the header
        for line in f:
            f_out.write(line)
        f.close()  # not really needed
    f_out.close()


def add_data2nc(nc_data, pd_df, data_title, nc_title, date_l='None', date_col='date_combined', sel='exact'):
    """ Adds data to a netCDF file
    :param nc_data: netCDF data set object
    :param pd_df: pandas data frame object
    :param data_title: title column for data to add as str
    :param nc_title: title of nc variable to add it to as str
    :param date_col: title of date column as str
    :return: nothing
    """
    var = nc_data.variables[nc_title]
    times = nc_data.variables['time']
    for x in xrange(len(pd_df[date_col])):
        try:
            if date_l != 'None':
                tm = round_time_nearest_10min(date_l[x])
            else:
                tm = round_time_nearest_10min(pd_df[date_col][x])
            # Find datetime index
            idx = nC.date2index(tm, times, calendar=times.calendar, select=sel)
        except TypeError:
            print x
            break
        except ValueError:
            print x
            break
        var[idx, 0, 0] = pd_df[data_title][x]


def add_excel_ah_obs(xls_file, nc_file, start_yr=1999, end_yr=2016, sel='exact'):
    years = np.arange(start_yr, end_yr)
    nc_data = open_netcdf(nc_file)
    nc_vars = ['air_temp', 'soil_temp', 'rg', 'co2_flux', 'qc_co2_flux', 'u_star', 'wind_dir', 'foot_print']
    for yr in years:
        print yr
        pd_df = open_xls_sheet(xls_file, str(yr))
        start_date = round_time_nearest_10min(pd_df['date_combined'][0])
        end_date = round_time_nearest_10min(pd_df['date_combined'][len(pd_df)-1])
        date_list = utils.create_date_list(start_date, end_date, del_t='half_hour')
        for var_title in nc_vars:
            print var_title
            add_data2nc(nc_data, pd_df, var_title, var_title, date_l=date_list, sel=sel)
    nc_data.close()
    return 'net_cdf file updated!'


def round_time_nearest_10min(datet):
    tm = datet  # datetime for var
    # Round datetime to nearest 10min mark
    discard = dt.timedelta(minutes=tm.minute % 10,
                     seconds=tm.second,
                     microseconds=tm.microsecond)
    tm -= discard
    if discard >= dt.timedelta(minutes=5):
        tm += dt.timedelta(minutes=10)
    return tm


def add_csv2nc(csv_name, nc_name):
    csv_data = ml.csv2rec(csv_name)
    nc_data = open_netcdf(nc_name)
    nc_vars = ['air_temp', 'soil_temp', 'rg', 'co2_flux', 'qc_co2_flux', 'u_star', 'wind_dir', 'foot_print']
    for var_title in nc_vars:
        nc_var = nc_data.variables[var_title]
        csv_var = csv_data[var_title]
        if len(csv_var) != len(nc_var):
            raise ValueError('Cannot project data of different shapes together')
        else:
            nc_var[:] = csv_var[:]
    nc_data.close()
    return 'All updated'


