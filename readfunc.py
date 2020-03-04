import h5py
import numpy as np

def readno2(InNO2FileName):
    
    data_l2 = h5py.File(InNO2FileName,'r')

    no2 = np.array(data_l2['/PRODUCT/nitrogendioxide_tropospheric_column'])
    qa  = np.array(data_l2['/PRODUCT/qa_value'])
    latno2 = np.array(data_l2['/PRODUCT/latitude'])
    lonno2 = np.array(data_l2['/PRODUCT/longitude'])

    qa = np.squeeze(qa); no2 = np.squeeze(no2); latno2 = np.squeeze(latno2); lonno2 = np.squeeze(lonno2)

    return qa, no2, latno2, lonno2


def readmet(InFileName):

    # h5dump -n [filename]
    data_l2 = h5py.File(InFileName,'r')

    uwind  = np.array(data_l2['/Meteorology/wind_u_profile_met'])
    vwind  = np.array(data_l2['/Meteorology/wind_v_profile_met'])
    metlat = np.array(data_l2['/SoundingGeometry/sounding_latitude'])
    metlon = np.array(data_l2['/SoundingGeometry/sounding_longitude'])
    metsid = np.array(data_l2['/SoundingGeometry/sounding_id'])
    spress = np.array(data_l2['/Meteorology/surface_pressure_met'])
    ak     = np.array(data_l2['/MeteorologyDiagnostics/ak'])
    bk     = np.array(data_l2['/MeteorologyDiagnostics/bk'])

    return uwind, vwind, metlat, metlon

