'''
This is a Python script to download publicly available materials required for generating figures.
usage: python fetchMaterials.py
author: Ryo Okuwaki (rokuwaki@geol.tsukuba.ac.jp) 2021-06-27
'''

import requests
import sys
import os
import pandas as pd
import io
import git

def main():

    outdir = '../materials/work'
    os.makedirs(outdir, exist_ok=True)

    # SPUD GCMT
    URL = 'http://ds.iris.edu/spudservice/momenttensor/bundle/quakeml?evtminlat=-40.0&evtmaxlat=-35.0&evtminlon=178.0&evtmaxlon=-178.0&evtmindepth=0.0&evtmaxdepth=700.0&evtminmag=0.0&evtmaxmag=10.0&evtstartdate=1950-01-01T00:00:00&evtenddate=2021-05-29T21:04:00&minstrike=0&maxstrike=360'
    response = requests.get(URL)
    with open(os.path.join(outdir, 'SPUD_QUAKEML_bundle.xml'), 'wb') as file:
        file.write(response.content)

    # USGS W-phase
    URL = 'https://earthquake.usgs.gov/archive/product/moment-tensor/us_7000dffl_mww/us/1620945817040/quakeml.xml'
    response = requests.get(URL)
    with open(os.path.join(outdir, 'quakeml_USGS_Wphase.xml'), 'wb') as file:
        file.write(response.content)
        
    # USGS Moment tensors
    URL = 'https://earthquake.usgs.gov/fdsnws/event/1/query.quakeml?starttime=1900-08-24%2000:00:00&endtime=2021-08-31%2023:59:59&maxlatitude=-36&minlatitude=-40&maxlongitude=182&minlongitude=177&minmagnitude=5&includeallorigins=true&includeallmagnitudes=true&orderby=time&producttype=moment-tensor'
    response = requests.get(URL)
    with open(os.path.join(outdir, 'quakeml_USGS_MomentTensor.xml'), 'wb') as file:
        file.write(response.content)
    
    # GeoNet CMT solution
    url = 'https://raw.githubusercontent.com/GeoNet/data/main/moment-tensor/GeoNet_CMT_solutions.csv'
    dfGeoNetCMTData = pd.read_csv(url)
    dfGeoNetCMTData.Date = pd.to_datetime(dfGeoNetCMTData.Date, format='%Y%m%d%H%M%S')
    dfGeoNetCMTData.Longitude[dfGeoNetCMTData.Longitude < 0] += 360
    dfGeoNetCMTData.to_csv(os.path.join(outdir, 'GeoNet_CMT_solutions.csv'))

    
    # NZ station information
    url = 'http://beta-service.geonet.org.nz/fdsnws/station/1/query?&format=text&level=station&network=NZ&starttime=2021-03-04T13:21:00'
    urlData = requests.get(url).content
    df = pd.read_csv(io.StringIO(urlData.decode('utf-8')), skipinitialspace=True, delimiter='|')
    df.to_csv(os.path.join(outdir, 'NZstation_info.txt'), sep='|', index=False)

    # Plate boundaries; Bird, 2003, doi:10.1029/2001GC000252
    try:
        git.Git(outdir).clone('https://github.com/fraxen/tectonicplates.git')
    except git.GitCommandError as e:
        pass

    # GeoNet seismicity
    url = 'https://quakesearch.geonet.org.nz/csv?bbox=178,-40,-179,-36&startdate=2011-03-04T13:26:35&enddate=2021-04-11T13:27:35'
    urlData = requests.get(url).content
    df = pd.read_csv(io.StringIO(urlData.decode('utf-8')), skipinitialspace=True)
    df.to_csv(os.path.join(outdir, 'geonetBackgroundSeismicity.csv'))

    # USGS seismicity
    url = 'https://earthquake.usgs.gov/fdsnws/event/1/query.csv?starttime=2011-03-04%2013:26:35&endtime=2021-04-11%2013:27:35&maxlatitude=-36&minlatitude=-40&maxlongitude=181&minlongitude=178&minmagnitude=0&orderby=time'
    urlData = requests.get(url).content
    df = pd.read_csv(io.StringIO(urlData.decode('utf-8')), skipinitialspace=True)
    df.to_csv(os.path.join(outdir, 'USGSBackgroundSeismicity.csv'))

    # Slab2 model
    url = 'https://www.sciencebase.gov/catalog/file/get/5aa318e1e4b0b1c392ea3f10?f=__disk__c6%2F8a%2F17%2Fc68a179673fdd9c737abf523f0c7c2675471955c'
    urlData = requests.get(url).content
    df = pd.read_csv(io.StringIO(urlData.decode('utf-8')), skipinitialspace=True)
    df.to_csv(os.path.join(outdir, 'ker_slab2_clp_02.24.18.csv'))

    URL = 'https://www.sciencebase.gov/catalog/file/get/5aa318e1e4b0b1c392ea3f10?f=__disk__01%2Fb6%2Ff7%2F01b6f7842985489f92723e43140d07039c89f049'
    response = requests.get(URL)
    with open(os.path.join(outdir, 'ker_slab2_dep_02.24.18.grd'), 'wb') as file:
        file.write(response.content)


if __name__ == '__main__':
    main()
