'''
This is a Python script to process bathymetry data using `gmt` commands.
usage: python prepBathy.py directory_to_strore_bathymetry_ncdata
author: Ryo Okuwaki (rokuwaki@geol.tsukuba.ac.jp) 2021-05-27
dependency: Generic Mappring Tools, Version 6.1.0 or higher
'''

import os
import sys
import numpy as np
import subprocess
import zipfile
import requests
import io

def main():

    args = sys.argv
    bathy_work_dir = '../materials/work/bathymetry'
    bathy_nc_dir = args[1]
    os.makedirs(bathy_work_dir, exist_ok=True)

    # Download bathymetry data
    print('Requesting data: NZ 250m ESRI ASCII grid [ZIP 563 MB] (this may take a while) ...')
    url = 'https://niwa.co.nz/static/bathymetry/NZBathy_DTM_2016_ascii_grid.zip'
    r = requests.get(url)

    print('Extracting zip: NZ 250m ESRI ASCII grid [ZIP 563 MB] (this may take a while) ...')
    z = zipfile.ZipFile(io.BytesIO(r.content))
    for f in z.filelist:
        if f.file_size == 0:
            print('Caution!! ', f.filename, ' has no data.')
        else:
            z.extract(f.filename, bathy_nc_dir)

    # trasform txt file into netCDF grid file
    command = ['gdalwarp','-s_srs','EPSG:3994','-t_srs','EPSG:4326','--config',
               'CENTER_LONG','180','-overwrite',
               bathy_nc_dir+'/nzbathymetry_2016_ascii-grid.txt',
               bathy_nc_dir+'/bathy.nc']
    process = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout,stderr)

    # grdcut
    command = ['gmt','grdcut',bathy_nc_dir+'/bathy.nc',
               '-R176/183/-40/-35','-G'+bathy_nc_dir+'/bathy_cut.nc']
    process = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout,stderr)


    # Base line along trench
    command = ['gmt','project','-C179.85/-37.466','-A200','-G0.5','-L0/50','-Q','>',
               bathy_work_dir+'/tmp.txt']
    process = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout,stderr)

    # cross section
    data = np.loadtxt(bathy_work_dir+'/tmp.txt')
    lons, lats = data[:,0], data[:,1]
    i = 1
    for lon, lat in zip(lons, lats):
        print(lon, lat)
        command0 = ['gmt','project','-C'+str(lon)+'/'+str(lat),'-A110','-G0.5',
                   '-L-70/0','-Q','>',bathy_work_dir+'/tmp']
        command1 = ['gmt','grdtrack', bathy_work_dir+'/tmp',
                   '-G'+bathy_nc_dir+'/bathy_cut.nc','>',bathy_work_dir+'/Xsection_'+str(i)+'.txt']
        for command in [command0,command1]:
            process = subprocess.Popen(command,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            print(stdout,stderr)

        i += 1

if __name__ == '__main__':
    main()
