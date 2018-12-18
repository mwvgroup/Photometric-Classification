#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module formats SDSS data for use with snpy and fits the SDSS light
curves"""


def write_snoopy_sdss_file(out_path, **kwargs):
    """Creates a snoopy input file <out_path>/<name>.txt

    Args:
        out_path         (str): Where to write the input file
        name             (str): The name of the target
        redshift       (float): The redshift of the target
        ra             (float): The ra of the target
        dec            (float): The dec of the target
        mjd            (float): The date of each observation
        <u, g, r, i, z> (dict): Dict with date ('mjd'), magnitude ('mag'),
                                  and magnitude error ('err')
    """

    file_text = ' '.join([kwargs["name"],
                          kwargs["redshift"],
                          kwargs["ra"],
                          kwargs["dec"]])

    for band in 'ugriz':
        filter_name = 'sdss_' + band
        if band not in kwargs:
            continue

        file_text += 'filter {}\n'.format(filter_name)
        band_data = kwargs[band]
        for mjd, mag, error in \
                zip(band_data['mjd'], band_data['mag'], band_data['err']):
            file_text += '{} {} {}\n'.format(mjd, mag, error)

    with open(out_path, 'w') as ofile:
        ofile.write(file_text)
