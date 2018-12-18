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
        <u, g, r, i, z> (dict): Dict with the date ('mjd') and band magnitude ('mag')
        <u_err, ... >   (dict): Dict with the date ('mjd') and band magnitude error ('mag')
    """

    file_text = f'{kwargs["name"]} {kwargs["redshift"]} {kwargs["ra"]} {kwargs["dec"]}\n'

    for band in 'ugriz':
        filter_name = f'sdss_{band}'
        if band not in kwargs:
            continue

        file_text += f'filter {filter_name}\n'

        iter_data = zip(kwargs[band]['mjd'], kwargs[band]['mag'], kwargs[band + '_err']['mag'])
        for mjd, mag, error in iter_data:
            file_text += f'{mjd} {mag} {error}\n'

    with open(out_path, 'w') as ofile:
        ofile.write(file_text)
