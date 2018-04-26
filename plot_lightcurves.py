#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Plot all supernova classified light curves from SDSS"""

from sdss_sn_analysis import plot_sn_light_curves

if __name__ == '__main__':
    plot_sn_light_curves(out_dir='./light_curves', verbose=True)
