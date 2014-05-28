#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: plot.py
# Copyright 2009, 2013-2014 Stefano Costa <steko@iosa.it>
#
# This file is part of IOSACal, the IOSA Radiocarbon Calibration Library.

# IOSACal is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# IOSACal is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with IOSACal.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from iosacal import hpd, util

COLORS = {
    'bgcolor': '#e5e4e5',
}

def single_plot(calibrated_age, oxcal=False, output=None, BP=True):

    calibrated_age = calibrated_age
    f_m = calibrated_age.radiocarbon_sample.date
    sigma_m = calibrated_age.radiocarbon_sample.sigma
    radiocarbon_sample_id = calibrated_age.radiocarbon_sample.id
    calibration_curve = calibrated_age.calibration_curve
    intervals68 = calibrated_age.intervals68
    intervals95 = calibrated_age.intervals95
    sample_interval = calibration_curve[:,0].copy() # for determination plot

    # adjust plot bounds
    min_year, max_year = (50000, -50000)
    min_x = min(calibrated_age[:,0])
    max_x = max(calibrated_age[:,0])
    if min_year > min_x:
        min_year = min_x
    if max_year < max_x:
        max_year = max_x

    # do not plot the part of calibration curve that is not visible
    # greatly reduces execution time \o/
    cutmin = calibration_curve[calibration_curve[:,0]>min_x]
    cutmax = cutmin[cutmin[:,0]<max_x]
    calibration_curve = cutmax

    if BP is False:
        if min_year < 0 and max_year > 0:
            ad_bp_label = "BC/AD"
        elif min_year < 0 and max_year < 0:
            ad_bp_label = "BC"
        elif min_year > 0 and max_year > 0:
            ad_bp_label = "AD"
    else:
        ad_bp_label = "BP"

    string68 = "".join(
        util.interval_to_string(
            itv, calibrated_age, BP
            ) for itv in intervals68
        )
    string95 = "".join(
        util.interval_to_string(
            itv, calibrated_age, BP
            ) for itv in intervals95
        )

    fig = plt.figure(figsize=(12,8))
    fig.clear()
    ax1 = plt.subplot(111)
    ax1.set_axis_bgcolor(COLORS['bgcolor'])
    plt.xlabel("Calibrated age (%s)" % ad_bp_label)
    plt.ylabel("Radiocarbon determination (BP)")
    plt.text(0.5, 0.95,r'%s: $%d \pm %d BP$' % (radiocarbon_sample_id, f_m, sigma_m),
         horizontalalignment='center',
         verticalalignment='center',
         transform = ax1.transAxes,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))
    plt.text(0.75, 0.80,'68.2%% probability\n%s\n95.4%% probability\n%s' \
                 % (string68, string95),
         horizontalalignment='left',
         verticalalignment='center',
         transform = ax1.transAxes,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))
    plt.text(0.0, 1.0,'IOSACal v0.2; %s' % calibration_curve.title,
         horizontalalignment='left',
         verticalalignment='bottom',
         transform = ax1.transAxes,
         size=10,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))

    # Calendar Age

    ax2 = plt.twinx()

    if oxcal is True:
        # imitate OxCal
        ax1.set_axis_bgcolor('white')
        ax2.fill(
            calibrated_age[:,0],
            calibrated_age[:,1] + max(calibrated_age[:,1])*0.3,
            'k',
            alpha=0.3,
            label='Calendar Age'
            )
        ax2.plot(
            calibrated_age[:,0],
            calibrated_age[:,1],
            'k',
            alpha=0
            )
    else:
        ax2.fill(
            calibrated_age[:,0],
            calibrated_age[:,1],
            'k',
            alpha=0.3,
            label='Calendar Age'
            )
        ax2.plot(
            calibrated_age[:,0],
            calibrated_age[:,1],
            'k',
            alpha=0
            )

    ax2.set_ybound(min(calibrated_age[:,1]),max(calibrated_age[:,1])*3)
    ax2.set_xbound(min(calibrated_age[:,0]),max(calibrated_age[:,0]))
    ax2.set_axis_off()

    # Radiocarbon Age
    sample_curve = mlab.normpdf(sample_interval, f_m, sigma_m)

    ax3 = plt.twiny(ax1)
    ax3.fill(
        sample_curve,
        sample_interval,
        '1.0',
        alpha=0.8
        )
    ax3.set_xbound(0,max(sample_curve)*4)
    ax3.set_axis_off()

    # Calibration Curve

    mlab_low = calibration_curve[:,1] - calibration_curve[:,2]
    mlab_high = calibration_curve[:,1] + calibration_curve[:,2]

    xs, ys = mlab.poly_between(calibration_curve[:,0],
                               mlab_low,
                               mlab_high)
    ax1.fill(xs, ys, fc='#000000', ec='none', alpha=0.15)
    ax1.plot(calibration_curve[:,0], calibration_curve[:,1], '#000000', alpha=0.5)

    # Confidence intervals

    if oxcal is True:
        for i in intervals68:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.05,
                ymax=0.07,
                facecolor='none',
                alpha=0.8)
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.068,
                ymax=0.072,
                facecolor='w',
                edgecolor='w',
                lw=2)
        for i in intervals95:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.025,
                ymax=0.045,
                facecolor='none',
                alpha=0.8)
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.043,
                ymax=0.047,
                facecolor='w',
                edgecolor='w',
                lw=2)
    else:
        for i in intervals68:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0,
                ymax=0.02,
                facecolor='k',
                alpha=0.5)
        for i in intervals95:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0,
                ymax=0.02,
                facecolor='k',
                alpha=0.5)

    # FIXME the following values 10 and 5 are arbitrary and could be probably
    # drawn from the f_m value itself, while preserving their ratio
    ax1.set_ybound(f_m - sigma_m * 15, f_m + sigma_m * 5)
    ax1.set_xbound(min(calibrated_age[:,0]),max(calibrated_age[:,0]))
    ax1.invert_xaxis()          # if BP == True

    if output:
        plt.savefig(output)


def multi_plot(calibrated_ages,name,oxcal=False):

    # Define the legend and descriptive text

    min_year, max_year = (50000, -50000)

    for calibrated_curve in calibrated_ages:
        if min_year < min(calibrated_curve.array[:,0]):
            pass
        else:
            min_year = min(calibrated_curve.array[:,0])
        if max_year > max(calibrated_curve.array[:,0]):
            pass
        else:
            max_year = max(calibrated_curve.array[:,0])

    if calibrated_ages[0].BP is False:
        if min_year < 0 and max_year > 0:
            ad_bp_label = "BC/AD"
        elif min_year < 0 and max_year < 0:
            ad_bp_label = "BC"
        elif min_year > 0 and max_year > 0:
            ad_bp_label = "AD"
    else:
        ad_bp_label = "BP"

    fig = plt.figure(1)
    plt.suptitle("%s" % name )
    plt.suptitle("Calibrated date (%s)" % ad_bp_label, y = 0.05)

    for n, calibrated_curve in enumerate(calibrated_ages):
        fignum = 1 + n
        numrows = len(calibrated_ages)
        ax1 = fig.add_subplot(numrows,1,fignum)

        # Calendar Age

        ax1.fill(
            calibrated_curve.array[:,0],
            calibrated_curve.array[:,1],
            'k',
            alpha=0.3,
            label='Calendar Age'
            )
        ax1.plot(
            calibrated_curve.array[:,0],
            calibrated_curve.array[:,1],
            'k',
            alpha=0
            )
        ax1.set_ybound(
            min(calibrated_curve.array[:,1]),
            max(calibrated_curve.array[:,1])*2
            )
        ax1.set_xbound(min_year, max_year)
        #ax1.set_axis_off()

        # Confidence intervals

        for i in calibrated_curve.intervals95:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.6,
                ymax=0.7,
                facecolor='k',
                alpha=0.5)
        for i in calibrated_curve.intervals68:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.6,
                ymax=0.7,
                facecolor='k',
                alpha=0.8)

    plt.savefig('image_%s.png' % name )
    fig = plt.gcf()
    fig.clear()
