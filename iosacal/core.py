#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: core.py
# Copyright 2016 Stefano Costa <steko@iosa.it>
# Copyright 2017 Mario Gutiérrez-Roig <mariogutierrezroig@gmail.com>
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

import pkg_resources
import sys
from csv import reader
from math import exp, pow, sqrt

import numpy as np

from iosacal.hpd import hpd_interval


def calibrate(f_m, sigma_m, f_t, sigma_t):
    r'''Calibration formula as defined by Bronk Ramsey 2008.

    .. math::

       P(t) \propto \frac{\exp \left[-\frac{(f_m - f(t))^2}{2 (\sigma^2_{fm} + \sigma^2_{f}(t))}\right]}{\sqrt{\sigma^2_{fm} + \sigma^2_{f}(t)}}

See doi: 10.1111/j.1475-4754.2008.00394.x for a detailed account.'''

    sigma_sum = pow(sigma_m, 2) + pow(sigma_t, 2)
    P_t = ( exp( - pow(f_m - f_t, 2 ) /
                   ( 2 * ( sigma_sum ) ) ) / sqrt(sigma_sum) )
    return P_t


def confint(boot_set):
    ''' Calculates the Confidence Intervals of a Bootstrap set.

    This function returns the 68.27% C.I. (1-sigma) and 95.44% C.I. (2-sigma)
    of a set of simulated curves.

    Arguments:
        boot_set -- List with simulated curves as (x,y).

    '''

    mini = np.int(np.min([np.min(x.T[0]) for x in boot_set]))
    maxi = np.int(np.max([np.max(x.T[0]) for x in boot_set]))
    Nboot = len(boot_set)

    yrange = []
    for spd in boot_set:
        y = np.array([[val[0], val[1]] for val in spd if (val[0] >= mini) and (val[0] <= maxi)])
        yy = np.lib.pad(y.T[1], (np.int(y.T[0][0] - mini), np.int(maxi - y.T[0][-1])), 'constant', constant_values=0)
        yrange.append(yy)

    yrange = np.array(yrange).T

    CIxrange = np.arange(mini, maxi + 1)
    CI95sup = []
    CI68sup = []
    CImed = []
    CI68inf = []
    CI95inf = []
    for val in yrange:
        CI95sup.append(np.sort(val)[np.int(0.9772 * Nboot)])
        CI68sup.append(np.sort(val)[np.int(0.8413 * Nboot)])
        CImed.append(np.sort(val)[np.int(0.5 * Nboot)])
        CI68inf.append(np.sort(val)[np.int(0.1587 * Nboot)])
        CI95inf.append(np.sort(val)[np.int(0.0228 * Nboot)])

    return np.array([CIxrange, CI95inf, CI68inf, CImed, CI68sup, CI95sup])


def spdsum(spdlist, norm=True):
    ''' Sums several SPDs stored in a list.

    Arguments:
        norm -- If 'True' final SPD will be normalized.

    '''

    spd = spdlist[0].copy()

    maxi = np.max([np.max(np.array(x).T[0]) for x in spdlist])
    mini = np.min([np.min(np.array(x).T[0]) for x in spdlist])

    # Creates a void vector where perform the sum
    xrange = np.arange(mini, maxi + 1, 1)
    yrange = np.zeros((int)(maxi) + 1 - (int)(mini))

    for d in spdlist:
        # Reshapes every array by adding zeros at the head and tail
        y = np.lib.pad(d.T[1], (np.int(d.T[0][0] - mini), np.int(maxi - d.T[0][-1])), 'constant', constant_values=0)

        # Summing over all histogram
        yrange += y

    # Normalizating the SPD in calendar scale
    if norm == True:
        yrange = yrange / np.sum(yrange)

    spd.resize(len(xrange), 2, refcheck=False)

    spd.T[0] = xrange
    spd.T[1] = yrange

    spd.ndates = np.sum([x.ndates for x in spdlist])

    return spd


class CalibrationCurve(np.ndarray):
    '''A radiocarbon calibration curve.

    Calibration data is loaded at runtime from source data files, and
    exposed a ``numpy.ndarray`` object.

    Implementation from
    http://docs.scipy.org/doc/numpy/user/basics.subclassing.html

    '''

    def __new__(cls, curve_filename):
        title = open(curve_filename, encoding='latin-1').readline().strip('#\n')
        _genfrom = np.genfromtxt(curve_filename, delimiter=',')
        # linear interpolation
        ud_curve = np.flipud(_genfrom)  # the sequence must be *increasing*
        curve_arange = np.arange(ud_curve[0,0],ud_curve[-1,0],1)
        values_interp = np.interp(curve_arange, ud_curve[:,0], ud_curve[:,1])
        stderr_interp = np.interp(curve_arange, ud_curve[:,0], ud_curve[:,2])
        ud_curve_interp = np.array([curve_arange, values_interp, stderr_interp]).transpose()
        _darray = np.flipud(ud_curve_interp)  # back to *decreasing* sequence
        # We cast _darray to be our class type
        obj = np.asarray(_darray).view(cls)
        # add the new attribute to the created instance
        obj.title = title
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.title = getattr(obj, 'title', None)

    def mixing(self, curve2_name, P, D=0, deltaR=0, err_deltaR=0):
        '''Transforms the calibration curve by mixing with another curve.

        The mathematical equation for the mixing is the same than in
        OxCal, as indicated in:

        http://c14.arch.ox.ac.uk/oxcal3/math_ca.htm#mix_curves

        Curve 1 (self): R_1 +/- E_1
        Curve 2 : R_2 +/- E_2
        Mixing : R_m +/- E_m

        where,

        R_m = (1 - P) * R_1 + P * R_2
        E_m = sqrt(((1 - P) * E_1)^2 + (P * E_2)^2 + (D * (R_1 - R_2))^2)

        The local corrections fo the reservoir effect are applicated over the
        calibration curve as:

        R_2_corrected(t) = R_2(t) + deltaR
        E_2_corrected(t) = sqrt((E_2(t))^2 + err_deltaR^2)

        Arguments:
            curve2_name -- The name of the second curve to mix
            P -- Proportion of the second curve
            D -- Error of the proportion
            deltaR -- Reservoir Effect
            err_deltaR -- Error in deltaR

        '''
        curve1 = self

        curve2_path = pkg_resources.resource_filename("iosacal", "data/%s.14c" % curve2_name)
        curve2 = CalibrationCurve(curve2_path)

        # Reservoir Effects
        if deltaR > 0:
            curve2.T[1] += deltaR
            curve2.T[2] = np.sqrt(np.power(curve2.T[2], 2) + np.power(err_deltaR, 2))

        self.T[0] = curve1.T[0]
        self.T[1] = (1. - P) * curve1.T[1] + P * curve2.T[1]
        self.T[2] = np.sqrt(np.power((1. - P) * curve1.T[2], 2) + np.power(P * curve2.T[2], 2) + np.power(
            D * (curve1.T[1] - curve2.T[1]), 2))

        # add the new attribute to the created instance
        self.title = "Mixed_curve"

    def __str__(self):
        return "CalibrationCurve( %s )" % self.title


class RadiocarbonDetermination(object):
    '''A radiocarbon determination as reported by the lab.'''

    def __init__(self, date, sigma, id):
        self.date  = date
        self.sigma = sigma
        self.id = id

    def calibrate(self, curve, norm=True, cutoff=5):
        '''Perform calibration, given a calibration curve.

        Arguments:
            curve -- Calibration curve
            norm -- Normalization in Calendar Scale
            cutoff -- How much of radiocarbon gaussian range are we considering (in sigmas)

        '''

        if not isinstance(curve, CalibrationCurve):
            curve_filename = pkg_resources.resource_filename("iosacal", "data/%s.14c" % curve)
            curve = CalibrationCurve(curve_filename)

        _calibrated_list = []

        # We apply the cutoff for the gaussian range in 14C scale
        idx = np.where((curve.T[1] > self.date - cutoff * self.sigma) & (curve.T[1] < self.date + cutoff * self.sigma))
        idxmin = np.min(idx)
        idxmax = np.max(idx)

        for i in curve[idxmin:idxmax]:
            f_t, sigma_t = i[1:3]
            ca = calibrate(self.date, self.sigma, f_t, sigma_t)
            _calibrated_list.append((i[0],ca))

        if norm == True:
            x = np.array(_calibrated_list).T[0]
            y = np.array(_calibrated_list).T[1] / np.sum(np.array(_calibrated_list).T[1])
            calibrated_curve = np.column_stack((x, y))
        else:
            calibrated_curve = np.array(_calibrated_list)

        cal_age = CalAge(calibrated_curve, self, curve)
        return cal_age

    def __str__(self):
        return "RadiocarbonSample( {id} : {date:.0f} ± {sigma:.0f} )".format(**self.__dict__)


class R(RadiocarbonDetermination):
    '''Shorthand for RadiocarbonDetermination.'''

    pass


class CalAge(np.ndarray):
    '''A calibrated radiocarbon age.

    It is expressed as a probability distribution on the calBP
    calendar scale.

    Implementation from
    http://docs.scipy.org/doc/numpy/user/basics.subclassing.html

    '''

    def __new__(cls, input_array, radiocarbon_sample, calibration_curve):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.radiocarbon_sample = radiocarbon_sample
        obj.calibration_curve = calibration_curve
        obj.intervals = {
            68: hpd_interval(obj,0.318),
            95: hpd_interval(obj,0.046)
        }
        obj.median = np.mean(hpd_interval(obj, 0.5))
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.radiocarbon_sample = getattr(obj, 'radiocarbon_sample', None)
        self.calibration_curve = getattr(obj, 'calibration_curve', None)

    def calendar(self):
        '''Return the calibrated age on the calAD calendar scale.

        This method returns a copy of the calBP array, leaving the
        main object untouched.

        '''

        calendarray = self.copy()
        calendarray[:,0] *= -1
        calendarray[:,0] += 1950
        return calendarray

    def __str__(self):
        return 'CalAge( {radiocarbon_sample.id} based on "{calibration_curve.title}" )'.format(**self.__dict__)


def combine(determinations):
    '''Combine n>1 determinations related to the same event.

    ``determinations`` is an iterable of tuples (mean, error).

    This covers case 1 as described by Ward and Wilson in their
    seminal 1978 paper (DOI: 10.1111/j.1475-4754.1978.tb00208.x)

    '''

    m, s, ids = zip(*[(d.date, d.sigma, d.id) for d in determinations])

    # pooled mean
    pool_m = sum(mi / si**2 for mi, si in zip(m, s)) / \
             sum(1 / si**2 for si in s)

    # standard error on the pooled mean
    pool_s = sqrt(1/sum(1/si**2 for si in s))

    # test statistic
    test = sum((mi - pool_m)**2 / si**2 for mi, si in zip(m, s))

    desc = 'Combined from {} with test statistic {:.3f}'.format(', '.join(ids), test)

    return R(pool_m, pool_s, desc)


class SPD(np.ndarray):
    ''' A Sum of Probability Distributions of calibrated dates '''

    def __new__(cls, caldates, norm=True):
        ''' Initializes the SPD class from calibrated dates

        Arguments:
            caldates -- List with all calibrated dates (after R.calibrate()) in a list
            norm -- If 'True' SPD will be normalized

        '''

        # Finds maximum and minimum value of calibrated dates
        maxi = np.max([np.max(np.array(x).T[0]) for x in caldates])
        mini = np.min([np.min(np.array(x).T[0]) for x in caldates])

        # Creates a void vector where perform the sum
        xrange = np.arange(mini, maxi + 1, 1)
        yrange = np.zeros((int)(maxi) + 1 - (int)(mini))

        for d in caldates:
            # Disentangles and reverses the arrays
            x = np.array(d).T[0][::-1]
            y = np.array(d).T[1][::-1]

            # Reshapes every array by adding zeros at the head and tail
            yy = np.lib.pad(y, (np.int(x[0] - mini), np.int(maxi - x[-1])), 'constant', constant_values=0)

            # Summing over all histogram
            yrange += yy

        # Normalizating the SPD in calendar scale
        if norm == True:
            yrange = yrange / np.sum(yrange)

        obj = np.asarray(np.column_stack((xrange, yrange))).view(cls)
        obj.ndates = len(caldates)

        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.ndates = getattr(obj, 'ndates', None)

    def normalize(self):
        ''' Normalizes the SPD '''

        y = np.array(self.T[1])

        self.T[1] /= np.sum(y)

    def taphcorr(self, function="Surovell", factor=1):
        ''' Performs a taphonomic correction over the SPD

        Taphonomic Correction can be made according two different curves
        indicated by the variable 'function': Surovell's or Williams' correction:

        Surovell's correction:

            n(t) = 5.726442 * 10^6 * (t + 2176.4)^-1.3925309

        Surovell, T. A., Finley, J. B., Smith, G. M., Brantingham, P. J.,
        & Kelly, R. (2009). Correcting temporal frequency distributions for
        taphonomic bias. Journal of Archaeological Science, 36(8), 1715-1724.

        Williams' correction:

            n(t) = 2.107 * 10^7 * (t + 2754)^-1.526

        Williams, A. N. (2012). The use of summed radiocarbon probability
        distributions in archaeology: a review of methods.
        Journal of Archaeological Science, 39(3), 578-589.

        The taphonomic correction is not applied homogeneously, but taking into
        account the percentage of open sites, rock shelters and caves samples
        in each age. This percentage correction is given by the vector
        "factor". If factor=None, then correction is applied homogeneously.

        Arguments:
            factor -- Percentage of correction between 0 and 1.
        '''

        corr = self.copy()

        if function == "Surovell":
            taph = 5726442 * np.power((self.T[0] + 2176.4), -1.3925309)
        elif function == "Williams":
            taph = 21070000 * np.power((self.T[0] + 2754), -1.526)
        else:
            print("ERROR: Invalid argument for 'function' variable")
            sys.exit()

        corr_y = [c / t if t > 0 else c for c, t in zip(np.array(self.T[1]), taph)]
        corr_y = np.array(corr_y) * np.sum(self.T[1]) / np.sum(corr_y)

        corr.T[1] = corr_y * factor + (1. - factor) * self.T[1]

        return corr

    def rollmean(self, m):
        ''' Performs a Rolling Mean of a 2*m+1 window (-m to +m)

        Note:
         - Returns a SPD of smaller size.
         - Variable m should always be greater than 0

        '''

        a = self.copy()
        n = 2 * m + 1
        ret = np.cumsum(a.T[1], dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        a.T[1][m:-m] = ret[n - 1:] / n

        return a[m:-m]

    def simdates(self, N, calcurve_name, c14std, seed=0):
        ''' Simulated dates drawn from SPD curve

        Performs a Bootstrap reseampling on the original SPD and generates a
        new list of simulated dates after backcalibrating. The errors are
        assigned at random from the original distribution.

        Arguments:
            N -- Number of dates to simulate
            seed -- Seed for random generator (if any)
            calcurve_name -- Calibration curve name for backcalibrating
            c14std -- List with all original errors for generating the simulated
                      std.
        '''

        np.random.seed = seed
        SPDcum = np.cumsum(self.T[1])

        calcurve_path = pkg_resources.resource_filename("iosacal", "data/%s.14c" % calcurve_name)
        calcurve = CalibrationCurve(calcurve_path)

        dates = []
        for i in range(N):
            ran = np.random.random()
            idx = (np.abs(SPDcum - ran)).argmin()
            date = self[idx][0]

            idx2 = (np.abs(calcurve.T[0] - date)).argmin()
            c14date = calcurve[idx2][1]

            std = np.random.choice(c14std)
            idname = "sim_" + str(i)
            dates.append([idname, c14date, std])

        return dates

    def simSPD(self, calcurve_name, c14std, seed=0):
        ''' Simulated SPD curve generated from the original SPD

        Performs a Bootstrap reseampling on the original SPD and generates a
        new list of simulated dates that are backcalibrated in the 14C axis,
        then those dates are calibrated again for constructing the SPD. The
        number of simulated dates is the same than the original curve.

        Arguments:
            calcurve_name -- Calibration curve name for backcalibrating
            c14std -- List with all original errors for generating the simulated
                      std.
            seed -- Seed for random generator (if any)
        '''

        simdateslist = self.simdates(self.ndates, calcurve_name, c14std, seed)
        simdates = [RadiocarbonDetermination(x[1], x[2], x[0]) for x in simdateslist]

        # Calibrated date-by-date in order to save memory
        caldates = []
        for x in simdates:
            cal_r = x.calibrate("intcal13", norm=True)
            del cal_r.calibration_curve
            del cal_r.radiocarbon_sample
            del cal_r.intervals
            del cal_r.median
            caldates.append(cal_r)

        simulatedSPD = SPD(caldates)

        return simulatedSPD


class FreqHist(np.ndarray):
    ''' A Frequency Histogram of calibrated dates '''

    def __new__(cls, caldates, bins):
        ''' Initializes the FreqHist class from calibrated dates

        Arguments:
            caldates -- List with all calibrated dates (after R.calibrate()) in a list
            binwidth -- Array with the bins

        '''

        meds = [x.median for x in caldates]

        freqs, bins = np.histogram(meds, bins)

        avbins = 0.5 * (bins[1:] + bins[:-1])

        obj = np.asarray(np.column_stack((avbins, freqs))).view(cls)

        return obj

    def taphcorr(self, function="Surovell", factor=None):
        ''' Performs a taphonomic correction over the Frequency Histrogram

        Taphonomic Correction can be made according two different curves
        indicated by the variable 'function': Surovell's or Williams' correction:

        Surovell's correction:

            n(t) = 5.726442 * 10^6 * (t + 2176.4)^-1.3925309

        Surovell, T. A., Finley, J. B., Smith, G. M., Brantingham, P. J.,
        & Kelly, R. (2009). Correcting temporal frequency distributions for
        taphonomic bias. Journal of Archaeological Science, 36(8), 1715-1724.

        Williams' correction:

            n(t) = 2.107 * 10^7 * (t + 2754)^-1.526

        Williams, A. N. (2012). The use of summed radiocarbon probability
        distributions in archaeology: a review of methods.
        Journal of Archaeological Science, 39(3), 578-589.

        The taphonomic correction is not applied homogeneously, but taking into
        account the percentage of open sites, rock shelters and caves samples
        in each age. This percentage correction is given by the vector
        "factor". If factor=None, then correction is applied homogeneously.

        Arguments:
            factor -- Percentage of correction in each bin.
        '''

        corr = self.copy()

        if function == "Surovell":
            taph = 5726442 * np.power((self.T[0] + 2176.4), -1.3925309)
        elif function == "Williams":
            taph = 21070000 * np.power((self.T[0] + 2754), -1.526)
        else:
            print("ERROR: Invalid argument for 'function' variable")
            sys.exit()

        if factor is None:
            taph = np.array(taph)
        else:
            taph = np.array(taph) * factor

        corr_y = [c / t if t > 0 else c for c, t in zip(np.array(self.T[1]), taph)]
        corr_y = np.array(corr_y) * np.sum(self.T[1]) / np.sum(corr_y)

        corr.T[1] = corr_y

        return corr