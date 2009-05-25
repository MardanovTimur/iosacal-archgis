# -*- coding: utf-8 -*-
# filename: util.py
# Copyright 2009 Stefano Costa <steko@iosa.it>
#
# This file is part of GNUCal.

# GNUCal is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# GNUCal is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNUCal.  If not, see <http://www.gnu.org/licenses/>.

from gnucal import hpd


def ad_bc_prefix(year, BP):
    '''Return a string with BC/AD prefix and the given year.'''
    if BP is False:
        if year > 0:
            return "AD %d" % year
        else:
            return "BC %d" % abs(year)
    else:
        return "BP %d" % year

def interval_to_string(interval, calibrated_curve, BP):
    '''Return a string describing the interval with probability percent.'''

    # TODO this should be rather a method of the core.ConfidenceInterval
    # object (to be written yet)

    i = [ad_bc_prefix(year, BP) for year in interval]
    percent = hpd.confidence_percent(interval, calibrated_curve) * 100
    return u' %s ‒ %s (%2.1f %%)\n' % (i[0], i[1], percent)
