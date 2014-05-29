# -*- coding: utf-8 -*-
# filename: text.py
# Copyright 2009-2010, 2013-2014 Stefano Costa <steko@iosa.it>
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

from string import Template
from iosacal import util


def text_dict(calibrated_age, BP=True):
    '''Return a dictionary with the meaningful pieces of information.

    This is useful for the web interface and for other custom output.'''

    calibrated_curve = calibrated_age
    f_m = calibrated_age.radiocarbon_sample.date
    sigma_m = calibrated_age.radiocarbon_sample.sigma
    rs_id = calibrated_age.radiocarbon_sample.id
    calibration_curve = calibrated_age.calibration_curve
    calibration_curve_title = calibrated_age.calibration_curve.title
    intervals68 = calibrated_age.intervals68
    intervals95 = calibrated_age.intervals95

    string68 = "".join(
        util.interval_to_string(
            itv, calibrated_curve, BP
            ) for itv in intervals68
        )
    string95 = "".join(
        util.interval_to_string(
            itv, calibrated_curve, BP
            ) for itv in intervals95
        )

    calibrated_data = {
        'calibrated_curve': calibrated_curve,
        'f_m': f_m,
        'sigma_m': sigma_m,
        'rs_id': rs_id,
        'calibration_curve': calibration_curve,
        'calibration_curve_title': calibration_curve_title,
        'intervals68': string68,
        'intervals95': string95,
        'BP': BP,
        }

    return calibrated_data


def single_text(calibrated_age, BP):
    '''Output calibrated age as text to the terminal.'''

    output = Template('''
============
IOSACal v0.1
============
    d = text_dict(calibrated_age, BP)

Radiocarbon sample
------------------

$rs_id: $f_m ± $sigma_m BP

Calibrated age
--------------

$calibration_curve_title

68.2% probability
$intervals68
95.4% probability
$intervals95
''')

    return output.substitute(d)
