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

from iosacal import util


def single_text(calibrated_age, BP):
    '''Output calibrated age as text to the terminal.'''

    output = '''
# IOSACal v0.2

Calibration of {0.radiocarbon_sample.id}: {0.radiocarbon_sample.date} ± {0.radiocarbon_sample.sigma} BP

## Calibrated age

{0.calibration_curve.title}

68.2% probability
{0.intervals[68]}
95.4% probability
{0.intervals[95]}
'''.format(calibrated_age)

    return output
