# -*- coding: utf-8 -*-
# filename: text.py
# Copyright 2018 Stefano Costa <steko@iosa.it>
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
import iosacal
from iosacal import ugettext as _, change_lang

from textwrap import indent

def single_text(calibrated_age, BP='bp', lang="en"):
    '''Output calibrated age as simple Markdown text to the terminal.'''
    _ = change_lang(lang)
    formatted_intervals = dict()
    for a, i in calibrated_age.intervals.items():
        formatted_intervals[a] = indent('{:{fmt}}'.format(i, fmt=BP), '* ')
    calibrated_age_max = max(list(calibrated_age.intervals[95]),
            key=lambda i: i.conf_perc)
    calibrated_age_average = abs(calibrated_age_max.from_year + calibrated_age_max.to_year) / 2
    radiocarbon_sample = calibrated_age.radiocarbon_sample.id
    if 'Combined' in radiocarbon_sample:
        radiocarbon_sample = ""
    output = '''
## {radiocarbon_sample_id} **{radiocarbon_sample_date} ± {sigma}**

{calibration_text}: {0.radiocarbon_sample.id}: {average_age} {years_ago_text}

### {calibrated_age_text}:

#### 68.2% {probability_text}

{1[68]}

#### 95.4% {probability_text}

{1[95]}
'''.format(calibrated_age,
        formatted_intervals,
        iosacal.__VERSION__,
        radiocarbon_sample_date=int(round(calibrated_age.radiocarbon_sample.date)),
        sigma=round(calibrated_age[0].radiocarbon_sample.sigma),
        calibration_text=_("Calibration"),
        calibrated_age_text=_("Calibrated calendar age"),
        probability_text=_("Probability"),
        years_ago_text=_("year ago"),
        cal_age_max=calibrated_age_max,
        average_age=int(round(calibrated_age_average)),
        radiocarbon_sample_id=radiocarbon_sample)
    return output
