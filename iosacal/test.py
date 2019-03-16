from . import *
from .text import *

if __name__ == '__main__':
    print("test radiocarbon")
    r = R(7505, 93, 'P-769')
    cal_r = r.calibrate('intcal13')
    print(single_text(cal_r))
