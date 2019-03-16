from iosacal import R, combine, iplot
from iosacal.plot import single_plot, stacked_plot
from iosacal.text import single_text
from matplotlib import pyplot as plt


if __name__ == '__main__':
    r = R(2012, 93.2, 'P-769')
    r2 = R(2013, 92, 'P-2123')
    r3 = R(2013, 92, 'P-2123')
    r4 = R(2013, 92, 'P-2123')
    combined = combine([r, r2, r3, r4] * 15)
    calibrated = combined.calibrate('intcal13')
    iplot(calibrated)
    plt.show()
    #  r = R(2012, 93.2, 'P-769')
    #  cal_r = r.calibrate('intcal13')
    #  plt.show()
