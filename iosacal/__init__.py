import gettext
from os.path import join, dirname, abspath
from iosacal.core import R, combine
from iosacal.plot import iplot


i18n_dir = abspath(join(abspath(dirname(__file__)), 'i18n'))
i18n = gettext.translation('iosacal', localedir=i18n_dir,
                           languages=['en'])
i18n.install()

ugettext = i18n.gettext


def change_lang(lang):
    i18n = gettext.translation('iosacal', localedir=i18n_dir,
                               languages=[lang])
    i18n.install()
    return i18n.gettext

__VERSION__ = '0.4.1'
