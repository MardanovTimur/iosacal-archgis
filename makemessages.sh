#!/bin/bash
find . -iname "*.py" | xargs xgettext -j --from-code=UTF-8 -d=iosacal -o i18n/en/LC_MESSAGES/iosacal.pot
find . -iname "*.py" | xargs xgettext -j --from-code=UTF-8 -d=iosacal -o i18n/ru/LC_MESSAGES/iosacal.pot
find . -iname "*.py" | xargs xgettext -j --from-code=UTF-8 -d=iosacal -o i18n/tat/LC_MESSAGES/iosacal.pot
msgfmt i18n/en/LC_MESSAGES/iosacal.pot -o i18n/en/LC_MESSAGES/iosacal.mo
msgfmt i18n/ru/LC_MESSAGES/iosacal.pot -o i18n/ru/LC_MESSAGES/iosacal.mo
msgfmt i18n/tat/LC_MESSAGES/iosacal.pot -o i18n/tat/LC_MESSAGES/iosacal.mo
