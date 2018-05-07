from setuptools import setup, find_packages

version = '0.4.0'

setup(name='iosacal',
      version=version,
      description="IOSACal is a radiocarbon (14C) calibration program",
      long_description="""\
""",
      classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering',
        ],
      keywords='radiocarbon calibration',
      author='Stefano Costa',
      author_email='steko@iosa.it',
      url='http://c14.iosa.it/',
      license='GNU GPLv3',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      package_data={'iosacal': ['data/*.14c']},
      include_package_data=True,
      zip_safe=False,
      install_requires=[
         # -*- Extra requirements: -*-
        'numpy >= 1.14.0',
        'scipy >= 1.1.0',
        'matplotlib >= 2.2.0'
      ],
      entry_points= {
        'console_scripts': [
            'iosacal = iosacal.cli:main',
            ]
        },
      )
