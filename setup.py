from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='emmetrop',
      version=version,
      description="A biological analysis of a schematic eye and natural images",
      long_description=""" """,
      install_requires=["Sphinx",
                        "numpy",
                        "scipy",
                        "matplotlib",
                        "tables"
                        ],
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='Brian Schmidt',
      author_email='bps10@uw.edu',
      url='',
      dependency_links = [
                          "python.org"
                          ],
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
#      package_data={'':'LICENSE'},
      include_package_data=True,
      zip_safe=True,
      entry_points="""
      # -*- Entry points: -*-
      """,
      use_2to3 = True,
      )
