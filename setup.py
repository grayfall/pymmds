from setuptools import setup
import sys


if sys.version_info < (3, 5):
    print('This package requires Python >= 3.5')
    sys.exit(1)


setup(
    name='pymmds',
    version='1.0b',
    packages=['mmds'],
    url='',
    license='MIT',
    author='Ilia Korvigo',
    author_email='ilia.korvigo@gmail.com',
    description='A Python package for active metric MDS.',
    install_requires=['numpy>=1.14.0', 'pandas>=0.22.0']
)