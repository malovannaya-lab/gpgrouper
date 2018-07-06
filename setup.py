import os
from setuptools import setup, find_packages
def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

_version_f = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                          'gpgrouper', '_version.py'
)
package_version = calculate_version(_version_f)

setup(
    name='gpGrouper',
    version=package_version,
    py_modules=['gpgrouper'],
    packages=find_packages(),
    # packages=['pygrouper',
    #           'pygrouper.pygrouper',
    #           'pygrouper.cli',
    #           'pygrouper.auto_grouper',
    #           'pygrouper.tests',
    # ],
    install_requires=[
        'Click',
        'RefProtDB>=0.1.1',
        'numpy',
        'pandas',
        'mock',

    ],
    entry_points="""
    [console_scripts]
    gpgrouper=gpgrouper.cli:cli
    """,
    # package_data={
    #     'pygrouper': ['base_config.ini']
    # }

    description='Gene level peptide grouping and quantification of bottom up proteomics data',
    license='GPL-3.0',
    classifiers=[
        'Intended Audience :: Proteomics',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    url='https://github.com/asalt/gpgrouper',
    author='Alexander Saltzman',
    author_email='saltzman@bcm.edu,'


)
