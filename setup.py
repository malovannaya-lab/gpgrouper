
from setuptools import setup, find_packages
def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

package_version = calculate_version('./pygrouper/cli.py')

setup(
    name='PyGrouper',
    version=package_version,
    #py_modules=['pygrouper'],
    packages=find_packages(),
    install_requires=[
        'Click', 'pandas',
    ],
    entry_points="""
    [console_scripts]
    pygrouper=pygrouper.cli:cli
    """,
    package_data={
        'pygrouper': ['base_config.ini']
    }
)
