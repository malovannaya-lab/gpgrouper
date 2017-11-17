
from setuptools import setup, find_packages
def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

package_version = calculate_version('./pygrouper/_version.py')

setup(
    name='PyGrouper',
    version=package_version,
    py_modules=['pygrouper'],
    packages=find_packages(),
    # packages=['pygrouper',
    #           'pygrouper.pygrouper',
    #           'pygrouper.cli',
    #           'pygrouper.auto_grouper',
    #           'pygrouper.tests',
    # ],
    install_requires=[
        'Click',
        'RefProtDB>=0.1.1'
    ],
    entry_points="""
    [console_scripts]
    pygrouper=pygrouper.cli:cli
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
    url='https://github.com/asalt/pygrouper',
    author='Alexander Saltzman',
    author_email='saltzman@bcm.edu,'


)
