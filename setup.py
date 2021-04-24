#!/usr/bin/env python

import pycsamt
import os 

# Check for setuptools package:

try:
    from setuptools import setup
except ImportError:
    setuptools = False
    from distutils.core import setup
else:
    setuptools = True

# LONG_DESCRIPTION = """
# pyCSAMT is a far field basic  open source software of controlled source audio-frequency magnetotellurics 
# for standard data processing , modeling and geophysical interpretation  enhancement.  
# """
with open(os.path.join(os.path.abspath('.'), 
                       'project_description.md'), 'r') as fm:
    LONG_DESCRIPTION ="""{}""".format(
        ' '.join([descp for descp in fm.readlines() ]))

# The advantage of setuptools is that EXE wrappers are created on Windows,
# which allows Tab-completion for the script names from the system Scripts
# folder.
# Add names of scripts here. You can also specify the function to call
# by adding :func_name after the module name, and the name of the script
# can be customized before the equals sign. Actually some scripts may not work 
# because scripts is under designing with Tkinter  graphical interface 
# scripts were pseudocreated for waiting their implementations 

setup_kwargs = {}
setup_kwargs['entry_points'] = {'console_scripts': 
                    ['occam2d_build_in = pycsamt.gui.oc2d_bdin:main',
                     'write_avg2edi= pycsamt.gui.wa2edi:main',
                     'write_avg2j= pycsamt.gui.wa2j:main',
                     'corrected_edi = pycsamt.gui.corrected_edi:main',
                     'plot_model_oc2d = pycsamt.gui.p_moc2d:main',
                     'plot_pseudolog = pycsamt.gui.p2log:main',
                     'write_occam2oasis= pycsamt.gui.oas2f:main',
                     'write_occam2golden= pycsamt.gui.gs2f:main'
                     'write_iter2dat = pycsamt.gui.wi2d:main',
                     'write_drillhole= pycsamt.gui.cmake_dh:main']}
                     
                        

# But many people will not have setuptools installed, so we need to handle
# the default Python installation, which only has Distutils:

if setuptools is False:
    # Different script specification style for ordinary Distutils:

    setup_kwargs['scripts'] = [
        s.split(' = ')[1].replace('.', '/').split(':')[0] + '.py' for s in 
        setup_kwargs['entry_points']['console_scripts']]
    del setup_kwargs['entry_points']

    # "You must explicitly list all packages in packages: the Distutils will not
    # recursively scan your source tree looking for any directory with an
    # __init__.py file"

setup_kwargs['packages'] = [ 
                            'pycsamt',
                            'pycsamt.etc',
                            'pycsamt.ff',
                            'pycsamt.ff.core',
                            'pycsamt.gui',
                            'pycsamt.ff.processing',
                            'pycsamt.geodrill',
                            'pycsamt.geodrill.geoCore',
                            'pycsamt.geodrill.geoDB',
                            'pycsamt.geodrill.geoDB.sql_utils',
                            'pycsamt.geodrill.geoDB.sql_utils.sql_DB',
                            'pycsamt.viewer',
                            'pycsamt.modeling',
                            'pycsamt.utils',
                            ]
# force install mtpy . Once mtpy is installed , pyyaml and pyproj 
# should already installed too. 
     
setup_kwargs['install_requires'] = ['numpy>=1.8.1',
                                     'scipy>=0.14.0',
                                     'matplotlib',
                                     'mtpy >=1.1.0',
                                     'pyyaml',
                                     'pyproj',
                                     'configparser']
                                     
setup_kwargs['python_requires'] ='>=3.6'

authors =["Kouadio K. Laurent", 'Rong Liu', 
          'Binbin Mi','Chum-ning Liu', 'Albert O. Malory']
authors_emails =['etanoyau@gmail.com', 'liurongkaoyang@126.com',
                'mibinbin@zju.edu.cn', 'lifuming001@163.com','amalory@zju.edu.cn']
setup(
	name="pycsamt",
	version=pycsamt.__version__,
	author=' '.join([aa for aa in authors]),
    author_email='kkouao@zju.edu.cn',
    maintainer="Kouadio K. Laurent",
    maintainer_email='etanoyau@gmail.com',
	description="A Python open-source toolkit for standard Controlled Source Audio-frequency MagnetoTellurics (CSAMT)",
	long_description=LONG_DESCRIPTION,
    url="https://github.com/WEgeophysics/pycsamt",
	#data_files=[('', ['pycsamt/utils/epsg.npy',]),], #this will install datafiles in wearied palce such as ~/.local/
	include_package_data=True,
	license="GNU GENERAL PUBLIC LICENSE v3",
	classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        ],
	**setup_kwargs)
