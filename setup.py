#!/usr/bin/env python
#https://github.com/pypa/sampleproject/blob/main/setup.py
# See:
# https://packaging.python.org/guides/distributing-packages-using-setuptools/
# https://github.com/pypa/sampleproject

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
# pyCSAMT is a far field basic  open source software of controlled 
# source audio-frequency magnetotellurics 
# for standard data processing , modeling and geophysical interpretation 
# enhancement.  
# """
with open(os.path.join(os.path.abspath('.'), 'README.md'), 'r') as fm:
    # LONG_DESCRIPTION ="""{}""".format(
    #     ' '.join([descp for descp in fm.readlines() ]))
    LONG_DESCRIPTION =fm.read()
try: 
    import pycsamt  # noqa
    VERSION = pycsamt.__version__
except: VERSION ='1.2.1'
# The advantage of setuptools is that EXE wrappers are created on Windows,
# which allows Tab-completion for the script names from the system Scripts
# folder.
# Add names of scripts here. You can also specify the function to call
# by adding :func_name after the module name, and the name of the script
# can be customized before the equals sign. Actually some scripts may not work 
# because scripts is under designing with Tkinter  graphical interface 
# scripts were pseudocreated for waiting their implementations 

setup_kwargs = {}
setup_kwargs['entry_points'] = {
                    'console_scripts':[
    'occam2d_build_in = pycsamt.gui.oc2d_bdin:main',
    'write_avg2edi= pycsamt.gui.wa2edi:main',
    'write_avg2j= pycsamt.gui.wa2j:main',
    'corrected_edi = pycsamt.gui.corrected_edi:main',
    'plot_model_oc2d = pycsamt.gui.p_moc2d:main',
    'plot_pseudolog = pycsamt.gui.p2log:main',
    'write_occam2oasis= pycsamt.gui.oas2f:main',
    'write_occam2golden= pycsamt.gui.gs2f:main',
    'write_iter2dat = pycsamt.gui.wi2d:main',
    'write_drillhole= pycsamt.gui.cmake_dh:main', 
    'correctedi = pycsamt.cli.correctedi:main', 
    'avg2edi= pycsamt.cli.avg2edi:main', 
    'j2edi=pycsamt.cli.j2edi:main', 
    'misfit2d=pycsamt.cli.misfit2d:main', 
    'nm=pycsamt.cli.nm:main', 
    'penetration1d=pycsamt.cli.penetration1d:main', 
    'penetration2d=pycsamt.cli.penetration2d:main', 
    'pseudostratigraphic=pycsamt.cli.pseudostratigraphic:main', 
    'rewriteedi=pycsamt.cli.rewriteedi:main', 
    'rms=pycsamt.cli.rms:main',
    'myconfigfile=pycsamt.cli.myconfigfile:main', 
    'pseudocrossresistivityandphase=pycsamt.cli.pseudocrossresistivityandphase:main', 
    'staticshift=pycsamt.cli.staticshift:main', 
    'fitforward=pycsamt.cli.fitforward:main', 
    'occambuildinputs=pycsamt.cli.occambuildinputs:main'
    ]
}                
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
    'pycsamt.core',
    'pycsamt.gui',
    'pycsamt.cli',
    'pycsamt.processing',
    'pycsamt.geodrill',
    'pycsamt.view',
    'pycsamt.modeling',
    'pycsamt.utils',
]
# force install pycsamt. Once pycsamt is installed , pyyaml and pyproj 
# should already installed too. 
     
setup_kwargs['install_requires'] = [
    'numpy>=1.8.1',
    'scipy>=0.14.0',
    'matplotlib',
    'mtpy >=1.1.0',
    'pyyaml',
    'pyproj',
    'configparser', 
    'tqdm'
]
                                     
setup_kwargs['python_requires'] ='>=3.7'
authors =["Kouadio K. Laurent, ", 'Rong Liu, ', 
          'Binbin Mi, ','Chum-ning Liu, ', 'Albert O. Malory.']
authors_emails =['etanoyau@gmail.com,', 'liurongkaoyang@126.com,',
                'mibinbin@zju.edu.cn,', 'lifuming001@163.com,','amalory@zju.edu.cn']
setup(
 	name="pycsamt",
 	version=VERSION,
 	author=' '.join([aa for aa in authors]),
    author_email='kkouao@zju.edu.cn',
    maintainer="Kouadio K. Laurent",
    maintainer_email='etanoyau@gmail.com',
 	description="A Python open-source toolkit Audio-frequency Magnetotelluric ",
 	long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/WEgeophysics/pyCSAMT",
    project_urls={
        "API Documentation"  : "https://pycsamt.readthedocs.io/en/master/",
        "Home page" : "https://github.com/WEgeophysics/pyCSAMT/wiki",
        "Bugs tracker": "https://github.com/WEgeophysics/pyCSAMT/issues",
        "Installation guide" : "https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux", 
        "User guide" : "https://github.com/WEgeophysics/pyCSAMT/blob/develop/docs/pyCSAMT%20User%20Guide.pdf",
        },
 	#data_files=[('', ['pycsamt/utils/epsg.npy',]),], #this will install datafiles in wearied palce such as ~/.local/
 	include_package_data=True,
 	license="GNU LESSER GENERAL PUBLIC LICENSE v3",
 	classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        # "Topic :: Software Development :: Build Tools",
        #"License :: OSI Approved :: GNU License",
        "Topic :: Software Development",
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        "Operating System :: OS Independent",
        ],
    keywords="hydrogeophysic, groundwater, exploration, csamt",
    #package_dir={"": "pyCSAMT"},  # Optional
    package_data={'pycsamt': [
                            'utils/p.configlog.yml', 
                            'utils/espg.npy',
                            'geodrill/_geocodes/*.csv', 
                            'geodrill/_geocodes/__memory.pkl', 
                            'geodrill/_geomemo/*.sq3', 
                            'metadata/e.g.data.json', 
                            'metadata/e.g.data.yml', 
                            '_loggerfiles/*.txt',
                            
                            ], 
                    "":[
                        'data/occam2d/*', 
                        'data/drill_example_files/*', 
                        'data/avg/K1.avg', 
                        'data/avg/K1.stn', 
                        'data/avg/K2.avg', 
                        'data/avg/K2.stn',
                        'project_description.md',
                        ]
                  },
    
 	**setup_kwargs
)
























