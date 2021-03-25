
import os
import shutil

from csamtpy.utils._csamtpylog import csamtpylog

TEST_pyCSAMT_ROOT = os.path.normpath(
    os.path.abspath(
        os.path.dirname(__file__)
#            os.path.dirname(__file__)
        )
    )
#)  # assume tests is on the root level of pyCSAMT

TEST_DIR = os.path.normpath(os.path.abspath(os.path.dirname(__file__)))
TEST_TEMP_DIR = os.path.normpath(os.path.join(TEST_DIR, "temp"))

if not os.path.isdir(TEST_TEMP_DIR):
    os.mkdir(TEST_TEMP_DIR)


def make_temp_dir(dir_name, base_dir=TEST_TEMP_DIR):
    _temp_dir = os.path.normpath(os.path.join(base_dir, dir_name))
    if os.path.isdir(_temp_dir):
        shutil.rmtree(_temp_dir)
    os.mkdir(_temp_dir)
    return _temp_dir

# Main Data directories

EDI_DATA_DIR = os.path.normpath(
    os.path.join(TEST_pyCSAMT_ROOT, 'data/edi'))
EDI_DATACOR_DIR = os.path.normpath(
    os.path.join(TEST_pyCSAMT_ROOT, 'data/correctedEDI'))
AVG_DATA_DIR = os.path.normpath(
    os.path.join(TEST_pyCSAMT_ROOT, 'data/avg'))
DRILL_PARSER_DIR = os.path.normpath(
    os.path.join(TEST_pyCSAMT_ROOT, 'data/drill_example_files'))  
OC2D_DIR = os.path.normpath(
    os.path.join(TEST_pyCSAMT_ROOT, 'data/occam2D'))

# set test logging configure
csamtpylog.load_configure(os.path.join('_logfile', "main_logging_configfile.yml"))
