import sys
import platform
import os
import shutil
import argparse
sys.path.append(os.getcwd())
from utility import DownloadUrl, UnTar, RunShellCmd, CloneGitRepo


parser = argparse.ArgumentParser(description="""
Install the HDF5 library manually for the OPS libary!\n
Please set the compiler environment before calling this tool, e.g.,
source /opt/intel/oneapi/setvars.sh for using the Intel oneAPI compiler.\n
Sudo might be needed when the INSTALLDIR is pointed to a non-user folder\n
""",formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-v", "--verbose", help="Display detailed output during compilation.", default=False,
                    action="store_true")
parser.add_argument("-p", "--parallel", help="Enable the parallel mode.", default=False,
                    action="store_true")
parser.add_argument("-f", "--fortran", help="Enable the Fortran library.", default=False,
                    action="store_true")
parser.add_argument("-I", "--InstallDir", type=str, default="/usr/local/",
                    help="Set the directory for installation.")
parser.add_argument("-B", "--BuildDir", type=str, default="BuildHDF5",
                    help="Set the directory for building.")
parser.add_argument("-T", "--test", default=False, action="store_true", help="Run testings.")
parser.add_argument("-D", "--download", default=True, action="store_false", help="Download the library.")
args = parser.parse_args()

system = platform.system()
parentPath = os.getcwd()
cpuNum = os.cpu_count()
buildPath = args.BuildDir
installDir = args.InstallDir
fortran = args.fortran
verbose = args.verbose
parallel = args.parallel
test = args.test
download = args.download
if ("Windows" in system):
    print("Windows is not supported at this moment")
try:
    if (not os.path.exists(buildPath)):
        os.mkdir(buildPath)
        print ("Successfully created the directory %s " % buildPath)
except OSError:
    print ("Creation of the directory %s failed" % buildPath)

os.chdir(buildPath)
hdf5Dir = "CMake-hdf5-1.10.7"
url="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/CMake-hdf5-1.10.7.tar.gz"
if download:
    UnTar(DownloadUrl(url))
os.chdir(hdf5Dir)
if not os.path.exists("build"):
    os.mkdir("build")
os.chdir("build")

# here the Z_LIB and SZIP support are disabled since CMake cannot find them.
baseCmd = "-DCMAKE_BUILD_TYPE:STRING=Release  -DHDF5_ENABLE_SZIP_SUPPORT:BOOL=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=OFF -DCMAKE_INSTALL_PREFIX='" + installDir
cmakeCmd = 'cmake ../hdf5-1.10.7'
cmakeCmd = cmakeCmd.split()
cmakeCmd.append('-G')
cmakeCmd.append("Unix Makefiles")
for para in baseCmd.split():
    cmakeCmd.append(para)
if parallel:
    cmakeCmd.append("-DHDF5_ENABLE_PARALLEL=ON")
    cmakeCmd.append("-DHDF5_BUILD_CPP_LIB=OFF")
if verbose:
    cmakeCmd.append("-DCMAKE_VERBOSE_MAKEFILE=ON")
if fortran:
    cmakeCmd.append("-DHDF5_BUILD_FORTRAN=ON")
print(cmakeCmd)
import subprocess
RunShellCmd(cmakeCmd)
RunShellCmd("cmake --build . --config Release -j "+ str(cpuNum))
if test:
    RunShellCmd("ctest . -C Release -j "+ str(cpuNum))
RunShellCmd("cmake --install . --config Release")
os.chdir(parentPath)
try:
    shutil.rmtree(buildPath)
except:
    raise
