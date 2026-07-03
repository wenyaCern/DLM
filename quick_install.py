#!/usr/bin/env python3
import sys, os, subprocess, time, multiprocessing

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def yes():
    answer = ''
    while True:
        answer = input()
        if answer == 'abort':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return False
        elif answer.lower() in ("y", "yes"):
            return True
        elif answer.lower() in ("n", "no"):
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return False
        else:
            print(' yes/no (y/n) answers only: ', end='')

def input_number_range(min: int, max: int):
    while True:
        answer = input()
        if answer == 'abort':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return min-1
        elif not answer.isnumeric():
            print(bcolors.WARNING+' The input needs to be numeric, try again: '+bcolors.ENDC, end='')
        else:
            v = int(answer)
            if v < min or v > max:
                print(bcolors.WARNING+' Unlisted version, try again: '+bcolors.ENDC, end='')
            else:
                return v

def is_gsl_include_okay(PATH_GSL_INC):
    required = ['gsl_sf_dawson.h', 'gsl_sf_coulomb.h', 'gsl_sf_bessel.h',
                'gsl_sf_legendre.h', 'gsl_sf_gamma.h']
    candidates = [PATH_GSL_INC, os.path.join(PATH_GSL_INC, 'gsl')]
    for p in candidates:
        try:
            files = os.listdir(p)
        except Exception:
            continue
        if all(f in files for f in required):
            return True
    return False

def is_gsl_lib_okay(PATH_GSL_LIB):
    if PATH_GSL_LIB == '':
        return False
    try:
        files = os.listdir(PATH_GSL_LIB)
    except Exception:
        return False
    for f in files:
        if f in ['libgsl.a', 'libgsl.so', 'libgsl.dylib', 'libgslcblas.a', 'libgslcblas.dylib', 'libgslcblas.so']:
            return True
    return False

def get_gsl_lib(PATH_GSL_LIB):
    libname = ''
    if '/libgsl.a' in PATH_GSL_LIB:
        libname = '/libgsl.a'
    elif '/libgsl.so' in PATH_GSL_LIB:
        libname = '/libgsl.so'
    elif '/libgsl.dylib' in PATH_GSL_LIB:
        libname = '/libgsl.dylib'
    if libname == '':
        return 'libgsl:'
    PATH_GSL_LIB = PATH_GSL_LIB.replace(libname,'')
    PATH_GSL_LIB = PATH_GSL_LIB.replace('\n','')
    return PATH_GSL_LIB

def brew_prefix(pkg: str) -> str:
    try:
        out = subprocess.check_output(["brew", "--prefix", pkg], stderr=subprocess.STDOUT).decode().strip()
        return out
    except Exception:
        return ""

def quick_install(argv):
    print(bcolors.BOLD+'CATS installer (macOS patched)'+bcolors.ENDC)
    print(' To work properly, this file has to be executed from the DLM repository!')
    print(' Type "abort" to terminate the script.')

    SCRIPT_NAME = argv[0]
    SCRIPT_REL_PATH = os.path.dirname(SCRIPT_NAME)
    SCRIPT_PATH = os.path.abspath(SCRIPT_REL_PATH)
    os.chdir(SCRIPT_PATH)

    allowed_types = ["basic","root","dev","compile","clean"]
    if len(argv)!=2 or argv[1] not in allowed_types:
        print('--------------------')
        print('Run this script using '+bcolors.BOLD+bcolors.OKCYAN+'"python3 quick_install_macos.py ARG"'+bcolors.ENDC)
        print('Available ARG: basic, root, dev, compile, clean')
        return

    if argv[1] not in ("compile","clean"):
        print(' Continue? (y/n) ', end='')
        bashrc_lines = []
        if not yes():
            return

    install_mode = argv[1]
    print(' Installation mode: '+install_mode)

    install_lvl = 0
    if install_mode == 'root':
        install_lvl = 1
    elif install_mode == 'dev':
        install_lvl = 2

    # DLM repository root
    PATH_DLM = os.path.abspath(os.getcwd())
    PATH_CATS = PATH_DLM + '/install'
    PATH_CATS_CMAKE = PATH_CATS + '/CMake'

    # compile / clean existing
    if install_mode in ('compile','clean'):
        if not os.path.exists(PATH_CATS_CMAKE):
            print(bcolors.FAIL+f'The installation folder "{PATH_CATS_CMAKE}" does not exist'+bcolors.ENDC)
            return
        os.chdir(PATH_CATS_CMAKE)
        if install_mode == 'compile':
            if os.system('make -j'+str(multiprocessing.cpu_count())):
                print(bcolors.FAIL+'"make" failed'+bcolors.ENDC); return
            if os.system('make install'):
                print(bcolors.FAIL+'"make install" failed'+bcolors.ENDC); return
            print(bcolors.OKGREEN+'CATS was successfully compiled!'+bcolors.ENDC)
        else:
            if os.system('make clean'):
                print(bcolors.FAIL+'"make clean" failed'+bcolors.ENDC); return
            print(bcolors.OKGREEN+'CATS was successfully removed!'+bcolors.ENDC)
        return

    # -------- macOS OpenMP (libomp) ----------
    IS_MAC = (sys.platform == "darwin")
    OMP_PREFIX = ""
    OMP_INC = ""
    OMP_LIB = ""
    if IS_MAC:
        OMP_PREFIX = brew_prefix("libomp")
        if not OMP_PREFIX:
            print(bcolors.FAIL+'Homebrew libomp not found!'+bcolors.ENDC)
            print('Please install: brew install libomp')
            return
        OMP_INC = OMP_PREFIX + "/include"
        OMP_LIB = OMP_PREFIX + "/lib"
        if not os.path.exists(os.path.join(OMP_INC, "omp.h")):
            print(bcolors.FAIL+f'omp.h not found in {OMP_INC}'+bcolors.ENDC)
            return
        print(bcolors.OKGREEN+' OpenMP (libomp) found: '+bcolors.ENDC+OMP_PREFIX)

    # ROOT only if requested
    PATH_ROOT = ""
    root_loaded = False
    if install_lvl >= 1:
        # prefer root-config if available, else ask for thisroot.sh folder
        rc = subprocess.getstatusoutput("command -v root-config")[1]
        if rc and os.path.exists(rc):
            # root prefix = .. from bin/root-config
            PATH_ROOT = os.path.abspath(os.path.join(os.path.dirname(rc), ".."))
            root_loaded = True
        else:
            print(bcolors.FAIL+' ROOT installation not found in PATH!'+bcolors.ENDC)
            print(' Please provide the directory that contains "thisroot.sh" (typically <ROOT>/bin):')
            print(' PATH_ROOT_BIN = ', end='')
            p = input()
            if p == 'abort': print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC); return
            if os.path.exists(os.path.join(p, "thisroot.sh")):
                PATH_ROOT = os.path.abspath(os.path.join(p, ".."))
                root_loaded = True
        if root_loaded:
            print(bcolors.OKGREEN+' ROOT location found: '+bcolors.ENDC+PATH_ROOT)
        else:
            print(bcolors.FAIL+' ROOT requested, but not found.'+bcolors.ENDC)
            return

    # -------- GSL (brew) ----------
    # Prefer brew prefix; fall back to your hardcoded Cellar path
    gsl_prefix = brew_prefix("gsl")
    if gsl_prefix:
        PATH_GSL_INC = gsl_prefix + "/include"
        PATH_GSL_LIB = gsl_prefix + "/lib"
    else:
        PATH_GSL_INC = "/usr/local/Cellar/gsl/2.7.1/include"
        PATH_GSL_LIB = "/usr/local/opt/gsl/lib"

    if not os.path.exists(PATH_GSL_INC):
        print(bcolors.FAIL+'GNU-GSL include folder not found!'+bcolors.ENDC)
        print('Try: brew install gsl')
        return

    GSL_INCLUDE_PATHS = PATH_GSL_INC + " " + os.path.join(PATH_GSL_INC, "gsl")
    print(bcolors.OKGREEN+'GSL include paths set to: '+bcolors.ENDC + GSL_INCLUDE_PATHS)

    if not is_gsl_include_okay(PATH_GSL_INC):
        print(bcolors.FAIL+' GNU-GSL headers not found in '+PATH_GSL_INC+bcolors.ENDC)
        return

    if not is_gsl_lib_okay(PATH_GSL_LIB):
        print(bcolors.FAIL+' GSL lib folder not valid: '+PATH_GSL_LIB+bcolors.ENDC)
        return
    print(bcolors.OKGREEN+' LIBGSL location found: '+bcolors.ENDC+PATH_GSL_LIB)

    # summary
    time.sleep(0.3)
    print('--------------------')
    print(bcolors.OKGREEN+'The automatic setup was successful!'+bcolors.ENDC)
    print(' DLM repository: '+PATH_DLM)
    print(' CATS to be installed in: '+PATH_CATS)
    print(' Using:')
    if install_lvl >= 1:
        print('  ROOT: '+PATH_ROOT)
    print('  GSL (include): '+PATH_GSL_INC)
    print('  GSL (lib): '+PATH_GSL_LIB)
    if IS_MAC:
        print('  OpenMP: '+OMP_PREFIX)
    print(' Are these setting OKAY? (y/n) ', end='')
    if not yes():
        return

    # sources / headers
    OnlyHeader = ['CATSconstants','DLM_Sort','DLM_Histo']

    IncDir = ['CATS', 'CATS_Extentions', 'DLM_CppTools', 'DLM_MathTools']
    if install_lvl >= 1:
        IncDir.append('DLM_RootTools')
    if install_lvl >= 2:
        IncDir.extend(['DLM_FemtoTools', 'CECA'])

    IncFile = {}
    IncFile['CATS'] = ['CATS', 'CATS3D', 'CATStools', 'CATSconstants']

    # IMPORTANT: for macOS basic, exclude ROOT-dependent DLM_WfModel.cpp
    # (your repo version includes TGraph.h)
    ext_files = ['DLM_Ck', 'DLM_ResponseMatrix', 'DLM_CkDecomp', 'DLM_Source', 
                 'DLM_StefanoPotentials', 'DLM_Potentials', 'DLM_CkModels']
    if install_lvl >= 1:
        ext_files.extend(['DLM_CkDecomposition', 'DLM_Fitters'])
        ext_files.append('DLM_WfModel')  # OK in root mode
    else:
        # basic mode: do NOT compile DLM_WfModel.cpp
        pass
    IncFile['CATS_Extentions'] = ext_files

    cpp_files = ['DLM_CppTools','DLM_Sort']
    if install_lvl >= 2:
        cpp_files.append('DLM_OmpTools')
    IncFile['DLM_CppTools'] = cpp_files

    math_files = ['DLM_Histo','DLM_Integration','DLM_Random','DLM_Bessel','DLM_MathFunctions','DLM_Fit','DLM_RootFinder']
    if install_lvl >= 1:
        math_files.append('DLM_Unfold')
    IncFile['DLM_MathTools'] = math_files

    if install_lvl >= 1:
        IncFile['DLM_RootTools'] = ['DLM_CRAB_PM','DLM_SubPads','DLM_DrawingTools','DLM_FSI_Wrapper',
                                    'DLM_HistoAnalysis','DLM_RootWrapper','DLM_MultiFit',
                                    'DLM_DecayMatrix','DLM_RootFit']
    if install_lvl >= 2:
        IncFile['DLM_FemtoTools'] = ['CommonAnaFunctions']
        IncFile['CECA'] = ['CECA','TREPNI']

    # recreate install dir
    if os.path.exists(PATH_CATS):
        if len(os.listdir(PATH_CATS)) != 0:
            print(bcolors.WARNING+' The installation folder is not empty, it will be DELETED!'+bcolors.ENDC)
            print('Proceed? (y/n) ', end='')
            if not yes(): return
        os.system('rm -rf '+PATH_CATS)
    os.mkdir(PATH_CATS)
    os.mkdir(PATH_CATS_CMAKE)

    cmakefile = open(PATH_CATS_CMAKE+'/CMakeLists.txt', 'w')

    cmakefile.write('cmake_minimum_required(VERSION 3.16)\n')
    cmakefile.write('project(CATS3 LANGUAGES CXX)\n\n')
    cmakefile.write('set(DLM_REPO "'+PATH_DLM+'")\n')
    cmakefile.write('set(CATS_INSTALL "'+PATH_CATS+'")\n')
    cmakefile.write('set(GSL_INCLUDE "'+PATH_GSL_INC+'")\n')
    cmakefile.write('set(GSL_LIB "'+PATH_GSL_LIB+'")\n')
    if install_lvl >= 1:
        cmakefile.write('set(ROOT_PATH "'+PATH_ROOT+'")\n')
    if IS_MAC:
        cmakefile.write('set(OMP_INCLUDE "'+OMP_INC+'")\n')
        cmakefile.write('set(OMP_LIB "'+OMP_LIB+'")\n')
    cmakefile.write('\n')

    # include dirs
    cmakefile.write('include_directories(' + GSL_INCLUDE_PATHS + ')\n')
    for d in IncDir:
        cmakefile.write('include_directories(${DLM_REPO}/'+d+')\n')
    if IS_MAC:
        cmakefile.write('include_directories(${OMP_INCLUDE})\n')

    # library
    cmakefile.write('\nadd_library(CATS SHARED\n')
    for d in IncDir:
        for f in IncFile[d]:
            if f in OnlyHeader:
                continue
            cmakefile.write('  ${DLM_REPO}/'+d+'/'+f+'.cpp\n')
    cmakefile.write(')\n\n')

    # compile + link flags
    cmakefile.write('target_compile_features(CATS PRIVATE cxx_std_11)\n')
    cmakefile.write('target_compile_options(CATS PRIVATE -O2 -Wno-vla-cxx-extension)\n')

    if IS_MAC:
        cmakefile.write('target_compile_options(CATS PRIVATE -Xpreprocessor -fopenmp)\n')
        cmakefile.write('target_link_directories(CATS PRIVATE ${OMP_LIB} ${GSL_LIB})\n')
        cmakefile.write('target_link_libraries(CATS PRIVATE gsl gslcblas omp)\n')
    else:
        cmakefile.write('target_link_libraries(CATS PRIVATE -L${GSL_LIB} -lgsl -lgslcblas)\n')

    if install_lvl >= 1:
        # ROOT flags via root-config (macOS/brew output often contains newlines/whitespace)
        cmakefile.write('\n# ===== ROOT (from root-config) =====\n')
    
        cmakefile.write('execute_process(\n')
        cmakefile.write('  COMMAND bash -c "source ${ROOT_PATH}/bin/thisroot.sh && root-config --cflags"\n')
        cmakefile.write('  OUTPUT_VARIABLE ROOT_CXXFLAGS_RAW\n')
        cmakefile.write(')\n')
        cmakefile.write('string(STRIP "${ROOT_CXXFLAGS_RAW}" ROOT_CXXFLAGS_RAW)\n')
        cmakefile.write('separate_arguments(ROOT_CXXFLAGS NATIVE_COMMAND "${ROOT_CXXFLAGS_RAW}")\n')
        cmakefile.write('target_compile_options(CATS PRIVATE ${ROOT_CXXFLAGS})\n')
    
        cmakefile.write('execute_process(\n')
        cmakefile.write('  COMMAND bash -c "source ${ROOT_PATH}/bin/thisroot.sh && root-config --libs"\n')
        cmakefile.write('  OUTPUT_VARIABLE ROOT_LIBS_RAW\n')
        cmakefile.write(')\n')
        cmakefile.write('string(STRIP "${ROOT_LIBS_RAW}" ROOT_LIBS_RAW)\n')
        cmakefile.write('separate_arguments(ROOT_LIBS NATIVE_COMMAND "${ROOT_LIBS_RAW}")\n')
        cmakefile.write('target_link_libraries(CATS PRIVATE ${ROOT_LIBS})\n')

    cmakefile.write('\nset_target_properties(CATS PROPERTIES\n')
    cmakefile.write('  ARCHIVE_OUTPUT_DIRECTORY "${CATS_INSTALL}/lib"\n')
    cmakefile.write('  LIBRARY_OUTPUT_DIRECTORY "${CATS_INSTALL}/lib"\n')
    cmakefile.write('  RUNTIME_OUTPUT_DIRECTORY "${CATS_INSTALL}/bin"\n')
    cmakefile.write(')\n\n')

    # install headers
    for d in IncDir:
        for f in IncFile[d]:
            cmakefile.write('install(FILES ${DLM_REPO}/'+d+'/'+f+'.h DESTINATION ${CATS_INSTALL}/include)\n')
    cmakefile.write('install(PROGRAMS ${DLM_REPO}/bin/cats-config DESTINATION ${CATS_INSTALL}/bin)\n')

    cmakefile.close()

    print()
    print(bcolors.OKGREEN+'READY TO INSTALL!'+bcolors.ENDC)
    print(' Proceed? (y/n) ', end='')
    if not yes():
        return

    os.chdir(PATH_CATS_CMAKE)
    if os.system('cmake -S . -B .'):
        print(bcolors.FAIL+'"cmake" failed'+bcolors.ENDC); return
    if os.system('cmake --build . -j'+str(multiprocessing.cpu_count())):
        print(bcolors.FAIL+'"make" failed'+bcolors.ENDC); return
    if os.system('cmake --install .'):
        print(bcolors.FAIL+'"install" failed'+bcolors.ENDC); return

    print(bcolors.OKGREEN+'CATS was successfully installed!'+bcolors.ENDC)
    print(' Installed in: '+PATH_CATS)
    print(' Lib: '+PATH_CATS+'/lib')
    print(' Include: '+PATH_CATS+'/include')
    if IS_MAC:
        print(bcolors.WARNING+'Note: runtime may need DYLD_LIBRARY_PATH='+bcolors.ENDC+PATH_CATS+'/lib')

if __name__ == "__main__":
    quick_install(sys.argv)
