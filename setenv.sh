#
# Modify with your own path to this directory
# Then source this file (Unix and Unix-like systems)
# 	$ source setenv.sh
#

install_folder=$HOME/path_to_Survival_directory/

## Linux Library Path:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$install_folder/ext_lib

## Mac Library Path:
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$install_folder/ext_lib

export DATA=$install_folder/data
export PATH=$PATH:$install_folder
