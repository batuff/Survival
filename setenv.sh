#
# Modify with your own path to this directory
# Then source this file (Unix and Unix-like systems)
# 	$ source setenv.sh
#

install_folder=$HOME/path_to_this_directory/

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$install_folder/lib
export DATA=$install_folder/data
export PATH=$PATH:$install_folder
