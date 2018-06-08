#!/bin/bash
# IMPORTANT NOTE: run in a subshell for cd commands: 
# . install.sh

usage() { echo "Argument missing for running unit tests. Usage: $0 [-t <1|0>] " 1>&2; exit 1; }

while getopts ":t:" opt; do
  case $opt in
    t)
	arg=${OPTARG}
	((arg == 1 || arg == 0)) || usage
        #echo "-t was triggered, Parameter: $OPTARG" >&2
        ;;
    \?)
	echo "Invalid option: -$OPTARG" >2
	usage
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
	usage
	exit 1
	;;
  esac
done

mkdir build/
mkdir build/debug
mkdir build/release
cd build/release
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release ../..
make
make install
cd ../..
#./unit_tests --show_progress=yes
#echo "arg = ${arg}"
if [ "$arg" == 1 ]; then
   ./unit_tests --log_level=message --show_progress=yes
fi

# FOR OLD SL
#cp SL_code/ilc_task.c ~/robolab/barrett/src/ilc_task.c

# FOR NEW SL COPY EVERYTHING INCLUDING CMAKELIST AND README
#cp SL_code/*.c ~/sl_xeno/barrett/src/learning_control/*.c
#cp SL_code/*.h ~/sl_xeno/barrett/include/*.h
#cp CMakeLists.txt ~/sl_xeno/barrett/src/learning_control/CMakeLists.txt
#cp readme.txt ~/sl_xeno/barrett/src/learning_control/readme.txt
