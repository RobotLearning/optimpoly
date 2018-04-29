# IMPORTANT NOTE: run in a subshell for cd commands: 
# . install.sh

mkdir build/
mkdir build/debug
mkdir build/release
cd build/release
cmake -DCMAKE_BUILD_TYPE=Release ../..
make
make install
cd ../..
./unit_tests --log_level=message --show_progress=yes

# FOR OLD SL
#cp SL_code/ilc_task.c ~/robolab/barrett/src/ilc_task.c

# FOR NEW SL COPY EVERYTHING INCLUDING CMAKELIST AND README
#cp SL_code/*.c ~/sl_xeno/barrett/src/learning_control/*.c
#cp SL_code/*.h ~/sl_xeno/barrett/include/*.h
#cp CMakeLists.txt ~/sl_xeno/barrett/src/learning_control/CMakeLists.txt
#cp readme.txt ~/sl_xeno/barrett/src/learning_control/readme.txt
