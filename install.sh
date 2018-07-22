#!/bin/bash

DEBUG=false
TEST=false
TEST_CMD="*" #--run_test={KIN,KF,OPT,TT,SL,SERVE}
BUILD="mpi/laptop"
while true; do
    case "$1" in
	-b | --build ) BUILD="$2"; shift 2 ;;
	-d | --debug ) DEBUG=true; shift ;;
	-t | --test  ) TEST=true; shift ;;
	--run_test   ) TEST_CMD="$2"; shift ;;
	* ) break ;;
    esac
done

if [ "$BUILD" = "robot" ]; then
    echo "Loading robot cmake files..."
    cp cmake_files/cmakelists_robot/CMakeLists.txt .
    cp cmake_files/cmakelists_robot/src/CMakeLists.txt src/
    cp cmake_files/cmakelists_robot/test/CMakeLists.txt test/
else
    cp cmake_files/cmakelists_mpi/CMakeLists.txt .
    cp cmake_files/cmakelists_mpi/src/CMakeLists.txt src/
    cp cmake_files/cmakelists_mpi/test/CMakeLists.txt test/
fi

run_test() {
    if $1; then
	../../unit_tests --log_level=message --show_progress=yes --run_test="$2" --color_output=yes
    fi
}

run_debug_test() {
    if $1; then
	test/unit_tests --log_level=message --show_progress=yes --run_test="$2" --color_output=yes
    fi
}

if $DEBUG; then
    echo "Building in debug mode..."
    mkdir -p build/debug/
    cd build/debug
    if $TEST; then
	cmake -Wno-dev -UCMAKE_BUILD_TYPE -DCMAKE_BUILD_TYPE=Debug -DBUILD_TEST=True ../..
    else
	cmake -Wno-dev -UCMAKE_BUILD_TYPE -DCMAKE_BUILD_TYPE=Debug -DBUILD_TEST=False ../..
    fi
    make && run_debug_test $TEST $TEST_CMD
#    if $TEST; then
#	test/unit_tests --log_level=message --show_progress=yes --run_test="$TEST_CMD" --color_output=yes
#    fi
    cd ../..
else
    echo "Building in release mode..."
    mkdir -p build/release
    cd build/release
    if $TEST; then
	cmake -Wno-dev -UCMAKE_BUILD_TYPE -DCMAKE_BUILD_TYPE=Release -DBUILD_TEST=True ../..
    else
	cmake -Wno-dev -UCMAKE_BUILD_TYPE -DCMAKE_BUILD_TYPE=Release -DBUILD_TEST=False ../..
    fi
    make && make install && run_test $TEST $TEST_CMD #only run tests if make passes
    cd ../..
#    if $TEST; then
#	./unit_tests --log_level=message --show_progress=yes --run_test="$TEST_CMD" --color_output=yes
#    fi
fi

# FOR OLD SL
#cp SL_code/ilc_task.c ~/robolab/barrett/src/ilc_task.c

# FOR NEW SL COPY EVERYTHING INCLUDING CMAKELIST AND README
#cp SL_code/*.c ~/sl_xeno/barrett/src/learning_control/*.c
#cp SL_code/*.h ~/sl_xeno/barrett/include/*.h
#cp CMakeLists.txt ~/sl_xeno/barrett/src/learning_control/CMakeLists.txt
#cp readme.txt ~/sl_xeno/barrett/src/learning_control/readme.txt
