#!/bin/bash

DEBUG=false
INSTALL_TESTS=false
RUN_TESTS=false
while true; do
    case "$1" in
	-d | --debug ) DEBUG=true; shift ;;
	-t | --test  ) INSTALL_TESTS=true; shift ;;
	--run_test   ) RUN_TESTS=true; shift ;;
	* ) break ;;
    esac
done

run_test() {
    if $1; then
	../../unit_tests --log_level=message --show_progress=yes --color_output=yes
    fi
}

run_debug_test() {
    if $1; then
	test/unit_tests --log_level=message --show_progress=yes --color_output=yes
    fi
}

if $DEBUG; then
    echo "Building in debug mode..."
    mkdir -p build/debug/
    cd build/debug
    if $INSTALL_TESTS; then
	cmake -Wno-dev -UCMAKE_BUILD_TYPE -DCMAKE_BUILD_TYPE=Debug -DBUILD_TEST=True ../..
    else
	cmake -Wno-dev -UCMAKE_BUILD_TYPE -DCMAKE_BUILD_TYPE=Debug -DBUILD_TEST=False ../..
    fi
    make && run_debug_test $RUN_TESTS
    cd ../..
else
    echo "Building in release mode..."
    mkdir -p build/release
    cd build/release
    if $INSTALL_TESTS; then
	cmake -Wno-dev -UCMAKE_BUILD_TYPE -DCMAKE_BUILD_TYPE=Release -DBUILD_TEST=True ../..
    else
	cmake -Wno-dev -UCMAKE_BUILD_TYPE -DCMAKE_BUILD_TYPE=Release -DBUILD_TEST=False ../..
    fi
    make && make install && run_test $RUN_TESTS #only run tests if make passes
    cd ../..
fi
