# copy code to both SL versions

# first to usual SL
cp SL_code/serve/CMakeLists.txt ~/projects/sl_xeno/barrett/src/serve/
cp SL_code/serve/serve_task.c ~/projects/sl_xeno/barrett/src/serve/
cp SL_code/grav_comp/grav_comp.c ~/projects/sl_xeno/barrett/src/grav_comp/
cp SL_code/grav_comp/CMakeLists.txt ~/projects/sl_xeno/barrett/src/grav_comp/
cp SL_code/player/table_tennis_task.c ~/projects/sl_xeno/barrett/src/table-tennis/

# then to BORIS SL with custom inverse dynamics
cp SL_code/serve/CMakeLists.txt ~/sebastian/sl_xeno/sl_wam/src/serve/
cp SL_code/serve/serve_task.c ~/sebastian/sl_xeno/sl_wam/src/serve/
cp SL_code/grav_comp/grav_comp.c ~/sebastian/sl_xeno/sl_wam/src/grav_comp/
cp SL_code/grav_comp/CMakeLists.txt ~/sebastian/sl_xeno/sl_wam/src/grav_comp/
cp SL_code/player/table_tennis_task.c ~/sebastian/sl_xeno/sl_wam/src/table-tennis/
