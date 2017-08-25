Do the following: 

1. initUserGraphics.c needs to have the display_ball_gun function!
2. initUserGraphics.c: add the function!
3. include the header files in the include/path 
(since initUserGraphics.c needs then the table.h header file!)
4. include the following line in initUserGraphics.c:
addToUserGraphics("table_tennis", "Display Table, 
Ball and Racket", &(display_table_tennis), 14*sizeof(double));
This was necessary as the simulation had problems with the 'standard' new-SL way
of writing graphics/sim.
