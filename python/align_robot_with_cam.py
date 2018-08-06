'''
1. detect start and finish of movement in joint space
2. check alignment with camera data 
camera 0: robot starts moving around 765 and stops around 7380
frame rate is around 180
for camera 1 it's the same
3. fix start stop frames for camera 0 and 1
4. for every camera frame in between run kinematics to get 3d position of racket
and also orientation
5. ball is x-cm behind robot centre (about 10cm?) along racket 
and a few cm on top
need to fix it?
6. run ransac to estimate one 3x4 matrix
7. predict positions for still balls and compare with old calibration
what is the predicted table length and width for instance?
