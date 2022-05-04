# halftime-show-simulation
Developed a half-time show using UAVs. This football half-time show 3D simulation uses MPI and OpenGL.




How to compile:
>module load mesa gcc mvapich2
>mpic++ *.cpp –lGLU –lglut –std=c++11

To run our code you will need to create an interactive session with 16 processes either using a remote VNC setup or using an X11 server.

In the folder where your executable is located, execute the following commands:
>module load mesa gcc mvapich2
>mpirun -np 16 ./name_of_your_exectable

In a few seconds, a window should appear displaying your OpenGL rendering.
