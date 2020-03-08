#ifndef _SEDONA_H
#define _SEDONA_H

// define real to choose either double or float precision
//typedef float real;
typedef double real;
typedef double  OpacityType;

#define DEFAULT_PARAM_FILE_NAME "param.lua"
#define MPI_PARALLEL 0

#define Max_MPI_Blocksize 1000000
#define TIME_PLOTFILE_NAME "integrated_quantities.dat"
#endif
