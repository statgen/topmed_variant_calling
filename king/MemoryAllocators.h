#ifndef __MEMORY_ALLOCATORS_H__
#define __MEMORY_ALLOCATORS_H__

char **  AllocateCharMatrix(int rows, int cols);
void     FreeCharMatrix(char ** & matrix, int rows);

float ** AllocateFloatMatrix(int rows, int cols);
void     FreeFloatMatrix(float ** & matrix, int rows);

int  **  AllocateIntMatrix(int rows, int cols);
void     FreeIntMatrix(int ** & matrix, int rows);

char *** AllocateCharCube(int n, int rows, int cols);
void     FreeCharCube(char *** & matrix, int n, int rows);

#endif

