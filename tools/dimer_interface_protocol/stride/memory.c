#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <ctype.h>

void *ckalloc(size_t bytes)
{
  register void *ret;
  void die(char *format, ... );
  
  if( !(ret = malloc(bytes)) ) die("Out of  memory\n");

  return ret;	
}

float **FloatMatrix(int M, int N)
{
  int m;
  float **Matrix;

  Matrix = (float **)ckalloc(M*sizeof(float *));

  for( m=0; m<M; m++ ) Matrix[m] = (float *)ckalloc(N*sizeof(float));

  return(Matrix);
}

float ***FloatCube(int M, int N, int K)
{
  int m, n, k;
  float ***Cube;

  Cube = (float ***)ckalloc(M*sizeof(float **));

  for( m=0; m<M; m++ ) {
    Cube[m] = (float **)ckalloc(N*sizeof(float *));
    for( n=0; n<N; n++ )
      Cube[m][n] = (float *)ckalloc(K*sizeof(float));
  }

  for( m=0; m<M; m++ )
    for( n=0; n<N; n++ )
      for( k=0; k<K; k++ )
	Cube[m][n][k] = 0.0;

  return(Cube);
}

float ****Float4Dim(int M, int N, int K, int L)
{
  int m, n, k, l;
  float ****FourDim;

  FourDim = (float ****)ckalloc(M*sizeof(float ***));


  for( m=0; m<M; m++ ) {
    FourDim[m] = (float ***)ckalloc(N*sizeof(float **));
    for( n=0; n<N; n++ ) {
      FourDim[m][n] = (float **)ckalloc(K*sizeof(float*));
      for( k=0; k<K; k++ )
	FourDim[m][n][k] = (float *)ckalloc(L*sizeof(float));
    }
  }

  for( m=0; m<M; m++ )
    for( n=0; n<N; n++ )
      for( k=0; k<K; k++ )
	for( l=0; l<L; l++ )
	  FourDim[m][n][k][l] = 0.0;

  return(FourDim);
}

void FreeFloatMatrix(float **Matrix, int M)
{
  int m;

  for( m=0; m<M; m++ ) free(Matrix[m]);

  free(Matrix);

}

int ***IntCube(int M, int N, int K)
{
  int m, n, k;
  int ***Cube;

  Cube = (int ***)ckalloc(M*sizeof(int **));

  for( m=0; m<M; m++ ) {
    Cube[m] = (int **)ckalloc(N*sizeof(int *));
    for( n=0; n<N; n++ ) Cube[m][n] = (int *)ckalloc(K*sizeof(int));
  }
  
  for( m=0; m<M; m++ )
    for( n=0; n<N; n++ )
      for( k=0; k<K; k++ )
	Cube[m][n][k] = 0;

  return(Cube);
}


int **IntMatrix(int M, int N)
{
  int m;
  int **Matrix;

  Matrix = (int **)ckalloc(M*sizeof(int *));

  for( m=0; m<M; m++ ) Matrix[m] = (int *)ckalloc(N*sizeof(int));

  return(Matrix);
}

int ****Int4Dim(int M, int N, int K, int L)
{
  int m, n, k, l;
  int ****FourDim;

  FourDim = (int ****)ckalloc(M*sizeof(int ***));


  for( m=0; m<M; m++ ) {
    FourDim[m] = (int ***)ckalloc(N*sizeof(int **));
    for( n=0; n<N; n++ ) {
      FourDim[m][n] = (int **)ckalloc(K*sizeof(int*));
      for( k=0; k<K; k++ )
	FourDim[m][n][k] = (int *)ckalloc(L*sizeof(int));
    }
  }

  for( m=0; m<M; m++ )
    for( n=0; n<N; n++ )
      for( k=0; k<K; k++ )
	for( l=0; l<L; l++ )
	  FourDim[m][n][k][l] = 0;

  return(FourDim);
}

void FreeIntMatrix(int **Matrix, int M)
{
  int m;

  for( m=0; m<M; m++ ) free(Matrix[m]);

  free(Matrix);

}

char **CharMatrix(int M, int N)
{
  int m;
  char **Matrix;

  Matrix = (char **)ckalloc(M*sizeof(char *));

  for( m=0; m<M; m++ ) Matrix[m] = (char *)ckalloc(N*sizeof(char));

  return(Matrix);
}

void FreeCharMatrix(char **Matrix, int M)
{
  int m;

  for( m=0; m<M; m++ ) free(Matrix[m]);

  free(Matrix);

}

void FreeIntCube(int ***Cube, int M, int N)
{
  int m, n;

  for( m=0; m<M; m++ ) {
    for( n=0; n<N; n++ )
      free(Cube[m][n]);
    free(Cube[m]);
  }

  free(Cube);
}

void FreeFloatCube(float ***Cube, int M, int N)
{
  int m, n;

  for( m=0; m<M; m++ ) {
    for( n=0; n<N; n++ )
      free(Cube[m][n]);
    free(Cube[m]);
  }

  free(Cube);
}


