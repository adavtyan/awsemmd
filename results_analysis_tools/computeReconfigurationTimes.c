/* Created by Nick Schafer on 3/15/12
This program takes as input a timeseries of pairwise distances
and computes a matrix of autocorrelation functions, computes
a matrix of reconfiguration times and averages
these reconfiguration times over sequence separation
to get reconfiguration times as a function of sequence separation. 

This was adapted from a python script of the same name that was
also written by Nick Schafer.
*/
      
#include <stdio.h> 

// functions for dynamically allocating 2d and 3d arrays
double ***alloc_data2d(size_t xlen, size_t ylen);   
void free_data2d(double **data, size_t xlen);
double ***alloc_data3d(size_t xlen, size_t ylen, size_t zlen);   
void free_data3d(double ***data, size_t xlen, size_t ylen);

// main program
int main(int argc, char *argv[])
{
  if ( argc != 8 && argc !=9 ) 
    {
      printf( "usage: %s datafilename pairdistdatafile maxlag reconfigmatrixoutfile sepaveragedoutfile meanreconfigoutfile meanvaluesfile variancevaluesfile [snapshotfrequency]\n", argv[0] );
      exit(1);
    }
  int maxlag=atoi(argv[2]);
  printf("maxlag to be computed: %d\n",maxlag);
  FILE *datafile;
  datafile=fopen(argv[1], "r");
  printf("Reading system size and number of snapshots...\n");
  int size=0;
  int totalsnapshots=0;
  int i=0;
  int lastwasspace=1;
  if ( datafile != NULL )
    {
      char line [ 1000 ];
      for (i=0; i<sizeof line; i++)
	{
	  line[i] = ' ';
	}

      while ( fgets ( line, sizeof line, datafile ) != NULL )
	{
	  if ( line[0] == 't' )
	    {
	      totalsnapshots++;
	    }
	}
      size=0;
      for (i=0; i<sizeof line; i++)
	{
	  if ( line[i] == ' ' )
	    {
	      lastwasspace=1;
	    }
	  else if ( lastwasspace == 1 )
	    {
	      lastwasspace=0;
	      if ( line[i] != '\n' )
		{
		  size++;
		}
	    }
	}

      fclose ( datafile );
    }
  
  int snapfreq = 1; // only read in every snapfreq snapshots
  if ( argc == 9 )
    {
      snapfreq = atoi(argv[8]);
    }
  printf("System size: %d\n",size);
  printf("Number of snapshots: %d\n",totalsnapshots);
  int minseqsep = 1;
  int maxseqsep = size-1;
  int row,column,snapshot;
  char timestepstr[8];
  int timestepnum;
  double pairdist;
  int snapshots;
  printf("Frequency at which to accept snapshots (snapfreq): %d\n",snapfreq);
  snapshots = totalsnapshots/snapfreq;
  printf("Number of snapshots to process (snapshots/snapfreq): %d\n",snapshots);
  if (maxlag < snapfreq)
    {
      printf("Cannot have a maxlag less than the snapshot frequency.\n");
      return;
    }
  maxlag /= snapfreq;
  printf("Effective maxlag (maxlag/snapfreq): %d\n",maxlag);
  if (maxlag > snapshots)
    {
      printf("Cannot have a maxlag greater than the number of snapshots.\n");
      return;
    }

  printf("Allocating and reading time series data...\n");
  double ***timeseries;
  timeseries = alloc_data3d(size,size,snapshots);
  datafile=fopen(argv[1], "r");
  for ( snapshot = 0; snapshot < totalsnapshots; snapshot++ )
    {  
      fscanf(datafile,"%s",timestepstr);
      fscanf(datafile,"%d",&timestepnum);
      for ( row = 0; row < size; row++)
	{
	  for ( column = 0; column < size; column++ )
	    {
	      fscanf(datafile,"%lf", &pairdist);
	      if( snapshot == 0 || (snapshot+1)%snapfreq == 0 )
		{
		  timeseries[row][column][snapshot/snapfreq]=pairdist;
		}
	    }
	}
    }
  
  printf("Allocating and initializing reconfiguration time matrix...\n");
  double **reconfigurationtimes;
  reconfigurationtimes = alloc_data2d(size,size);
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  reconfigurationtimes[row][column] = 0.0;
	}
    }

  printf("Allocating and initializing average values matrix...\n");
  double **averagevalues;
  averagevalues = alloc_data2d(size,size);
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  averagevalues[row][column] = 0.0;
	}
    }

  printf("Calculating average values matrix...\n");
  double average = 0.0;
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  average = 0.0;
	  for ( snapshot = 0; snapshot < snapshots; snapshot++ )
	    { 	      
	      average += timeseries[row][column][snapshot];
	    }
	  averagevalues[row][column] = average/(double)snapshots;
	}
    }

  printf("Subtracting average values from time series data...\n");
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  average = 0.0;
	  for ( snapshot = 0; snapshot < snapshots; snapshot++ )
	    { 	      
	      timeseries[row][column][snapshot] -= averagevalues[row][column];
	    }
	}
    }  

  printf("Allocating and initializing variance values matrix...\n");
  double **variancevalues;
  variancevalues = alloc_data2d(size,size);
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  variancevalues[row][column] = 0.0;
	}
    }

  printf("Calculating variance values matrix...\n");
  double variance = 0.0;
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  variance = 0.0;
	  for ( snapshot = 0; snapshot < snapshots; snapshot++ )
	    { 	      
	      variance += pow(timeseries[row][column][snapshot],2);
	    }
	  variancevalues[row][column] = variance/((double)snapshots-1);
	}
    }  

  printf("Allocating and calculating autocorrelation values...\n");
  double ***autocorrelationvalues;
  autocorrelationvalues = alloc_data3d(size,size,maxlag);
  int tau = 0;
  int t = 0;
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  if ( abs(row-column) > maxseqsep || abs(row-column) < minseqsep )
	    {
	      continue;
	    }
	  for ( tau = 0; tau < maxlag; tau++ )
	    { 	  
	      for ( t = 0; t < snapshots-tau; t++ )
		{
		  autocorrelationvalues[row][column][tau] += timeseries[row][column][t]*timeseries[row][column][t+tau];
		}
	    }
	  autocorrelationvalues[row][column][tau] /= (double)(maxlag-tau);
	}
    }  

  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  if ( abs(row-column) > maxseqsep || abs(row-column) < minseqsep )
	    {
	      continue;
	    }
	  if ( variancevalues[row][column] == 0.0 )
	    {
	      printf("WARNING: variancevalues[%d][%d] is zero. The corresponding autocorrelation function has been artificially set to zero. \n",row,column);
	      for ( tau =0; tau < maxlag; tau++ )
		{
		  autocorrelationvalues[row][column][tau] = 0;
		}
	    }
	  else
	    {
	      for ( tau = 0; tau < maxlag; tau++ )
		{ 	  
		  if ( tau == 0 )
		    {
		      autocorrelationvalues[row][column][tau] = 1.0;
		      continue;
		    }
		  autocorrelationvalues[row][column][tau] /= snapshots*variancevalues[row][column];
		}	      
	    }
	}
    }  

  printf("Calculating reconfiguration times...\n");  
  int reachedzero = 0;
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  reachedzero = 0;
	  if ( abs(row-column) > maxseqsep || abs(row-column) < minseqsep )
	    {
	      continue;
	    }
	  for ( tau = 0; tau < maxlag; tau++ )
	    { 	  
	      if ( autocorrelationvalues[row][column][tau] < 0 )
		{
		  reachedzero = 1;
		  break;
		}
	      reconfigurationtimes[row][column] += autocorrelationvalues[row][column][tau];
	    }
	  if ( !reachedzero )
	    {
	      printf("WARNING: reconfigurationtimes[%d][%d] is likely underestimated. The corresponding autocorrelation function didn't cross zero within the effective maxlag (%d) time units. \n",row,column,maxlag);
	    }
	}
    }  

  printf("Writing output files...\n");  
  FILE *reconfigmatoutfile;
  reconfigmatoutfile=fopen(argv[3], "w");
  printf(argv[3]);
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  fprintf(reconfigmatoutfile,"%f ",reconfigurationtimes[row][column]*snapfreq);
	}
      fprintf(reconfigmatoutfile,"\n");
    }  
  fclose(reconfigmatoutfile);

  FILE *averageoutfile;
  averageoutfile=fopen(argv[6], "w");
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  fprintf(averageoutfile,"%f ",averagevalues[row][column]);
	}
      fprintf(averageoutfile,"\n");
    }  
  fclose(averageoutfile);

  FILE *varianceoutfile;
  varianceoutfile=fopen(argv[7], "w");
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  fprintf(varianceoutfile,"%f ",variancevalues[row][column]);
	}
      fprintf(varianceoutfile,"\n");
    }  
  fclose(varianceoutfile);
 
  double *sepaveragedreconfigtimes = malloc(size * sizeof *sepaveragedreconfigtimes);
  for ( i = 0; i < size; i++ )
    {
      sepaveragedreconfigtimes[i] = 0.0;
    }
  for ( row = 0; row < size; row++)
    {
      for ( column = row; column < size; column++ )
	{
	  sepaveragedreconfigtimes[abs(row-column)] += reconfigurationtimes[row][column];
	}
    }  
  int sep = 0;
  FILE *sepaveragedoutfile;
  sepaveragedoutfile=fopen(argv[4], "w");
  for ( sep = 0; sep < size; sep++ )
    {
      sepaveragedreconfigtimes[sep] /= (double)(size-sep);
      fprintf(sepaveragedoutfile,"%f \n",sepaveragedreconfigtimes[sep]*snapfreq);
    }
  fclose(sepaveragedoutfile);

  double meanreconfigtime = 0.0;
  for ( row = 0; row < size; row++)
    {
      for ( column = row; column < size; column++ )
	{
	  meanreconfigtime += reconfigurationtimes[row][column];
	}
    }
  
  meanreconfigtime /= ((double)(size*(size-1)))/2.0;
  FILE *meanreconfigoutfile;
  meanreconfigoutfile=fopen(argv[5], "w");
  fprintf(meanreconfigoutfile,"%f ",meanreconfigtime*snapfreq);
  fclose(meanreconfigoutfile);

  return 0;
}

double ***alloc_data3d(size_t xlen, size_t ylen, size_t zlen)
{
    double ***p;
    size_t i, j;

    if ((p = malloc(xlen * sizeof *p)) == NULL) {
        perror("malloc 1");
        return NULL;
    }

    for (i=0; i < xlen; ++i)
        p[i] = NULL;

    for (i=0; i < xlen; ++i)
        if ((p[i] = malloc(ylen * sizeof *p[i])) == NULL) {
            perror("malloc 2");
            free_data3d(p, xlen, ylen);
            return NULL;
        }

    for (i=0; i < xlen; ++i)
        for (j=0; j < ylen; ++j)
            p[i][j] = NULL;

    for (i=0; i < xlen; ++i)
        for (j=0; j < ylen; ++j)
            if ((p[i][j] = malloc(zlen * sizeof *p[i][j])) == NULL) {
                perror("malloc 3");
                free_data3d(p, xlen, ylen);
                return NULL;
            }

    return p;
}

void free_data3d(double ***data, size_t xlen, size_t ylen)
{
    size_t i, j;

    for (i=0; i < xlen; ++i) {
        if (data[i] != NULL) {
            for (j=0; j < ylen; ++j)
                free(data[i][j]);
            free(data[i]);
        }
    }
    free(data);
}

double ***alloc_data2d(size_t xlen, size_t ylen)
{
    double **p;
    size_t i;

    if ((p = malloc(xlen * sizeof *p)) == NULL) {
        perror("malloc 1");
        return NULL;
    }

    for (i=0; i < xlen; ++i)
        p[i] = NULL;

    for (i=0; i < xlen; ++i)
        if ((p[i] = malloc(ylen * sizeof *p[i])) == NULL) {
            perror("malloc 2");
            free_data2d(p, xlen);
            return NULL;
        }

    return p;
}

void free_data2d(double **data, size_t xlen)
{
    size_t i;

    for (i=0; i < xlen; ++i) {
        if (data[i] != NULL) {
	  free(data[i]);
        }
    }
    free(data);
}
