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
  if ( argc != 6 && argc !=7 ) 
    {
      printf( "usage: %s datafilename pairdistdatafile maxlag reconfigmatrixoutfile sepaveragedoutfile meanreconfigoutfile [snapshotfrequency]\n", argv[0] );
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
  if ( datafile != NULL )
    {
      char line [ 1000 ];
      for (i=0; i<sizeof line; i++)
	{
	  line[i] = NULL;
	}

      while ( fgets ( line, sizeof line, datafile ) != NULL )
	{
	  if ( line[0] == 't' )
	    {
	      totalsnapshots++;
	    }
	  size=0;
	  for (i=0; i<sizeof line; i++)
	    {
	      if ( line[i] == '.' )
		{
		  size++;
		}
	    }
	  //	  fputs ( line, stdout ); /* write the line */
	}
      fclose ( datafile );
    }

  int snapfreq = 1; // only read in every snapfreq snapshots
  if ( argc == 7 )
    {
      snapfreq = atoi(argv[6]);
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
  maxlag /= snapfreq;
  printf("Effective maxlag (maxlag/snapfreq): %d\n",maxlag);
  if (maxlag > snapshots)
    {
      printf("Cannot have a maxlag greater than the number of snapshots.\n");
      return;
    }
  double ***timeseries;
  timeseries = alloc_data3d(size,size,snapshots);
  datafile=fopen(argv[1], "r");
  for ( snapshot = 0; snapshot < totalsnapshots; snapshot++ )
    {  
      fscanf(datafile,"%s",timestepstr);
      fscanf(datafile,"%d",&timestepnum);
      //printf("timestep: %s, number: %d \n",timestepstr,timestepnum);
      for ( row = 0; row < size; row++)
	{
	  for ( column = 0; column < size; column++ )
	    {
	      fscanf(datafile,"%lf", &pairdist);
	      if( snapshot == 0 || (snapshot+1)%snapfreq == 0 )
		{
		  //printf("%d\n",snapshot);
		  timeseries[row][column][snapshot/snapfreq]=pairdist;
		}
	    }
	}
    }
  
    /* for ( snapshot = 0; snapshot < snapshots; snapshot++ ) */
    /* { */
    /*   for ( row = 0; row < size; row++) */
    /* 	{ */
    /* 	  for ( column = 0; column < size; column++ ) */
    /* 	    { */
    /* 	      printf("%f \n", timeseries[row][column][snapshot]); */
    /* 	    } */
    /* 	} */
    /* } */
  double **reconfigurationtimes;
  reconfigurationtimes = alloc_data2d(size,size);
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  reconfigurationtimes[row][column] = 0.0;
	}
    }

  double **averagevalues;
  averagevalues = alloc_data2d(size,size);
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  averagevalues[row][column] = 0.0;
	}
    }
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
  double **variancevalues;
  variancevalues = alloc_data2d(size,size);
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  variancevalues[row][column] = 0.0;
	}
    }
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
  FILE *reconfigmatoutfile;
  reconfigmatoutfile=fopen(argv[3], "w");
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  fprintf(reconfigmatoutfile,"%f ",reconfigurationtimes[row][column]*snapfreq);
	}
      fprintf(reconfigmatoutfile,"\n");
    }  
  fclose(reconfigmatoutfile);
  
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
