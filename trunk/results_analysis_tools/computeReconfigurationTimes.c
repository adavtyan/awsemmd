
/*# Created by Nick Schafer on 2/24/12
# This script takes as input a timeseries of pairwise distances
# and computes a matrix of autocorrelation functions and averages
# these autocorrelation functions over sequence separation
# to get reconfiguration times as a function of sequence separation.

# Features to add:
# * save data to files
# * flexibly choose what data to save
# * plot reconfiguration time versus sequence separation
# * plot autocorrelation functions
# * compute uncertainties for acfs and reconfiguration times
# * optionally specify number of snapshots to process
# * flexibly specify separation range (e.g. range(1,5,20))

# import necessary libraries
import numpy
import sys
import commands

# check commad line arguments
if len(sys.argv) != 3:
    print "Usage: python computeReconfigurationTimes.py maxlag"
    sys.exit()

# set system parameters
maxlag=int(sys.argv[2]) # maximum lag time at which to compute autocorrelation
print "maximum lag time to be computed: " + str(maxlag)

# open the data file
f=open(sys.argv[1],'r')

# read in system size
print "Reading system size..."
size=0
for line in f:
    line=line.split()
    if line[0] == 'timestep':
        continue
    
    for element in line:
        size += 1
    
    break

print "size:" + str(size)

# set min and max sequence separation
minseqsep=1 # at least 1
maxseqsep=size-1 # at most size-1

# read in number of snapshots
print "Reading in number of snapshots..."
snapshots=int(commands.getstatusoutput('grep timestep ' + sys.argv[1] + ' | wc -l')[1])
print "snapshots:" + str(snapshots)

if maxlag > snapshots:
    print "Cannot have a maxlag greater than the number of snapshots."
    sys.exit()

# close and reopen data file to read over again
f.close()
f=open(sys.argv[1],'r')

# for debugging purposes
# testrow=10 
# testcolumn=15

# build timeseries matrix
print "Building timeseries matrix..."
timeseries=[]
for row in range(size):
    timeseries.append([])
    for column in range(size):
        timeseries[row].append([])
        for snapshot in range(snapshots):
            timeseries[row][column].append([])

# read in timeseries values
print "Reading timeseries values..."
row=0
column=0
snapshot=-1

for line in f:
    line=line.split()
    if line[0] == 'timestep':
        snapshot=snapshot+1
        row=0
        continue
    
    for element in line:
        timeseries[row][column][snapshot]=element
        column=column+1

    row=row+1
    column=0

# print timeseries[testrow][testcolumn]

# build reconfiguration times matrix
print "Building reconfiguration times matrix..."
reconfigurationtimes=[]
for row in range(size):
    reconfigurationtimes.append([])
    for column in range(size):
        reconfigurationtimes[row].append(0.0)

# build average matrix
print "Building average matrix..."
averagevalues=[]
for row in range(size):
    averagevalues.append([])
    for column in range(size):
        averagevalues[row].append([])

# calculate average matrix
print "Calculating average matrix..."
for row in range(size):
    for column in range(row,size):
        average = 0
        for snapshot in range(snapshots):
            average += float(timeseries[row][column][snapshot])

        averagevalues[row][column] = average/float(snapshots)

# subtract all averages from all time series
print "Subtracting average values..."
for row in range(size):
    for column in range(row,size):
        for snapshot in range(snapshots):
            timeseries[row][column][snapshot] = float(timeseries[row][column][snapshot])-averagevalues[row][column]

# build variance matrix
print "Building variance matrix..."
variancevalues=[]
for row in range(size):
    variancevalues.append([])
    for column in range(size):
        variancevalues[row].append([])

# calculate variance matrix
print "Calculating variance matrix..."
for row in range(size):
    for column in range(row,size):
        variance = 0
        for snapshot in range(snapshots):
            variance += pow(timeseries[row][column][snapshot],2)

        variancevalues[row][column] = variance/(snapshots-1)

# build autocorrelation matrix
print "Building autocorrelation matrix..."
autocorrelationvalues=[]
for row in range(size):
    autocorrelationvalues.append([])
    for column in range(size):
        autocorrelationvalues[row].append([])
        for tau in range(maxlag):
            autocorrelationvalues[row][column].append(0.0)

# calculate autocorrelation matrix
print "Calculating autocorrelation matrix..."
for row in range(size):
    for column in range(row,size):
        if abs(row-column) > maxseqsep or abs(row-column) < minseqsep:
            continue

        for tau in range(maxlag):
            for t in range(snapshots-tau):
                autocorrelationvalues[row][column][tau] += timeseries[row][column][t]*timeseries[row][column][t+tau]

        autocorrelationvalues[row][column][tau] /= float(maxlag-tau)

for row in range(size):
    for column in range(row,size):
        if abs(row-column) > maxseqsep or abs(row-column) < minseqsep:
            continue

        for tau in range(maxlag):
            if variancevalues[row][column] == 0:
                print "Zero variance at row " + str(row) + ", column " + str(column)
                print "Exiting..."
                sys.exit()
            if tau == 0:
                autocorrelationvalues[row][column][tau]=1.0
                continue

            autocorrelationvalues[row][column][tau] /= snapshots*variancevalues[row][column]

# calculate reconfiguration time matrix
print "Calculating reconfiguration time matrix..."
for row in range(size):
    for column in range(row,size):
        if abs(row-column) > maxseqsep or abs(row-column) < minseqsep:
            continue

        for tau in range(maxlag):
            if autocorrelationvalues[row][column][tau] < 0:
                break
            
            reconfigurationtimes[row][column] += autocorrelationvalues[row][column][tau]

# build autocorrelation array averaged over separation
print "Averaging autocorrelation functions over sequence separations..."
sepaveragedautocorrarray=[]
for sep in range(size):
    sepaveragedautocorrarray.append([])
    for tau in range(maxlag):
        sepaveragedautocorrarray[sep].append(0.0)

# calculate autocorrelation array averaged over separation    
for tau in range(maxlag):
    for row in range(size):
        for column in range(row,size):
            sepaveragedautocorrarray[abs(row-column)][tau] += autocorrelationvalues[row][column][tau] 

for tau in range(maxlag):
    for sep in range(size):
        sepaveragedautocorrarray[sep][tau] /= float(size-sep)
            
# print sepaveragedautocorrarray[5]
# print timeseries[testrow][testcolumn]
# print averagevalues[testrow][testcolumn]
# print variancevalues[testrow][testcolumn]
# print autocorrelationvalues[testrow][testcolumn]

# calculate reconfiguration times for all seps between minsep and maxsep
# (sum all values of averaged autocor function before it goes to zero)

# build decay time array
print "Building decay time arrays..."
decaytimevalues=[]
for sep in range(0,size):
    decaytimevalues.append(0.0)

# calculate decay times
reachedzero=0
print "Calculating decay times..."
for sep in range(minseqsep,maxseqsep+1):
    decaytime=0.0
    for tau in range(maxlag):
        if sepaveragedautocorrarray[sep][tau] < 0:
            reachedzero=1
            break

        decaytime += sepaveragedautocorrarray[sep][tau]
    
    decaytimevalues[sep] = decaytime

print decaytimevalues
print reconfigurationtimes

if reachedzero == 0:
    print "Warning: maxlag may be too small"
*/
          
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
  if ( argc != 6 ) 
    {
      printf( "usage: %s datafilename pairdistdatafile maxlag reconfigmatrixoutfile sepaveragedoutfile meanreconfigoutfile\n", argv[0] );
      exit(1);
    }
  int maxlag=atoi(argv[2]);
  printf("maxlag to be computed: %d\n",maxlag);
  FILE *datafile;
  datafile=fopen(argv[1], "r");
  printf("Reading system size and number of snapshots...\n");
  int size=0;
  int snapshots=0;
  int i=0;
  if ( datafile != NULL )
    {
      char line [ 1000 ]; /* or other suitable maximum line size */
      for (i=0; i<sizeof line; i++)
	{
	  line[i] = NULL;
	}

      while ( fgets ( line, sizeof line, datafile ) != NULL ) /* read a line */
	{
	  if ( line[0] == 't' )
	    {
	      snapshots++;
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
  printf("System size: %d\n",size);
  printf("Number of snapshots: %d\n",snapshots);
  int minseqsep = 1;
  int maxseqsep = size-1;
  if (maxlag > snapshots)
    {
      printf("Cannot have a maxlag greater than the number of snapshots.\n");
      return;
    }
  double ***timeseries;
  timeseries = alloc_data3d(size,size,snapshots);
  int row,column,snapshot;
  char timestepstr[8];
  int timestepnum;
  double pairdist;
  datafile=fopen(argv[1], "r");
  for ( snapshot = 0; snapshot < snapshots; snapshot++ )
    {  
      fscanf(datafile,"%s",timestepstr);
      fscanf(datafile,"%d",&timestepnum);
      //printf("timestep: %s, number: %d \n",timestepstr,timestepnum);
      for ( row = 0; row < size; row++)
	{
	  for ( column = 0; column < size; column++ )
	    {
	      fscanf(datafile,"%lf", &pairdist);
	      timeseries[row][column][snapshot]=pairdist;
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
	      printf("WARNING: reconfigurationtimes[%d][%d] is likely underestimated. The corresponding autocorrelation function didn't cross zero within maxlag (%d) time units. \n",row,column,maxlag);
	    }
	}
    }  
  FILE *reconfigmatoutfile;
  reconfigmatoutfile=fopen(argv[3], "w");
  for ( row = 0; row < size; row++)
    {
      for ( column = 0; column < size; column++ )
	{
	  fprintf(reconfigmatoutfile,"%f ",reconfigurationtimes[row][column]);
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
      fprintf(sepaveragedoutfile,"%f \n",sepaveragedreconfigtimes[sep]);
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
  fprintf(meanreconfigoutfile,"%f ",meanreconfigtime);
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
