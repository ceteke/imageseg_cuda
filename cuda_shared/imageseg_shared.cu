#include <stdio.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <math.h>
#include <sys/time.h>
#include <iostream>
#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include <omp.h>
using namespace std;

#define PIXEL_WIDTH 3
#define TRESHOLD 4
#define BLOCK_SIZE 16

#define SIZE    (1024*1024*1024)
#define ELEMENTS    (SIZE / sizeof(unsigned int))
#define HASH_ENTRIES     1024

struct Entry {
	unsigned int    key;
	int            value;
	Entry           *next;
};

struct Table {
	size_t  count;
	Entry   **entries;
	Entry   *pool;
	Entry   *firstFree;
};

size_t hash( unsigned int key, size_t count ) {
	return key % count;
}

void initialize_table( Table &table, int entries, int elements ) {
	table.count = entries;
	table.entries = (Entry**)calloc( entries, sizeof(Entry*) );
	table.pool = (Entry*)malloc( elements * sizeof( Entry ) );
	table.firstFree = table.pool;
}

void free_table( Table &table ) {
	free( table.entries );
	free( table.pool );
}

void add_to_table( Table &table, unsigned int key, int value ) {
	size_t hashValue = hash( key, table.count );
	Entry *location = table.firstFree++;
	location->key = key;
	location->value = value;
	location->next = table.entries[hashValue];
	table.entries[hashValue] = location;
}

int getFromTable( Table &table, unsigned int key) {
	size_t hashValue = hash( key, table.count );

	Entry  *current = table.entries[hashValue];
	while (current != NULL) {
		if (hash( current->key, table.count ) != hashValue){
			current = current->next;
		} else {
			return current->value;
		}
	}
	return 0;
}

static const double kMicro = 1.0e-6;
double getTime()
{
	struct timeval TV;
	struct timezone TZ;

	const int RC = gettimeofday(&TV, &TZ);
	if(RC == -1) {
		printf("ERROR: Bad call to gettimeofday\n");
		return(-1);
	}

	return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}

void doHashTableStuff(Table red, Table green, Table blue, int x, int y, int *labels, unsigned char *data, int seed, int num_threads);
void cmdLine(int argc, char *argv[], char* i, char* output, int& d, int& num_threads, int &s);

__global__ void devicePhase(int *d_labels, unsigned char *d_data, int x, int y){
	int i = (blockDim.y-2) * blockIdx.y + threadIdx.y;
	int j = (blockDim.x-2) * blockIdx.x + threadIdx.x;

	if (i >= y || j >= x)return;

	int s_i = threadIdx.y;
	int s_j = threadIdx.x;

	int s_i3 = s_i*PIXEL_WIDTH;
	int s_j3 = s_j*PIXEL_WIDTH;

	int idx = i * x + j;
	int idx3 = idx*PIXEL_WIDTH;

	__shared__ int grid[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ int s_data[BLOCK_SIZE*PIXEL_WIDTH][BLOCK_SIZE];
	grid[s_i][s_j] = d_labels[idx];
	s_data[s_i3][s_j3] = (int) d_data[idx3];
	__syncthreads();

	int ll = grid[s_i][s_j];
	if (ll != 0){
		int dat = s_data[s_i3][s_j3];

		if(s_i != BLOCK_SIZE-1 && i != y-1 && abs(s_data[s_i3+PIXEL_WIDTH][s_j3] - dat) < TRESHOLD)
			grid[s_i][s_j] = max(grid[s_i][s_j], grid[s_i+1][s_j]);

		if(s_i != 0 && i != 0 && abs(s_data[s_i3-PIXEL_WIDTH][s_j3] - dat) < TRESHOLD)
			grid[s_i][s_j]  = max(grid[s_i][s_j] , grid[s_i-1][s_j]);

		if(s_i != BLOCK_SIZE-1 && s_j != BLOCK_SIZE-1 && i != y - 1 && j != x - 1 && abs(s_data[s_i3+PIXEL_WIDTH][s_j3+PIXEL_WIDTH] - dat) < TRESHOLD)
			grid[s_i][s_j]  = max(grid[s_i][s_j] , grid[s_i+1][s_j+1]);

		if(s_i != 0 && s_j != BLOCK_SIZE-1 && i != 0 && j != x - 1 && abs(s_data[s_i3-PIXEL_WIDTH][s_j3+PIXEL_WIDTH] - dat) < TRESHOLD)
			grid[s_i][s_j]  = max(grid[s_i][s_j] , grid[s_i-1][s_j+1]);

		if(s_i != BLOCK_SIZE-1 && s_j!= 0 && i != y - 1 && j != 0 && abs(s_data[s_i3+PIXEL_WIDTH][s_j3-PIXEL_WIDTH] - dat) < TRESHOLD)
			grid[s_i][s_j]  = max(grid[s_i][s_j] , grid[s_i+1][s_j-1]);

		if(s_i != 0 && s_j != 0 && i != 0 && j != 0 && abs(s_data[s_i3-PIXEL_WIDTH][s_j3-PIXEL_WIDTH] - dat) < TRESHOLD)
			grid[s_i][s_j]  = max(grid[s_i][s_j] , grid[s_i-1][s_j-1]);

		if (s_j != 0 && j != 0 && abs(s_data[s_i3][s_j3-PIXEL_WIDTH] - dat) < TRESHOLD)
			grid[s_i][s_j]  = max(grid[s_i][s_j] , grid[s_i][s_j-1]);

		if (s_j != BLOCK_SIZE-1 && j != x-1 && abs(s_data[s_i3][s_j3+PIXEL_WIDTH] - dat) < TRESHOLD)
			grid[s_i][s_j]  = max(grid[s_i][s_j] , grid[s_i][s_j+1]);
	}

	int label = grid[s_i][s_j];

	if (ll < label) {
		if (d_labels[ll - 1] < label)
			d_labels[ll - 1] = label;
	}

	if (label != 0) {
		d_labels[idx] = max(label, d_labels[label - 1]);
	}
}

int main(int argc,char **argv)
{
	int display=0;
	int num_threads=1;
	int seed=time(NULL);
	char *i = (char *)malloc(100*sizeof(char));
	char *output = (char *)malloc(100*sizeof(char));
	i = strcpy(i, "input.png");
	output = strcpy(output, "output.png");

	cmdLine( argc, argv, i, output, display, num_threads,seed);
	int x,y,n;

	printf("Reading image %s...\n", i);
	unsigned char *data = stbi_load(i, &x, &y, &n, 0);
	if (!data) {
		fprintf(stderr, "Couldn't load image.\n");
		return (-1);
	}

	int *labels = (int *)malloc(sizeof(int)*x*y);

	int *d_labels;
	unsigned char *d_data;

	cudaMalloc((void**)&d_labels, sizeof(int)*x*y);
	cudaMalloc((void**)&d_data, sizeof(unsigned char)*x*y*PIXEL_WIDTH);

	for(int i = 0; i < y; i++){
		for(int j = 0; j < x; j++){
			int idx = (i*x+j);
			int idx3 = idx*PIXEL_WIDTH;
			labels[idx] = 0;
			if((int)data[idx3] == 0) continue;
			labels[idx] = idx + 1;
		}
	}

	Table red;
	Table green;
	Table blue;

	initialize_table( red, HASH_ENTRIES, ELEMENTS );
	initialize_table( green, HASH_ENTRIES, ELEMENTS );
	initialize_table( blue, HASH_ENTRIES, ELEMENTS );

	int maxN = max(x,y);
	int phases = (int) ceil(log(maxN)/log(2)) + 1;

	cudaMemcpy(d_data, data, sizeof(unsigned char)*x*y*PIXEL_WIDTH, cudaMemcpyHostToDevice);
	cudaMemcpy(d_labels, labels, sizeof(int)*x*y, cudaMemcpyHostToDevice);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1); //Block size
	int linGrid_x = (int)ceil(x/(float)(BLOCK_SIZE-2)); //Threads per block -2 since we are loading ghost cells
	int linGrid_y = (int)ceil(y/(float)(BLOCK_SIZE-2));
	dim3 gridSize(linGrid_x,linGrid_y,1);

	printf("Applying segmentation...\n");
	cudaEventRecord(start);
	for(int pp = 0; pp <= phases; pp++){
		devicePhase<<<gridSize, blockSize>>>(d_labels, d_data, x, y);
		if(display){
			cudaMemcpy(labels, d_labels, sizeof(int)*x*y, cudaMemcpyDeviceToHost);
			doHashTableStuff(red, green, blue, x, y, labels, data, seed, num_threads);
			//Burada display
		}
	}
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	double t0, t1;
	if(!display){
		t0 = getTime();
		cudaMemcpy(labels, d_labels, sizeof(int)*x*y, cudaMemcpyDeviceToHost);
		doHashTableStuff(red, green, blue, x, y, labels, data, seed, num_threads);
		t1 = getTime();
	}
	
	cudaFree(d_labels);
	cudaFree(d_data);

	float kernelTime = 0;
	cudaEventElapsedTime(&kernelTime, start, stop);

	printf("Kernel time: %f\n", kernelTime / 1000);
	printf("Coloring time: %f\n", t1-t0);
	printf("Total segmentation time: %f\n", (t1-t0) + (kernelTime / 1000));
	printf("Writing segmented image %s...\n", output);

	printf("returned %d\n",stbi_write_png(output, x, y, n, data, 0));

	free_table( red );
	free_table( blue );
	free_table( green );

	stbi_image_free(data);
	return(0);
}

void cmdLine(int argc, char *argv[], char* i, char* output, int& d, int& num_threads, int &s){
/// Command line arguments
 // Default value of the domain sizes
 static struct option long_options[] = {
        {"i", required_argument, 0, 'i'},
        {"o", required_argument, 0, 'o'},
        {"display", required_argument, 0, 'd'},
        {"numthreads", required_argument, 0, 't'},
        {"seed", required_argument, 0, 's'},
 };
    // Process command line arguments
 int ac;
 for(ac=1;ac<argc;ac++) {
    int c;
    while ((c=getopt_long(argc,argv,"i:o:d:t:s:",long_options,NULL)) != -1){
        switch (c) {

	    // Name of input image
            case 'i':
            	strcpy(i, optarg);
                break;

	    // Nuber of threads
            case 't':
                num_threads = atoi(optarg);
                break;

	    // Turn on display
            case 'd':
                d = atoi(optarg);
                break;

	    // Output file name
            case 'o':
                output = strcpy(output, optarg);
                break;

        //Random seed
             case 's':
              	s = atoi(optarg);
               	break;

	    // Error
            default:
                printf("Usage: a.out [-i <input image name with format>] [-t <number of threads>]\n\t [-d turn on display]\n\t[-o <output file name with format>] [-s <random seed>]\n");
                exit(-1);
            }
    }
 }
}

void doHashTableStuff(Table red, Table green, Table blue, int x, int y, int *labels, unsigned char *data, int seed, int num_threads){
	srand(seed);
	int i, j;

	for (i = 0; i < y; i++) {
		for (j = 0; j < x; j++) {
			int label = labels[i*x+j];
			add_to_table( red, label, (int)rand()*255);
			add_to_table( green, label, (int)rand()*255);
			add_to_table( blue, label, (int)rand()*255);
		}
	}
	
	#pragma omp parallel for num_threads(num_threads) private(i,j) schedule(dynamic)
	for (i = 0; i < y; i++) {
		for (j = 0; j < x; j++) {
			int idx = i*x+j;
			int idx3 = idx*PIXEL_WIDTH;
			int label = labels[idx];
			data[idx3+0] = (char)getFromTable(red,label);
			data[idx3+1] = (char)getFromTable(blue,label);
			data[idx3+2] = (char)getFromTable(green,label);

		}
	}
}
