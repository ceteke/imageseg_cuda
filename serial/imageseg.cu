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
using namespace std;

#define PIXEL_WIDTH 3
#define TRESHOLD 4




#define SIZE    (100*1024*1024)
#define ELEMENTS    (SIZE / sizeof(unsigned int))
#define HASH_ENTRIES     2048

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

void cmdLine(int argc, char *argv[], char* i, char* output, int& d, int& num_threads, int& seed);
void doHashTableStuff(Table red, Table green, Table blue, int x, int y, int *labels, unsigned char *data, int seed);

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
		fprintf(stderr, "Couldn't load image.\n%s image not found.",i);
		return (-1);
	}

	printf("Running segmentation...\n");
	int *labels = (int *)malloc(sizeof(int)*x*y);
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
	double t0 = getTime();
	for(int pp = 0; pp <= phases; pp++){
		for(int i = 0; i < y; i++){
			for(int j = 0; j < x; j++){
				int idx = i*x + j;
				int idx3 = idx*PIXEL_WIDTH;
				if (labels[idx] == 0) continue;
				int ll = labels[idx];

				if(i != y-1 && abs((int)data[((i+1)*x + j)*PIXEL_WIDTH] - (int)data[idx3]) < TRESHOLD)
					labels[idx] = max(labels[idx], labels[(i+1)*x + j]);

				if(i != 0 && abs((int)data[((i-1)*x + j)*PIXEL_WIDTH] - (int)data[idx3]) < TRESHOLD)
					labels[idx] = max(labels[idx], labels[(i-1)*x + j]);

				if(i != y-1 && j != x-1 && abs((int)data[((i+1)*x + j + 1)*PIXEL_WIDTH] - (int)data[idx3]) < TRESHOLD)
					labels[idx] = max(labels[idx], labels[(i+1) * x + j + 1]);

				if(i != 0 && j != x-1 && abs((int)data[((i-1)*x + j + 1)*PIXEL_WIDTH] - (int)data[idx3]) < TRESHOLD)
					labels[idx] = max(labels[idx], labels[(i-1) * x + j + 1]);

				if(i != y-1 && j!= 0 && abs((int)data[((i+1)*x + j - 1)*PIXEL_WIDTH] - (int)data[idx3]) < TRESHOLD)
					labels[idx] = max(labels[idx], labels[(i+1) * x + j - 1]);

				if(i != 0 && j != 0 && abs((int)data[((i-1)*x + j - 1)*PIXEL_WIDTH] - (int)data[idx3]) < TRESHOLD)
					labels[idx] = max(labels[idx], labels[(i-1) * x + j - 1]);

				if (j != 0 && abs((int)data[(i*x + j - 1)*PIXEL_WIDTH] - (int)data[idx3]) < TRESHOLD)
					labels[idx] = max(labels[idx], labels[i*x + j - 1]);

				if (j != x-1 && abs((int)data[(i*x + j + 1)*PIXEL_WIDTH] - (int)data[idx3]) < TRESHOLD)
					labels[idx] = max(labels[idx], labels[i*x + j + 1]);

				if (ll < labels[idx]) {
					if (labels[ll - 1] < labels[idx])
						labels[ll - 1] = labels[idx];
				}
			}

		}

		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				int idx = i*x + j;
				if (labels[idx] != 0) {
					labels[idx] = max(labels[idx], labels[labels[idx] - 1]);
				}
			}
		}
		if(display){
			doHashTableStuff(red, green, blue, x, y, labels, data, seed);
			//Do the displaying here!
		}
	}
	
	double t1 = getTime();

	double t2, t3;
	if(!display){
		t2 = getTime();
		doHashTableStuff(red, green, blue, x, y, labels, data, seed);
		t3 = getTime();
	}
	
	printf("Labeling time: %f\n",(t1-t0));
	printf("Coloring time: %f\n",(t3-t2));
	printf("Total time: %f\n", (t1-t0)+(t3-t2));
	printf("Writing segmented image...\n");

	printf("returned %d\n",stbi_write_png(output, x, y, n, data, 0));

	free_table( red );
	free_table( blue );
	free_table( green );
	free(labels);
	stbi_image_free(data);
	return(0);
}

void doHashTableStuff(Table red, Table green, Table blue, int x, int y, int *labels, unsigned char *data, int seed){
	srand(seed);
	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			int label = labels[i*x+j];

			add_to_table( red, label, (int) random()*255 );
			add_to_table( green, label, (int) random()*255 );
			add_to_table( blue, label, (int) random()*255 );

		}
	}

	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			int label = labels[i*x+j];
			data[(i*x+j)*PIXEL_WIDTH+0] = (char)getFromTable(red,label);
			data[(i*x+j)*PIXEL_WIDTH+1] = (char)getFromTable(blue,label);
			data[(i*x+j)*PIXEL_WIDTH+2] = (char)getFromTable(green,label);

		}
	}
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
