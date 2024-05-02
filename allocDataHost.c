//
// Created by emichel on 01/05/24.
//

#include "allocDataHost.h"
#include "kernel.h"

void allocDataHost(uint8_t ***data, int numantenna, int numchannels, int numffts, int nbit, int nPol, int iscomplex, int *subintbytes)
{
    int cfactor;

    if (iscomplex)
    {
        cfactor = 1;
    }
    else
    {
        cfactor = 2; // If real data FFT size twice size of number of frequency channels
    }

    *subintbytes = numchannels * cfactor * (numffts + 1) * nbit / 8 * nPol;
    printf("subintbytes : %d \n", *subintbytes);
    printf("Allocating %d MB per antenna per subint\n", *subintbytes / (1024 * 1024));
    printf("           %d MB total\n", *subintbytes * numantenna / (1024 * 1024));

    *data = (uint8_t **)malloc(numantenna * sizeof(uint8_t *));
    if (*data == NULL)
    {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        exit(1);
    }

    for (int a = 0; a < numantenna; a++)
    {
        (*data)[a] = (uint8_t *)malloc(*subintbytes);
        if ((*data)[a] == NULL)
        {
            fprintf(stderr, "Unable to allocate %d bytes. Quitting\n", *subintbytes);
            exit(1);
        }
    }
}

void checkDataAlloc(uint8_t *** data, int *subintbytes, int *numantenna) {
    for(int i = 0; i < *numantenna; i++) {
        for(int j = 0; j < *subintbytes; j++) {
            (*data)[i][j] = i;
        }
    }
}

void freeMem(char ***antenna, char ***antFiles, int *nant, double **antFileOffsets, FILE ***antStream) {
    for (int i = 0; i < 20; i++) {
        //free(*antenna[i]);
        free(*antFiles[i]);
    }

    free(*antenna);
    free(*antFiles);
    free(*antFileOffsets);

    for (int i = 0; i < *nant; i++) {
        if (*antStream[i] != NULL) {
            fclose(*antStream[i]);
        }
    }
    free(*antStream);
}