//
// Created by emichel on 01/05/24.
//

#include "openFiles.h"

void openFiles(int numantennas, char **antFiles, FILE **antStream) {
    for (int i = 0; i < numantennas; i++) {
        antStream[i] = fopen(antFiles[i], "rb");
        if (antStream[i] == NULL) {
            fprintf(stderr, "Problem with file %s - does it exist?\n", antFiles[i]);
        }
    }
}

int readdata(int bytestoread, FILE **antStream, uint8_t **inputdata, int numStreams) {
    printf("Size of antStream : %d \n", numStreams);
    for (int i = 0; i < numStreams; i++) {
        size_t bytesRead = fread(inputdata[i], sizeof(uint8_t), bytestoread, antStream[i]);
        if (bytesRead != bytestoread) {
            if (feof(antStream[i])) {
                return 2;
            } else {
                fprintf(stderr, "Error: Problem reading data (Read %zu bytes)\n", bytesRead);
                return 1;
            }
        }
    }
    return 0;
}
