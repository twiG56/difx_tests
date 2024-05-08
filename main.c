#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "allocDataHost.h"
#include "openFiles.h"
#include "kernel.h"

#define MAX_ANTENNAS 20
#define MAX_ANTENNA_NAME 20
#define MAX_ANTENNA_FILE_NAME 50

void parseConfig2(char *configFileName, int *nBit, int *nPol, int *isComplex, int *nChan, int *nant, double *lo,
                  double *bandwidth, int *numFFT, char **antenna, char **antFiles, double **delays,
                  double *antFileOffsets) {
    FILE *fconfig = fopen(configFileName, "r");

    if (fconfig == NULL) {
        printf("Error Opening File %s", configFileName);
        return;
    }

    char line[1024];
    int antToRead = 0;
    int iant = 0;

    // Set some defaults
    *nPol = 2;
    *isComplex = 0;
    *nBit = 2;

    // Read the config file
    while (fgets(line, 1024, fconfig)) {
        char keyword[1024];
        char *token = strtok(line, " \n");
        if (!token) {
            fprintf(stderr, "Error: Could not parse \"%s\"\n", line);
            exit(1);
        }
        strcpy(keyword, token);
        if (antToRead) {
            char thisfile[1024];
            sscanf(keyword, "%s", thisfile);
            strcpy(antenna[iant], keyword);
            token = strtok(NULL, " \n");
            strcpy(antFiles[iant], token);
            token = strtok(NULL, " \n");
            for (int i = 0; i < 3; i++) {
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for antenna %d\n", iant);
                    exit(1);
                }
                sscanf(token, "%lf", &((delays)[iant][i]));
                token = strtok(NULL, " \n");
            }
            if (token == NULL) {
                fprintf(stderr, "Error: Insufficient data for antenna %d\n", iant);
                exit(1);
            }
            sscanf(token, "%lf", &antFileOffsets[iant]); // Error checking needed
            iant++;
            antToRead--;
        } else {
            if (strcasecmp(keyword, "COMPLEX") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for COMPLEX\n");
                    exit(1);
                }
                sscanf(token, "%d", isComplex); // Should error check
            } else if (strcasecmp(keyword, "NBIT") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NBIT\n");
                    exit(1);
                }
                sscanf(token, "%d", nBit); // Should error check
            } else if (strcasecmp(keyword, "NPOL") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NPOL\n");
                    exit(1);
                }
                sscanf(token, "%d", nPol); // Should error check
            } else if (strcasecmp(keyword, "NCHAN") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NCHAN\n");
                    exit(1);
                }
                sscanf(token, "%d", nChan); // Should error check
            } else if (strcasecmp(keyword, "LO") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for LO\n");
                    exit(1);
                }
                sscanf(token, "%lf", lo); // Should error check
            } else if (strcasecmp(keyword, "BANDWIDTH") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for BANDWIDTH\n");
                    exit(1);
                }
                sscanf(token, "%lf", bandwidth); // Should error check
            } else if (strcasecmp(keyword, "NUMFFTS") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NUMFFTS\n");
                    exit(1);
                }
                sscanf(token, "%d", numFFT); // Should error check
            } else if (strcasecmp(keyword, "NANT") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NANT\n");
                    exit(1);
                }
                sscanf(token, "%d", nant); // Should error check
                antToRead = *nant;
                iant = 0;
            } else {
                fprintf(stderr, "Error: Unknown keyword \"%s\"\n", keyword);
            }
        }
    }
    fclose(fconfig);
}

void parseConfig(char *configFileName, int *nBit, int *nPol, int *isComplex, int *nChan, int *nant, double *lo,
                 double *bandwidth, int *numFFT, char **antenna, char **antFiles, double (*delays)[3],
                 double *antFileOffsets) {
    FILE *fconfig = fopen(configFileName, "r");

    if (fconfig == NULL) {
        printf("Error Opening File %s", configFileName);
        return;
    }

    char line[1024];
    int antToRead = 0;
    int iant = 0;

    // Set some defaults
    *nPol = 2;
    *isComplex = 0;
    *nBit = 2;

    // Read the config file
    while (fgets(line, 1024, fconfig)) {
        char keyword[1024];
        char *token = strtok(line, " \n");
        if (!token) {
            fprintf(stderr, "Error: Could not parse \"%s\"\n", line);
            exit(1);
        }
        strcpy(keyword, token);
        if (antToRead) {
            char thisfile[1024];
            sscanf(keyword, "%s", thisfile);
            strcpy(antenna[iant], keyword);
            token = strtok(NULL, " \n");
            strcpy(antFiles[iant], token);
            token = strtok(NULL, " \n");
            for (int i = 0; i < 3; i++) {
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for antenna %d\n", iant);
                    exit(1);
                }
                sscanf(token, "%lf", &delays[iant][i]); // Corrected allocation and access
                token = strtok(NULL, " \n");
            }
            if (token == NULL) {
                fprintf(stderr, "Error: Insufficient data for antenna %d\n", iant);
                exit(1);
            }
            sscanf(token, "%lf", &antFileOffsets[iant]); // Error checking needed
            iant++;
            antToRead--;
        } else {
            if (strcasecmp(keyword, "COMPLEX") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for COMPLEX\n");
                    exit(1);
                }
                sscanf(token, "%d", isComplex); // Should error check
            } else if (strcasecmp(keyword, "NBIT") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NBIT\n");
                    exit(1);
                }
                sscanf(token, "%d", nBit); // Should error check
            } else if (strcasecmp(keyword, "NPOL") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NPOL\n");
                    exit(1);
                }
                sscanf(token, "%d", nPol); // Should error check
            } else if (strcasecmp(keyword, "NCHAN") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NCHAN\n");
                    exit(1);
                }
                sscanf(token, "%d", nChan); // Should error check
            } else if (strcasecmp(keyword, "LO") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for LO\n");
                    exit(1);
                }
                sscanf(token, "%lf", lo); // Should error check
            } else if (strcasecmp(keyword, "BANDWIDTH") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for BANDWIDTH\n");
                    exit(1);
                }
                sscanf(token, "%lf", bandwidth); // Should error check
            } else if (strcasecmp(keyword, "NUMFFTS") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NUMFFTS\n");
                    exit(1);
                }
                sscanf(token, "%d", numFFT); // Should error check
            } else if (strcasecmp(keyword, "NANT") == 0) {
                token = strtok(NULL, " \n");
                if (token == NULL) {
                    fprintf(stderr, "Error: Insufficient data for NANT\n");
                    exit(1);
                }
                sscanf(token, "%d", nant); // Should error check
                antToRead = *nant;
                iant = 0;
            } else {
                fprintf(stderr, "Error: Unknown keyword \"%s\"\n", keyword);
            }
        }
    }
    fclose(fconfig);
}


void readConfig(int *nbit, int *nPol, int *iscomplex, int *nchan, int *nant, double *lo, double *bandwidth,
                int *numffts, char **antenna, char **antFiles, double **delays, double *antfileoffsets) {
    printf("nbit: %d\n", *nbit);
    printf("nPol: %d\n", *nPol);
    printf("iscomplex: %s\n", *iscomplex ? "true" : "false");
    printf("nchan: %d\n", *nchan);
    printf("nant: %d\n", *nant);
    printf("lo: %lf\n", *lo);
    printf("bandwidth: %lf\n", *bandwidth);
    printf("numffts: %d\n", *numffts);

    // Print antenna names
    printf("Antenna names:\n");
    for (int i = 0; i < *nant; i++) {
        printf("%s\n", antenna[i]);
    }

    // Print antenna files
    printf("Antenna files:\n");
    for (int i = 0; i < *nant; i++) {
        printf("%s\n", antFiles[i]);
    }

    // Print delays
    printf("Delays:\n");
    for (int i = 0; i < *nant; i++) {
        printf("Antenna %d:\n", i);
        for (int j = 0; j < 3; j++) {
            printf("%e ", delays[i][j]);
        }
        printf("\n");
    }

    // Print antfileoffsets
    printf("Antenna file offsets:\n");
    for (int i = 0; i < *nant; i++) {
        printf("%e\n", antfileoffsets[i]);
    }
}

void initMem(char ***antenna, char ***antFiles, double ***delays, double **antFileOffsets) {
    // Initialize variables before passing to parseConfig
    *antenna = (char **) malloc(MAX_ANTENNAS * sizeof(char *));
    *antFiles = (char **) malloc(MAX_ANTENNAS * sizeof(char *));
    *delays = (double **) malloc(MAX_ANTENNAS * sizeof(double *));
    *antFileOffsets = (double *) malloc(MAX_ANTENNAS * sizeof(double));

    for (int i = 0; i < MAX_ANTENNAS; i++) {
        (*antenna)[i] = (char *) malloc(MAX_ANTENNA_NAME * sizeof(char));
        (*antFiles)[i] = (char *) malloc(MAX_ANTENNA_FILE_NAME * sizeof(char));
        (*delays)[i] = (double *) malloc(3 * sizeof(double));
    }
}

int main() {
    char *configFileName = "./test16.conf"; // Provide the path to your config file
    int nBit, nPol, isComplex, nChan, nant, numFFT;
    double lo, bandwidth;
    char **antenna = NULL, **antFiles = NULL;
    double *antFileOffsets; // Changed to a regular array
    double **delays = NULL; // Changed to a 2D array

    initMem(&antenna, &antFiles, &delays, &antFileOffsets);

    // Call parseConfig function
    parseConfig2(configFileName, &nBit, &nPol, &isComplex, &nChan, &nant, &lo, &bandwidth, &numFFT, antenna, antFiles,
                 delays,
                 antFileOffsets);

    u_int8_t **inputdata = (u8 **) malloc(nant * sizeof(u8 *)); /**< the input data [numstations][subintbytes] */
    int subIntBytes = 0;

    // Output the parsed config

    readConfig(&nBit, &nPol, &isComplex, &nChan, &nant, &lo, &bandwidth, &numFFT, antenna, antFiles, delays,
               antFileOffsets);

    // Mem Alloc
    allocDataHost(&nChan, &numFFT, &nBit, &nPol, &isComplex, &subIntBytes);

    // File open

    FILE **antStream;

    antStream = (FILE **) malloc(nant * sizeof(FILE *));
    if (antStream == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        return 1;
    }

    openFiles(&nant, antFiles, antStream);

    readdata(&subIntBytes, antStream, inputdata, &nant, nant);

    // Checkpoint for timing
    struct timespec *starttime = (struct timespec *) malloc(sizeof(struct timespec));
    char *starttimestring = (char *) malloc(64 * sizeof(char));

    //FxKernel declarations
    int fftchannels, stridesize, substridesize, fractionalLoFreq, nbaselines;
    double sampletime;
    cf32 ***unpacked = (cf32 ***) malloc(nant * sizeof(cf32 **));
    cf32 ***channelised = (cf32 ***) malloc(nant * sizeof(cf32 **));
    cf32 ***conjchannels = (cf32 ***) malloc(nant * sizeof(cf32 **));
    cf32 ***visibilities = (cf32 ***) malloc((nant * (nant - 1) / 2) * sizeof(cf32 **));

    int *baselineCount = (int *) malloc((nant * (nant - 1) / 2) * sizeof(int));

    init_FxKernel(&nant, &nChan, &nBit, &lo, &bandwidth, starttime,
                  starttimestring, &fftchannels, &sampletime, &stridesize,
                  &substridesize, &fractionalLoFreq, unpacked, channelised, conjchannels, visibilities, &nbaselines, baselineCount);

    // Should be a SPLIT

    // END OF SPLIT

    long long diff_ms;
    int split_size = 2;

    cf32 ***vis_to_norm = (cf32 ***) malloc((nant * (nant - 1) / 2) * sizeof(cf32 **));
    int *baselineCount_to_norm = (int *) malloc((nant * (nant - 1) / 2) * sizeof(int));

        processAntennasAndBaseline(&nant, &numFFT, &fftchannels, antFileOffsets, &sampletime, inputdata, unpacked,
                                   channelised,
                                   conjchannels, delays, &nBit, &substridesize, &stridesize, &fractionalLoFreq, &lo,
                                   &nChan,
                                   visibilities, baselineCount, &bandwidth, vis_to_norm, baselineCount_to_norm);


    printf("nBaselines : %d \n", nbaselines);
    int nbaselines_split = nbaselines / split_size;
    cf32 ***vis_out = (cf32 ***) malloc((nant * (nant - 1) / 2) * sizeof(cf32 **));

    for (int i = 0; i < split_size; i++) {
        processNormalize(&nbaselines_split, &(baselineCount_to_norm[nbaselines_split * i]),
                         &(vis_to_norm[nbaselines_split * i]), &nChan, &(vis_out[nbaselines_split * i]));
    }

    endTiming(starttime, &diff_ms);

    saveVisibilities(&nbaselines, &nChan, vis_out, &bandwidth);

    saveLog(&diff_ms, starttimestring);

    return 0;
}
