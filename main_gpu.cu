#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint-gcc.h>
#include <strings.h>

extern "C" {
#include "openFiles.h"
}

#include "gkernel.cuh"

#define MAX_ANTENNAS 20
#define MAX_ANTENNA_NAME 20
#define MAX_ANTENNA_FILE_NAME 50
#define NTHREADS 256

void parseConfig2(char *configFileName, int *nBit, int *nPol, int *isComplex, int *nChan, int *nant, double *lo,
                  double *bandwidth, int *numFFT, char **antenna, char **antFiles, double *delays,
                  double *antFileOffsets) {
    FILE *fconfig = fopen(configFileName, "r");

    if (fconfig == NULL) {
        printf("Error Opening File %s\n", configFileName);
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
                sscanf(token, "%lf", &((delays)[4 * iant + i])); // TODO CHECK
                token = strtok(NULL, " \n");
            }
            if (token == NULL) {
                fprintf(stderr, "Error: Insufficient data for antenna %d\n", iant);
                exit(1);
            }
            sscanf(token, "%lf", &((delays)[4 * iant + 3])); // Error checking needed
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

    // Check that the number of FFTs is a valid number
    if (*numFFT % 8) {
        printf("Error: numffts must be divisible by 8");
        exit(1);
    }
}

void initMem(char ***antenna, char ***antFiles, double ***delays, double **antFileOffsets) {
    // Initialize variables before passing to parseConfig
    *antenna = (char **) malloc(MAX_ANTENNAS * sizeof(char *));
    *antFiles = (char **) malloc(MAX_ANTENNAS * sizeof(char *));
    *antFileOffsets = (double *) malloc(MAX_ANTENNAS * sizeof(double));

    for (int i = 0; i < MAX_ANTENNAS; i++) {
        (*antenna)[i] = (char *) malloc(MAX_ANTENNA_NAME * sizeof(char));
        (*antFiles)[i] = (char *) malloc(MAX_ANTENNA_FILE_NAME * sizeof(char));
    }
}

void allocDataHost(uint8_t **data, int numantenna, int numchannels, int numffts, int nbit, int nPol, int iscomplex,
                   int *subintbytes) {
    int cfactor;

    if (iscomplex) {
        cfactor = 1;
    } else {
        cfactor = 2; // If real data FFT size twice size of number of frequecy channels
    }

    *subintbytes = numchannels * cfactor * (numffts + 1) * nbit / 8 * nPol;
    printf("Allocating %d MB per antenna per subint\n", *subintbytes / 1024 / 1024);
    printf("          %d MB total\n", (*subintbytes * numantenna) / 1024 / 1024);

    printf("subintbytes %d \n", *subintbytes);

    *data = (uint8_t *) malloc(numantenna * *subintbytes * sizeof(uint8_t));
    if (*data == NULL) {
        // Handle allocation failure
        fprintf(stderr, "Memory allocation failed\n");
    }

//    for (int a = 0; a < numantenna; a++) {
//        cudaError_t status = cudaHostAlloc((void **) &(*data)[a], *subintbytes * sizeof(uint8_t), cudaHostAllocDefault);
//        if (status != cudaSuccess) {
//            fprintf(stderr, "Unable to allocate %d bytes. Quitting\n", *subintbytes);
//            exit(1);
//        }
//    }
}

int readdata(int bytestoread, FILE **antStream, uint8_t *inputdata, int numantennas) {
    printf("Size antStream : %d \n", numantennas);
    for (int i = 0; i < numantennas; i++) {
        size_t bytes_read = fread(&inputdata[i * bytestoread], 1, bytestoread, antStream[i]);
        if (bytes_read != bytestoread) {
            if (feof(antStream[i])) {
                return 2; // End of file reached
            } else if (ferror(antStream[i])) {
                perror("Error: Problem reading data");
                return 1; // Error reading file
            }
        }
    }
    return 0; // Successful read
}

int main() {
    // Init Mem
    int numchannels, numantennas, nbaseline, numffts, nbit, iscomplex, nPol, samplegranularity, cfactor, subintbytes;
    double lo, bandwidth, sampletime, subinttime;
    char *configFileName = "./test4.conf";
    char **antennas = NULL, **antFiles = NULL;
    double *delays;
    double *antFileOffsets;
    uint8_t *inputdata;
    cudaEvent_t start_exec, stop_exec;
    cudaEventCreate(&start_exec);
    cudaEventCreate(&stop_exec);

    init_2bitLevels();

    initMem(&antennas, &antFiles, nullptr, &antFileOffsets);

    delays = (double *) malloc(MAX_ANTENNAS * 4 *  sizeof(double));

    // Parse Config
    printf("Parsing Config\n");

    parseConfig2(configFileName, &nbit, &nPol, &iscomplex, &numchannels, &numantennas, &lo, &bandwidth, &numffts,
                 antennas, antFiles, delays, antFileOffsets);

    printf("Del : %e %e %e %e \n", delays[0], delays[1], delays[2], delays[3]);
    printf("Del : %e %e %e %e \n", delays[4], delays[5], delays[6], delays[7]);
    printf("Del : %e %e %e %e \n", delays[8], delays[9], delays[10], delays[11]);
    printf("Del : %e %e %e %e \n", delays[12], delays[13], delays[14], delays[15]);


    samplegranularity = 8 / (nbit * nPol);
    if (samplegranularity < 1) {
        samplegranularity = 1;
    }
    nbaseline = numantennas * (numantennas - 1) / 2;
    if (iscomplex) {
        cfactor = 1;
    } else {
        cfactor = 2; // If real data FFT size twice size of number of frequecy channels
    }

    int fftsamples = numchannels * cfactor;
    int subintsamples = numffts * fftsamples;  // Number of time samples - need to factor # channels (pols) also
    printf("Subintsamples = %d\n", subintsamples);
    sampletime = 1.0 / bandwidth;

    if (!iscomplex) sampletime /= 2.0;
    subinttime = subintsamples * sampletime;
    printf("Subint = %f msec \n", subinttime * 1000);

    // Unpack
    int numkernelexecutions, executionsperthread;

    int unpackThreads;
    if (nbit == 2 && !iscomplex) {
        numkernelexecutions = fftsamples / 2;
    } else if (nbit == 8 && iscomplex) {
        numkernelexecutions = fftsamples * 2;
    } else {
        printf("Error: Unsupported number if bits/complex (%d/%d)", nbit, iscomplex);
        exit(1);
    }

    if (numkernelexecutions <= NTHREADS) {
        unpackThreads = numkernelexecutions;
        executionsperthread = 1;
    } else {
        unpackThreads = NTHREADS;
        executionsperthread = numkernelexecutions / (NTHREADS);
        if (numkernelexecutions % (NTHREADS)) {
            printf("Error: NTHREADS not divisible into number of kernel executions for unpack");
            exit(1);
        }
    }
    dim3 unpackBlocks = dim3(executionsperthread, numffts);


    // Set the number of blocks for when calculating the delay and phase info
    int delayPhaseThreads;
    numkernelexecutions = numffts;

    if (numkernelexecutions <= NTHREADS) {
        delayPhaseThreads = numkernelexecutions;
        executionsperthread = 1;
    } else {
        delayPhaseThreads = NTHREADS;
        executionsperthread = numkernelexecutions / NTHREADS;
        if (numkernelexecutions % NTHREADS) {
            printf("Error: NTHREADS not divisible into numkernelexecutions for delay and phase calculations");
            exit(1);
        }
    }
    dim3 delayPhaseBlocks = dim3(executionsperthread, numantennas);

    // Fractional Delay
    int fracDelayThreads;
    numkernelexecutions = numchannels;

    if (numkernelexecutions <= NTHREADS) {
        fracDelayThreads = numkernelexecutions;
        executionsperthread = 1;
    } else {
        fracDelayThreads = NTHREADS;
        executionsperthread = numkernelexecutions / NTHREADS;
        if (numkernelexecutions % NTHREADS) {
            printf("Error: NTHREADS not divisible into numkernelexecutions for fractional sample correction");
            exit(1);
        }
    }
    dim3 fracDelayBlocks = dim3(executionsperthread, numffts, numantennas);

    // Cross Correlation

    int targetThreads = 50e4;  // This seems a *lot*
    int parallelAccum = (int) ceil(targetThreads / numchannels + 1); // I suspect this has failure modes
    //cout << "Initial parallelAccum=" << parallelAccum << endl;
    while (parallelAccum && numffts % parallelAccum) parallelAccum--;
    if (parallelAccum == 0) {
        printf("Error: Could not determine block size for Cross Correlation \n");
        exit(1);
    }

    printf("Allocate Memory\n");
    // Allocate space in the buffers for the data and the delays
    printf("Allocating host data\n");
    allocDataHost(&inputdata, numantennas, numchannels, numffts, nbit, nPol, iscomplex, &subintbytes);

    // Allocate space on the GPU
    printf("Allocating GPU data\n");

    int8_t ***packedData;
    float **rotationPhaseInfo;
    float **fractionalSampleDelays;
    int32_t **sampleShifts;
    double **gpuDelays;
    COMPLEX **unpackedData, **channelisedData;
    cuComplex **baselineData;
    cufftHandle plan;
    FILE **antStream = (FILE **) malloc(numantennas * sizeof(FILE *));
    // Open Antenna Files
    if (antStream == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        return 1;
    }


    allocDataGPU(&packedData, &unpackedData, &channelisedData,
                 &baselineData, &rotationPhaseInfo, &fractionalSampleDelays, &sampleShifts,
                 &gpuDelays, numantennas, subintsamples,
                 nbit, nPol, iscomplex, numchannels, numffts, parallelAccum, 1, subintbytes);


    openFiles(&numantennas, antFiles, antStream);

    // Configure CUFFT
    if (cufftPlan1d(&plan, fftsamples, CUFFT_C2C, nPol * numantennas * numffts) != CUFFT_SUCCESS) {
        printf("CUFFT error: Plan creation failed");
        return (0);
    }

    // Read Data
    printf("Reading data\n");
    int status = readdata(subintbytes, antStream, inputdata, numantennas);
    printf("Status : %d \n", status);
    if (status) exit(1);
    init_2bitLevels();

    // Record the start time
    cudaEventRecord(start_exec, 0);

    int stream = 0;
    // Copy data to GPU
    printf("Copy data to GPU \n");
    for (int i = 0; i < numantennas; i++) {
        gpuErrchk(cudaMemcpyAsync(packedData[stream][i], &inputdata[i * subintbytes], subintbytes, cudaMemcpyHostToDevice));
    }

    // Copy delays to GPU
    printf("Copy delays to GPU \n");

    gpuErrchk(cudaMemcpyAsync(gpuDelays[stream], delays, 4 * numantennas * sizeof(double), cudaMemcpyHostToDevice));

    // Use the delays to calculate fringe rotation phases and fractional sample delays for each FFT //

    calculateDelaysAndPhases<<<delayPhaseBlocks, delayPhaseThreads>>>(gpuDelays[stream], lo,
                                                                                          sampletime, fftsamples,
                                                                                          numchannels,
                                                                                          samplegranularity,
                                                                                          rotationPhaseInfo[stream],
                                                                                          sampleShifts[stream],
                                                                                          fractionalSampleDelays[stream]);

    CudaCheckError();

    // Unpack the data

    printf("Unpack data \n");
    for (int i = 0; i < numantennas; i++) {
        if (nbit == 2 && !iscomplex) {
            unpack2bit_2chan_rotate<<<unpackBlocks, unpackThreads>>>(
                    &unpackedData[stream][2 * i * subintsamples], packedData[stream][i],
                    &rotationPhaseInfo[stream][i * numffts * 2], &(sampleShifts[stream][numffts * i]), fftsamples);
        } else if (nbit == 8 && iscomplex) {
            unpack8bitcomplex_2chan_rotate<<<unpackBlocks, unpackThreads>>>(
                    &unpackedData[stream][2 * i * subintsamples], packedData[stream][i],
                    &rotationPhaseInfo[stream][i * numffts * 2], &(sampleShifts[stream][numffts * i]), fftsamples);
        }
        CudaCheckError();
    }


    // FFT
    printf("Do FFT \n");
    //cufftSetStream(plan[stream], streams[stream]);

    if (cufftExecC2C(plan, unpackedData[stream], channelisedData[stream], CUFFT_FORWARD) != CUFFT_SUCCESS) {
        printf("CUFFT error: ExecC2C Forward failed");
        return(0);
    }


    // Fractional Delay Correction
    printf("Frac Delay Correct \n");

    //FracSampleCorrection<<<fracDelayBlocks,fracDelayThreads>>>(channelisedData[stream], fractionalSampleDelays[stream], numchannels, fftsamples, numffts, subintsamples);

    // TEST 1D

    // Set fracDelayThreads to numchannels
    fracDelayBlocks = dim3(numffts * numantennas, 1, 1);
    fracDelayThreads = numchannels;

    // Launch the kernel with updated grid and block configuration
    FracSampleCorrection1D<<<fracDelayBlocks, fracDelayThreads>>>(channelisedData[stream], fractionalSampleDelays[stream], numchannels, fftsamples, numffts, subintsamples);


    // END TEST

    CudaCheckError();

    // Cross correlate
    printf("Cross correlate \n");

    //copyAndIterate(channelisedData, 0, numantennas, nPol, subintsamples);

    gpuErrchk(cudaMemsetAsync(baselineData[stream], 0, nbaseline*4*numchannels*sizeof(cuComplex)));

    int ccblock_width = 128;
    //int nantxp = numantennas*2;
    int nantxp = numantennas;
    dim3 ccblock(1+(numchannels-1)/ccblock_width, nantxp-1, nantxp-1);
    CCAH3<<<ccblock, ccblock_width>>>(baselineData[stream], channelisedData[stream], numantennas, numffts, numchannels, fftsamples);
    //CrossCorrAccumHoriz<<<ccblock, ccblock_width, 0, streams[stream]>>>(baselineData[stream], channelisedData[stream], numantennas, numffts, numchannels, fftsamples);

    //copyAndIterateBaseline(baselineData, nbaseline, numchannels, parallelAccum, subintsamples);

    float dtime;
    cudaEventRecord(stop_exec, 0);
    cudaEventSynchronize(stop_exec);
    cudaEventElapsedTime(&dtime, start_exec, stop_exec);

    printf("Total Exec Time : %f \n", dtime);
    float rate = (float)subintsamples * numantennas * (2./cfactor) * nPol * nbit /(dtime/1000.)/1e9;
    printf("Processed %f sec of data (%f Gbps) \n", subinttime, rate);

    saveVisibilities("vis.out", baselineData[0], nbaseline, numchannels, numchannels, bandwidth);

    cudaDeviceSynchronize();
    cudaDeviceReset();
}
