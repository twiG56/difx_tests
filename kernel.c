//
// Created by emichel on 01/05/24.
//

#include "kernel.h"

// Initialization function
void init_FxKernel(int *nant, int *nchan, int *nbit, double *lo, double *bw, struct timespec *starttime,
              char *starttimestring, int *fftchannels, double *sampletime, int *stridesize,
              int *substridesize, int *fractionalLoFreq, cf32 ***unpacked, cf32 ***channelised, cf32 ***conjchannels,
              cf32 ***visibilities, int *nbaselines, int *baselineCount, int split) {

    *fftchannels = 2 * *nchan;
    *sampletime = 1.0 / (2.0 * *bw);
    int iscomplex = 0; // Allow for further generalization later
    int cfact = iscomplex ? 2 : 1;

    // Check for consistency and initialize lookup tables, if required
    if (*nbit == 2) {
        if (*fftchannels % 2) {
            fprintf(stderr, "Error: FFT length must be divisible by 2 for 2 bit data. Aborting\n");
            exit(1);
        }
        initLUT2bitReal();
    } else if (*nbit != 8) {
        fprintf(stderr, "Error: Do not support %d bits. Aborting\n", nbit);
        exit(1);
    }

    // Figure out the array stride size
    *stridesize = (int)sqrt(*nchan);
    if (*stridesize * *stridesize != *nchan) {
        fprintf(stderr, "Please choose a number of channels that is a square\n");
        exit(1);
    }
    *substridesize = (2 / cfact) * *stridesize;

    // Check if LO frequency has a fractional component
    *fractionalLoFreq = (*lo - (int)(*lo)) > TINY;

    // Allocate memory for unpacked array
    for (int i = 0; i < split * *nant; i++) {
        unpacked[i] = (cf32**)malloc(2 * sizeof(cf32*));
        if (unpacked[i] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }
        for (int j = 0; j < 2; j++) {
            unpacked[i][j] = vectorAlloc_cf32(*fftchannels);
            if (unpacked[i][j] == NULL) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                exit(1);
            }
        }
    }

    // Allocate memory for FFT'ed data and initialize FFT
    int order = 0;
    while((*fftchannels) >> order != 1)
    {
        order++;
    }

// Allocate memory for each element of channelised and conjchannels arrays
    for(int i = 0; i < split * *nant; i++)
    {
        channelised[i] = (cf32**)malloc(2 * sizeof(cf32*));
        if (channelised[i] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }

        conjchannels[i] = (cf32**)malloc(2 * sizeof(cf32*));
        if (conjchannels[i] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }

        for(int j = 0; j < 2; j++)
        {
            channelised[i][j] = (cf32*)vectorAlloc_cf32(*fftchannels);
            if (channelised[i][j] == NULL) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                exit(1);
            }

            conjchannels[i][j] = (cf32*)vectorAlloc_cf32(*nchan);
            if (conjchannels[i][j] == NULL) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                exit(1);
            }
        }
    }

    // Visibilities
    *nbaselines = *nant * (*nant - 1) / 2;
    for (int i = 0; i < split * *nbaselines; i++) {
        visibilities[i] = (cf32**)malloc(4 * sizeof(cf32*));
        if (visibilities[i] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }
        for (int j = 0; j < 4; j++) {
            visibilities[i][j] = vectorAlloc_cf32(*nchan);
            if (visibilities[i][j] == NULL) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                exit(1);
            }
        }
    }

    startTiming(starttime, starttimestring);

    resetBeforeProcess(nbaselines, visibilities, nchan, baselineCount, split);
}

void allocKernel(FxKernel *kernel, int *substridesize, int*stridesize, int* fftchannels, int* numchannels, double *sampletime, double *bw) {
    kernel->subtoff  = vectorAlloc_f64(*substridesize);
    kernel->subtval  = vectorAlloc_f64(*substridesize);
    kernel->subxoff  = vectorAlloc_f64(*substridesize);
    kernel->subxval  = vectorAlloc_f64(*substridesize);
    kernel->subphase = vectorAlloc_f64(*substridesize);
    kernel->subarg   = vectorAlloc_f32(*substridesize);
    kernel->subsin   = vectorAlloc_f32(*substridesize);
    kernel->subcos   = vectorAlloc_f32(*substridesize);
    kernel->steptoff  = vectorAlloc_f64(*stridesize);
    kernel->steptval  = vectorAlloc_f64(*stridesize);
    kernel->stepxoff  = vectorAlloc_f64(*stridesize);
    kernel->stepxval  = vectorAlloc_f64(*stridesize);
    kernel->stepphase = vectorAlloc_f64(*stridesize);
    kernel->steparg   = vectorAlloc_f32(*stridesize);
    kernel->stepsin   = vectorAlloc_f32(*stridesize);
    kernel->stepcos   = vectorAlloc_f32(*stridesize);
    kernel->stepcplx  = vectorAlloc_cf32(*stridesize);
    kernel->complexrotator = vectorAlloc_cf32(*fftchannels);

    // populate the fringe rotation arrays that can be pre-populated
    for(int i=0;i<*substridesize;i++)
    {
        kernel->subxoff[i] = ((double)i)/((double)*fftchannels);
        kernel->subtoff[i] = i**sampletime;
    }
    for(int i=0;i<*stridesize;i++)
    {
        kernel->stepxoff[i] = ((double)i* (double)*stridesize)/ ((double)*fftchannels);
        kernel->steptoff[i] = i* *stridesize * *sampletime;
    }

    kernel->subfracsamparg = vectorAlloc_f32(*stridesize);
    kernel->subfracsampsin = vectorAlloc_f32(*stridesize);
    kernel->subfracsampcos = vectorAlloc_f32(*stridesize);
    kernel->subchannelfreqs = vectorAlloc_f32(*stridesize);

    kernel->stepfracsamparg = vectorAlloc_f32(*stridesize);
    kernel->stepfracsampsin = vectorAlloc_f32(*stridesize);
    kernel->stepfracsampcos = vectorAlloc_f32(*stridesize);
    kernel->stepfracsampcplx = vectorAlloc_cf32(*stridesize);
    kernel->stepchannelfreqs = vectorAlloc_f32(*stridesize);
    kernel->fracsamprotator = vectorAlloc_cf32(*numchannels);

    // populate the channel frequency arrays
    for(int i=0;i<*stridesize;i++)
    {
        kernel->subchannelfreqs[i] = (float)((TWO_PI*i* *bw)/(double)*numchannels);
        kernel->stepchannelfreqs[i] = (float)((TWO_PI*i* *stridesize* *bw)/(double)*numchannels);
    }

    // FFTs
    // Get the size of, and initialise, the FFT
    int order = 0;
    while((*fftchannels) >> order != 1)
    {
        order++;
    }


    int sizeFFTSpec, sizeFFTInitBuf, wbufsize;
    u8 *fftInitBuf, *fftSpecBuf;
    ippsFFTGetSize_C_32fc(order, vecFFT_NoReNorm, vecAlgHintFast, &sizeFFTSpec, &sizeFFTInitBuf, &wbufsize);

    fftSpecBuf = ippsMalloc_8u(sizeFFTSpec);
    fftInitBuf = ippsMalloc_8u(sizeFFTInitBuf);
    kernel->fftbuffer = ippsMalloc_8u(wbufsize);

    ippsFFTInit_C_32fc(&kernel->pFFTSpecC, order, vecFFT_NoReNorm, vecAlgHintFast, fftSpecBuf, fftInitBuf);
    if (fftInitBuf) ippFree(fftInitBuf);
}

void initLUT2bitReal ()
{
    static const float HiMag = (float)3.3359;  // Optimal value
    const float lut4level[4] = {-HiMag, (float)-1.0, (float)1.0, HiMag};

    int l;
    for (int b = 0; b < 256; b++)	{
        for (int i = 0; i < 4; i++) {
            l = (b >> (2*i)) & 0x03;
            lut2bit[b][i].re = lut4level[l];
            lut2bit[b][i].im = 0;
        }
    }
}


void processAntennas(FxKernel *kernel, int *nant, int *numffts, int *fftchannels, double *antFileOffsets,
                     double *sampletime, int *antValid, u8 ** inputData, cf32 *** unpacked, cf32 ***channelised,
                     cf32 *** conjchannels, double **delays, int *nbit, int *substridesize, int *stridesize,
                     int *fractionalLoFreq, double *lo, int *nchan) {
    double meandelay, delaya, delayb, fractionaldelay, netdelaysamples_f;
    int netdelaysamples, offset;
    for(int i=0;i<*numffts;i++) {
        for (int j = 0; j < *nant; j++) {
            int maxoffset = (*numffts - 1) * *fftchannels;
            // unpack
            getStationDelay(j, i, &meandelay, &delaya, &delayb, delays);
            netdelaysamples_f = (meandelay - antFileOffsets[j]) / *sampletime;
            netdelaysamples = (int) (netdelaysamples_f + 0.5);

            fractionaldelay = -(netdelaysamples_f - netdelaysamples) * *sampletime;  // seconds
            offset = i * *fftchannels - netdelaysamples;
            if (offset == -1) // can happen due to changing geometric delay over the subint
            {
                ++(offset);
                fractionaldelay += *sampletime;
            }
            if (offset == maxoffset + 1) // can happen due to changing geometric delay over the subint
            {
                --(offset);
                fractionaldelay -= *sampletime;
            }
            if (offset < 0 || offset > maxoffset) {
                antValid[j] = 0;
                return;
            }
            antValid[j] = 1;
            unpack(inputData[j], unpacked[j], offset, *nbit, *fftchannels);

            // fringe rotate - after this function, each unpacked array has been fringe-rotated in-place
            fringerotate(kernel, unpacked[j], delaya, delayb, *substridesize, *stridesize, *lo, *fftchannels,
                         *fractionalLoFreq);

            // Channelise
            dofft(kernel, unpacked[j], channelised[j]);

            // If original data was real voltages, required channels will fill n/2+1 of the n channels, so move in-place

            // Fractional sample correct
            fracSampleCorrect(kernel, channelised[j], fractionaldelay, *stridesize, *nchan);

            // Calculate complex conjugate once, for efficency
            conjChannels(channelised[j], conjchannels[j], *nchan);
        }
    }
}

void processBaseline(int *nant, int *antValid, cf32 ***channelised, cf32 *** conjchannels, cf32 *** visibilities, int *nchan,
                     int *baselineCount) {
        int b = 0; // Baseline counter
        for (int j = 0; j < *nant - 1; j++) {
            if (!antValid[j]) {
                b += (*nant - (j + 1));
                continue;
            }
            for (int k = j + 1; k < *nant; k++) {
                if (!antValid[k]) {
                    (b)++;
                    continue;
                }

                for (int l = 0; l < 2; l++) {
                    for (int m = 0; m < 2; m++) {
                        // cross multiply + accumulate
                        vectorAddProduct_cf32(channelised[j][l], conjchannels[k][m],
                                              visibilities[b][2 * l + m], *nchan);
                    }
                }
                baselineCount[b]++;
                (b)++;
            }
        }
}

void processAntennasAndBaseline(int *nant, int *numffts, int *fftchannels, double *antFileOffsets,
                     double *sampletime, u8 ** inputData, cf32 *** unpacked, cf32 ***channelised,
                     cf32 *** conjchannels, double **delays, int *nbit, int *substridesize, int *stridesize,
                     int *fractionalLoFreq, double *lo, int *nchan, cf32 *** visibilities, int *baselineCount, double *bw,
                     cf32 *** visibilities_out, int *baselineCount_out) {
    FxKernel kernel;
    allocKernel(&kernel, substridesize, stridesize, fftchannels, nchan, sampletime, bw);
    
    double meandelay, delaya, delayb, fractionaldelay, netdelaysamples_f;
    int *antValid = (int *) malloc(*nant * sizeof(int));

    int netdelaysamples, offset;
    for(int i=0;i<*numffts;i++) {
        for (int j = 0; j < *nant; j++) {
            int maxoffset = (*numffts - 1) * *fftchannels;
            // unpack
            getStationDelay(j, i, &meandelay, &delaya, &delayb, delays);
            netdelaysamples_f = (meandelay - antFileOffsets[j]) / *sampletime;
            netdelaysamples = (int) (netdelaysamples_f + 0.5);

            fractionaldelay = -(netdelaysamples_f - netdelaysamples) * *sampletime;  // seconds
            offset = i * *fftchannels - netdelaysamples;
            if (offset == -1) // can happen due to changing geometric delay over the subint
            {
                ++(offset);
                fractionaldelay += *sampletime;
            }
            if (offset == maxoffset + 1) // can happen due to changing geometric delay over the subint
            {
                --(offset);
                fractionaldelay -= *sampletime;
            }
            if (offset < 0 || offset > maxoffset) {
                antValid[j] = 0;
                return;
            }
            antValid[j] = 1;
            unpack(inputData[j], unpacked[j], offset, *nbit, *fftchannels);

            printf("TEST \n", unpacked[j][0][0].re);



            // fringe rotate - after this function, each unpacked array has been fringe-rotated in-place
            fringerotate(&kernel, unpacked[j], delaya, delayb, *substridesize, *stridesize, *lo, *fftchannels,
                         *fractionalLoFreq);

            // Channelise
            dofft(&kernel, unpacked[j], channelised[j]);

            // If original data was real voltages, required channels will fill n/2+1 of the n channels, so move in-place

            // Fractional sample correct
            fracSampleCorrect(&kernel, channelised[j], fractionaldelay, *stridesize, *nchan);

            // Calculate complex conjugate once, for efficency
            conjChannels(channelised[j], conjchannels[j], *nchan);
        }
        processBaseline(nant, antValid, channelised, conjchannels, visibilities, nchan, baselineCount);
    }
    free(antValid);
    memcpy(visibilities_out, visibilities, (*nant * (*nant - 1) / 2) * sizeof(cf32***));
    memcpy(baselineCount_out, baselineCount, (*nant * (*nant - 1) / 2) * sizeof(int*));
}

void processNormalize(int *nbaselines, int *baselineCount,  cf32 *** visibilities, int *nchan, cf32 *** visibilities_out) {
    for (int i=0; i<*nbaselines; i++) {
        if (baselineCount[i] == 0) continue; // Really should flag data
        cf32 norm;
        norm.re = baselineCount[i];
        norm.im = 0;
        for (int j = 0; j < 4; j++) {
            vectorDivC_cf32_I(norm, visibilities[i][j], *nchan);
        }
    }
    memcpy(visibilities_out, visibilities, *nbaselines * sizeof(cf32***));
}

void resetBeforeProcess(int *nbaselines, cf32 *** visibilities, int *nchan, int *baselineCount, int split) {
    // Zero visibilities
    for (int i=0;i< split * *nbaselines; i++)
    {
        for (int j=0; j<4; j++)
        {
            vectorZero_cf32(visibilities[i][j], *nchan);
        }
    }
    memset(baselineCount, 0, split * sizeof(int)* *nbaselines); // Reset baselinecount
}

void getStationDelay(int antenna, int fftindex, double *meandelay, double *a, double *b, double **delays) {
    double *interpolator = delays[antenna];
    double d0, d1, d2;

    // Calculate values at the beginning, middle, and end of this FFT
    d0 = interpolator[0] * fftindex * fftindex + interpolator[1] * fftindex + interpolator[2];
    d1 = interpolator[0] * (fftindex + 0.5) * (fftindex + 0.5) + interpolator[1] * (fftindex + 0.5) + interpolator[2];
    d2 = interpolator[0] * (fftindex + 1.0) * (fftindex + 1.0) + interpolator[1] * (fftindex + 1.0) + interpolator[2];

    // Use these to calculate a linear interpolator across the FFT, as well as a mean value
    *a = d2 - d0;
    *b = d0 + (d1 - (*a * 0.5 + d0)) / 3.0;
    *meandelay = *a * 0.5 + *b;
}

void unpack(u8 * inputdata, cf32 ** unpacked, int offset, int nbit, int fftchannels)
{
    if (nbit==2) {
        unpackReal2bit(inputdata, unpacked, offset, fftchannels);
    } else {
        fprintf(stderr, "Unsupported number of bits!!!");
    }
}

void unpackReal2bit(u8 *inputdata, cf32 **unpacked, int offset, int nsamp) {
    cf32 *fp;
    int o = 0;
    u8 *byte = &inputdata[offset / 2];

    if (offset % 2) {
        fp = lut2bit[*byte];
        unpacked[0][o] = fp[2];
        unpacked[1][o] = fp[3];
        o++;
        byte++;
        nsamp--;
    }

    for (; o < nsamp - 1; o++) { // 2 time samples/byte
        fp = lut2bit[*byte];  // pointer to vector of 4 complex floats
        byte++;               // move pointer to next byte for next iteration of loop
        unpacked[0][o] = fp[0];
        unpacked[1][o] = fp[1];
        o++;
        unpacked[0][o] = fp[2];
        unpacked[1][o] = fp[3];
    }

    if (nsamp % 2) {
        fp = lut2bit[*byte];
        unpacked[0][o] = fp[0];
        unpacked[1][o] = fp[1];
    }
}

void fringerotate(FxKernel *kernel, cf32 ** unpacked, f64 a, f64 b, int substridesize, int stridesize, double lo, int fftchannels, double fractionalLoFreq)
{
    int integerdelay;
    int status;

    // subtract off any integer delay (whole seconds) present (of course, this shouldn't be present).
    integerdelay = (int)b;
    b -= integerdelay;

    // Fill in the delay values, using a and b and the precomputed offsets
    status = vectorMulC_f64(kernel->subxoff, a, kernel->subxval, substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate, subval multiplication\n");
    status = vectorMulC_f64(kernel->stepxoff, a, kernel->stepxval, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate, stepval multiplication\n");
    status = vectorAddC_f64_I(b, kernel->subxval, substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate, subval addition!!!\n");

    // Turn delay into turns of phase by multiplying by the lo
    status = vectorMulC_f64(kernel->subxval, lo, kernel->subphase, substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate lofreq sub multiplication!!!\n");
    status = vectorMulC_f64(kernel->stepxval, lo, kernel->stepphase, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate lofreq step multiplication!!!\n");
    if(fractionalLoFreq)
    {
        status = vectorAddC_f64_I((lo-(int)(lo))*(double)(integerdelay), kernel->subphase, substridesize);
        if(status != vecNoErr)
            fprintf(stderr, "Error in linearinterpolate lofreq non-integer freq addition!!!\n");
    }

    // Convert turns of phase into radians and bound into [0,2pi), then take sin/cos and assemble rotator vector
    for(int i=0;i<substridesize;i++)
    {
        kernel->subarg[i] = -TWO_PI*(kernel->subphase[i] - (int)(kernel->subphase[i]));
    }
    for(int i=0;i<stridesize;i++)
    {
        kernel->steparg[i] = -TWO_PI*(kernel->stepphase[i] - (int)(kernel->stepphase[i]));
    }
    status = vectorSinCos_f32(kernel->subarg, kernel->subsin, kernel->subcos, substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in sin/cos of sub rotate argument!!!\n");
    status = vectorSinCos_f32(kernel->steparg, kernel->stepsin, kernel->stepcos, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in sin/cos of step rotate argument!!!\n");
    status = vectorRealToComplex_f32(kernel->subcos, kernel->subsin, kernel->complexrotator, substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error assembling sub into complex!!!\n");
    status = vectorRealToComplex_f32(kernel->stepcos, kernel->stepsin, kernel->stepcplx, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error assembling step into complex!!!\n");
    for(int i=1;i<stridesize;i++)
    {
        status = vectorMulC_cf32(kernel->complexrotator, kernel->stepcplx[i], &kernel->complexrotator[i*substridesize], substridesize);
        if(status != vecNoErr)
            fprintf(stderr, "Error doing the time-saving complex multiplication!!!\n");
    }

    // Actually apply the fringe rotation to each polarisation in turn
    for (int i=0;i<2;i++)
    {
        status = vectorMul_cf32_I(kernel->complexrotator, unpacked[i], fftchannels);
        if (status != vecNoErr)
            fprintf(stderr, "Error in complex fringe rotation" );
    }
}

void dofft(FxKernel *kernel, cf32 ** unpacked, cf32 ** channelised) {
    // Do a single FFT on the 2 pols for a single antenna
    vecStatus status;

    for (int i=0; i<2; i++) {
        status = vectorFFT_CtoC_cf32(unpacked[i], channelised[i], kernel->pFFTSpecC, kernel->fftbuffer);
        if(status != vecNoErr) {
            fprintf(stderr, "Error calling FFT");
            exit(1);
        }
    }
}

void fracSampleCorrect(FxKernel *kernel, cf32 ** channelised, f64 fracdelay, int stridesize, int nchan)
{
    int status;

    // Create an array of phases for the fractional sample correction
    status = vectorMulC_f32(kernel->subchannelfreqs, fracdelay, kernel->subfracsamparg, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, arg generation (sub)!!! %d \n", status);

    status = vectorMulC_f32(kernel->stepchannelfreqs, fracdelay, kernel->stepfracsamparg, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, arg generation (step)!!! %d \n", status);

    //create the fractional sample correction array
    status = vectorSinCos_f32(kernel->subfracsamparg, kernel->subfracsampsin, kernel->subfracsampcos, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, sin/cos (sub)!!! %d \n", status);

    status = vectorSinCos_f32(kernel->stepfracsamparg, kernel->stepfracsampsin, kernel->stepfracsampcos, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, sin/cos (step)!!! %d \n", status);

    status = vectorRealToComplex_f32(kernel->subfracsampcos, kernel->subfracsampsin, kernel->fracsamprotator, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, real to complex (sub)!!! %d \n", status);

    status = vectorRealToComplex_f32(kernel->stepfracsampcos, kernel->stepfracsampsin, kernel->stepfracsampcplx, stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, real to complex (step)!!! %d \n", status);
    for(int i=1;i<stridesize;i++)
    {
        status = vectorMulC_cf32(kernel->fracsamprotator, kernel->stepfracsampcplx[i], &(kernel->fracsamprotator[i*stridesize]), stridesize);
        if(status != vecNoErr)
            fprintf(stderr, "Error doing the time-saving complex multiplication in frac sample correction!!! \n");
    }

    // Apply the fractional sample correction to each polarisation in turn
    for(int i=0;i<2;i++)
    {
        status = vectorMul_cf32_I(kernel->fracsamprotator, channelised[i], nchan);
        if(status != vecNoErr)
            fprintf(stderr, "Error in frac sample correction!!! %d \n", status);
    }
}

void conjChannels(cf32 ** channelised, cf32 ** conjchannels, int nchan) {
    // To avoid calculating this multiple times, generate the complex conjugate of the channelised data
    vecStatus status;

    for (int i=0; i<2; i++) {
        status = vectorConj_cf32(channelised[i], conjchannels[i], nchan); // Assumes USB and throws away 1 channel for real data
        if(status != vecNoErr) {
            fprintf(stderr, "Error calling vectorConj \n");
            exit(1);
        }
    }
}

void saveVisibilities(int *nbaselines, int *nchan, cf32 *** visibilities, double *bw) {
    f32 ***amp, ***phase;

    FILE *fvis = fopen("vis.out", "w");
    if (!fvis) {
        fprintf(stderr, "Error opening file for writing.\n");
        return;
    }

    amp = (f32***)malloc(*nbaselines * sizeof(f32**));
    phase = (f32***)malloc(*nbaselines * sizeof(f32**));
    if (!amp || !phase) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        return;
    }

    for (int i = 0; i < *nbaselines; i++) {
        amp[i] = (f32**)malloc(4 * sizeof(f32*)); // 2 pols and crosspol
        phase[i] = (f32**)malloc(4 * sizeof(f32*)); // 2 pols and crosspol
        if (!amp[i] || !phase[i]) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            return;
        }
        for (int j = 0; j < 4; j++) {
            amp[i][j] = vectorAlloc_f32(*nchan);
            phase[i][j] = vectorAlloc_f32(*nchan);
            if (!amp[i][j] || !phase[i][j]) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                return;
            }
            vectorMagnitude_cf32(visibilities[i][j], amp[i][j], *nchan);
            vectorPhase_cf32(visibilities[i][j], phase[i][j], *nchan);
        }
    }

    //printf("Test : %e \n", visibilities[0][0][0].re);
    //printf("Test : %e \n", visibilities[0][0][0].im);

    for (int c = 0; c < *nchan; c++) {
        fprintf(fvis, "%5d %11.6f", c, (c + 0.5) / *nchan * *bw / 1e6);
        for (int i = 0; i < *nbaselines; i++) {
            for (int j = 0; j < 4; j++) {
                fprintf(fvis, " %11.6f %11.6f", visibilities[i][j][c].re, visibilities[i][j][c].im);
                fprintf(fvis, " %11.6f %10.6f", amp[i][j][c], phase[i][j][c]);
            }
        }
        fprintf(fvis, "\n");
    }
    fclose(fvis);

    for (int i = 0; i < *nbaselines; i++) {
        for (int j = 0; j < 4; j++) {
            vectorFree(amp[i][j]);
            vectorFree(phase[i][j]);
        }
        free(amp[i]);
        free(phase[i]);
    }
    free(amp);
    free(phase);
}

void saveLog(long long *diff_ms, char *starttimestring) {
    // Open the file for writing
    FILE *ftiming = fopen("runtime.fxkernel.log", "w");
    if (ftiming == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        // Handle error as needed
    }

// Write to the file
    fprintf(ftiming, "Run time was %lld milliseconds\n", *diff_ms);
    fprintf(ftiming, "Start time was %s\n", starttimestring);

// Close the file
    fclose(ftiming);

}

void startTiming(struct timespec *starttime, char *starttimestring) {
    clock_gettime(CLOCK_MONOTONIC, starttime);
    time_t time_now_t = time(NULL);
    strftime(starttimestring, 64*sizeof(char), "%c", localtime(&time_now_t));
    starttimestring[strlen(starttimestring) - 1] = '\0'; // Removing the newline character at the end
}

void endTiming(struct timespec *starttime, long long *diff_ms) {
    struct timespec endtime;
    clock_gettime(CLOCK_MONOTONIC, &endtime);
    long long diff_ns = (endtime.tv_sec - (*starttime).tv_sec) * 1000000000LL + (endtime.tv_nsec - (*starttime).tv_nsec);
    *diff_ms = diff_ns / 1000000LL;

    printf("Run time was %lld milliseconds\n", *diff_ms);

    *diff_ms = diff_ns / 1000LL;

    printf("Run time was %lld microsec\n", *diff_ms);
}

void split_array(int *arr, int length, int **arr1, int **arr2) {
    int mid = length / 2;
    *arr1 = (int*)malloc(mid * sizeof(int));
    *arr2 = (int*)malloc((length - mid) * sizeof(int));

    // Copy elements to the first array
    memcpy(*arr1, arr, mid * sizeof(int));

    // Copy elements to the second array
    memcpy(*arr2, arr + mid, (length - mid) * sizeof(int));
}

void manageAntennasAndBaseline(int *nant, int *numffts, int *fftchannels, double *antFileOffsets,
                               double *sampletime, u8 ** inputData, cf32 *** unpacked, cf32 ***channelised,
                               cf32 *** conjchannels, double **delays, int *nbit, int *substridesize, int *stridesize,
                               int *fractionalLoFreq, double *lo, int *nchan, cf32 *** visibilities, int *baselineCount, double *bw,
                               cf32 *** visibilities_out, int *baselineCount_out) {

    // SPLIT TRY

    int *iter = (int *)malloc(*numffts * sizeof(int));
    int *iter_out = (int *)malloc(*numffts * sizeof(int));


    for(int i = 0; i < *numffts; i ++){
        iter[i] = i;
    }


    int split_size = 2;
    int numffts_split = *numffts / split_size;

    for (int i = 0; i < split_size; i++) {
        processAntennasAndBaselinePara(nant, &numffts_split, fftchannels, antFileOffsets, sampletime, inputData, unpacked,
                                       channelised, conjchannels, delays, nbit, substridesize, stridesize, fractionalLoFreq,
                                       lo, nchan, visibilities, baselineCount, &(iter[numffts_split*i]), bw, visibilities_out,
                                       baselineCount_out, 2);
    }
}

void processAntennasAndBaselinePara(int *nant, int *numffts, int *fftchannels, double *antFileOffsets,
                                double *sampletime, u8 ** inputData, cf32 *** unpacked, cf32 ***channelised,
                                cf32 *** conjchannels, double **delays, int *nbit, int *substridesize, int *stridesize,
                                int *fractionalLoFreq, double *lo, int *nchan, cf32 *** visibilities, int *baselineCount,
                                int *iteration, double *bw, cf32 *** visibilities_out, int *baselineCount_out, int split) {



    double meandelay, delaya, delayb, fractionaldelay, netdelaysamples_f;
    int *antValid = (int *) malloc(*nant * sizeof(int));

    clock_t start_station, start_unpack, start_fringe, start_dofft, start_frac, start_conj, start_base;
    long diff_station = 0L, diff_unpack = 0L, diff_fringe = 0L, diff_dofft = 0L, diff_frac = 0L, diff_conj = 0L, diff_base = 0L;

    int netdelaysamples, offset;

    for(int k = 0; k < *numffts; k++) {
        for (int j = 0; j < *nant; j++) {
            FxKernel kernel;
            allocKernel(&kernel, substridesize, stridesize, fftchannels, nchan, sampletime, bw);

            int maxoffset = (*numffts * split - 1) * *fftchannels;
            // unpack
            start_station = clock();
            getStationDelay(j, iteration[k], &meandelay, &delaya, &delayb, delays);
            diff_station += clock() - start_station;

            netdelaysamples_f = (meandelay - antFileOffsets[j]) / *sampletime;
            netdelaysamples = (int) (netdelaysamples_f + 0.5);

            fractionaldelay = -(netdelaysamples_f - netdelaysamples) * *sampletime;  // seconds
            offset = iteration[k] * *fftchannels - netdelaysamples;

            if (offset == -1) // can happen due to changing geometric delay over the subint
            {
                ++(offset);
                fractionaldelay += *sampletime;
            }
            if (offset == maxoffset + 1) // can happen due to changing geometric delay over the subint
            {
                --(offset);
                fractionaldelay -= *sampletime;
            }
            if (offset < 0 || offset > maxoffset) {
                antValid[j] = 0;
                continue;
            }
            antValid[j] = 1;
            start_unpack = clock();
            unpack(inputData[j], unpacked[j], offset, *nbit, *fftchannels);

            diff_unpack += clock() - start_unpack;

            // fringe rotate - after this function, each unpacked array has been fringe-rotated in-place
            start_fringe = clock();
            fringerotate(&kernel, unpacked[j], delaya, delayb, *substridesize, *stridesize, *lo, *fftchannels,
                         *fractionalLoFreq);

            diff_fringe += clock() - start_fringe;

            // Channelise
            start_dofft = clock();
            dofft(&kernel, unpacked[j], channelised[j]);


            diff_dofft += clock() - start_dofft;

            // If original data was real voltages, required channels will fill n/2+1 of the n channels, so move in-place

            // Fractional sample correct
            start_frac = clock();
            fracSampleCorrect(&kernel, channelised[j], fractionaldelay, *stridesize, *nchan);

            diff_frac += clock() - start_frac;

            // Calculate complex conjugate once, for efficency
            start_conj = clock();
            conjChannels(channelised[j], conjchannels[j], *nchan);

            diff_conj += clock() - start_conj;
        }

//        printf("test : %e \n", channelised[0][0][0].re);
//        printf("test : %e \n", channelised[0][0][1].re);
//
//        printf("test : %e \n", channelised[1][0][0].re);
//        printf("test : %e \n", channelised[1][0][1].re);
//
//        printf("test : %e \n", channelised[1][1][0].re);
//        printf("test : %e \n", channelised[1][1][1].re);
//        exit(0);

        start_base = clock();
        processBaselinePara(nant, antValid, channelised, conjchannels, visibilities, nchan, baselineCount);

        diff_base += clock() - start_base;
    }

    printf("test : %e \n", channelised[0][0][0].re);
    printf("test : %e \n", channelised[0][0][1].re);
    printf("test : %e \n", channelised[1][0][0].re);
    printf("test : %e \n", channelised[1][0][1].re);



    printf("Time taken station in ms : %f \n", ((double) (1000 * diff_station)) / (CLOCKS_PER_SEC) );
    printf("Time taken unpack in ms : %f \n", ((double) (1000 * diff_unpack)) / ((CLOCKS_PER_SEC) * 1));
    printf("Time taken fringe in ms : %f \n", ((double) (1000 * diff_fringe)) / (CLOCKS_PER_SEC) );
    printf("Time taken fft in ms : %f \n", ((double) (1000 * diff_dofft)) / (CLOCKS_PER_SEC) );
    printf("Time taken frac in ms : %f \n", ((double) (1000 * diff_frac)) / (CLOCKS_PER_SEC) );
    printf("Time taken conj in ms : %f \n", ((double) (1000 * diff_conj)) / (CLOCKS_PER_SEC) );

    printf("Time taken base in ms : %f \n", ((double) (1000 * diff_base)) / (CLOCKS_PER_SEC) );


    free(antValid);
    memcpy(visibilities_out, visibilities, (*nant * (*nant - 1) / 2) * sizeof(cf32***));
    memcpy(baselineCount_out, baselineCount, (*nant * (*nant - 1) / 2) * sizeof(int));
}

void processBaselinePara(int *nant, int *antValid, cf32 ***channelised, cf32 *** conjchannels, cf32 *** visibilities, int *nchan,
                     int *baselineCount) {
    int b = 0; // Baseline counter
    for (int j = 0; j < *nant - 1; j++) {
        if (!antValid[j]) {
            b += (*nant - (j + 1));
            continue;
        }
        for (int k = j + 1; k < *nant; k++) {
            if (!antValid[k]) {
                (b)++;
                continue;
            }

            for (int l = 0; l < 2; l++) {
                for (int m = 0; m < 2; m++) {
                    // cross multiply + accumulate
                    vectorAddProduct_cf32(channelised[j][l], conjchannels[k][m],
                                          visibilities[b][2 * l + m], *nchan);
                }
            }
            baselineCount[b]++;
            (b)++;
        }
    }
}

void merge(int split, int nbaselines, cf32 ***visibilities, int *nant, int *nChan, int *baselineCount, cf32 ***visibilities_out, int *baselineCount_out) {
    printf("Test : %f \n", visibilities[0][0][0].re);
    //printf("Test : %f \n", visibilities[120][0][0].re);
    //printf("Test : %f \n", visibilities[240][0][0].re);
    //printf("Test : %f \n", visibilities[360][0][0].re);
    if(split != 1) {
        for(int i = 0; i < nbaselines; i++) {
            for(int j = 0; j < 4; j++) {
                for(int k = 1; k < split; k++) {
                    vectorAdd_cf32_I(visibilities[i + k * (*nant * (*nant - 1) / 2)][j], visibilities[i][j], *nChan);
                }
            }
            for(int k = 1; k < split; k++) {
                baselineCount[i] += baselineCount[i + k * nbaselines];
            }
        }
    }
    printf("Test B : %d \n", baselineCount[0]);

    printf("Test : %f \n", visibilities[0][0][0].re);
    memcpy(visibilities_out, visibilities, (*nant * (*nant - 1) / 2) * sizeof(cf32***));
    memcpy(baselineCount_out, baselineCount, (*nant * (*nant - 1) / 2) * sizeof(int));
}


void unpackImpl(cf32** unpacked, int split, int *nant, int *fftchannels, u8* inputData, int *offset, int *nbit){
    unpacked = (cf32**)malloc(2 * sizeof(cf32*));
    if (unpacked == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        exit(1);
    }
    for (int j = 0; j < 2; j++) {
        unpacked[j] = vectorAlloc_cf32(*fftchannels);
        if (unpacked[j] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }
    }
    unpack(inputData, unpacked, *offset, *nbit, *fftchannels);
}

void fringeRotateImpl(FxKernel *kernel, cf32 ** unpacked, f64 a, f64 b, int substridesize, int stridesize, double lo,
                      int fftchannels, double fractionalLoFreq, cf32** unpacked_out, FxKernel *kernel_out) {

    fringerotate(kernel, unpacked, a, b, substridesize, stridesize, lo, fftchannels, fractionalLoFreq);

    memcpy(unpacked_out, unpacked, 2 * sizeof(cf32*));
    memcpy(kernel_out, kernel, sizeof(FxKernel));
}

void doFFTImpl(FxKernel * kernel, cf32** unpacked, cf32** channelised, int *fftchannels) {
    for(int j = 0; j < 2; j++) {
        channelised[j] = (cf32 *) vectorAlloc_cf32(*fftchannels);
        if (channelised[j] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }
    }
    dofft(kernel, unpacked, channelised);
}

void fracSampleCorrectImpl(FxKernel * kernel, cf32** channelised, double *fractionaldelay, int *stridesize, int *nchan,
                           cf32** channelised_out) {
    fracSampleCorrect(kernel, channelised, *fractionaldelay, *stridesize, *nchan);

    memcpy(channelised_out, channelised, 2 * sizeof(cf32*));
}

void conjChannelsImpl(cf32** channelised, cf32** conjchannels, int *nchan, int *fftchannels, cf32** channelised_out) {
    for(int j = 0; j < 2; j++) {
        channelised[j] = (cf32 *) vectorAlloc_cf32(*fftchannels);
        if (channelised[j] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }
    }

    conjChannels(channelised, conjchannels, *nchan);

    memcpy(channelised_out, channelised, 2 * sizeof(cf32*));
}

void stationDelayAndOffset(double ** delays, int *iteration, double *antFileOffsets, int *antValid, int *offset,
                           double *fracDelay, double *delaya, double* delayb, double *sampletime, int *fftchannels,
                           int numffts) {
    double meandelay, netdelaysamples_f;
    int netdelaysamples;
    int maxoffset = (numffts - 1) * *fftchannels;
    int j = 0; // TODO : change
    getStationDelay(j, *iteration, &meandelay, delaya, delayb, delays);

    netdelaysamples_f = (meandelay - antFileOffsets[j]) / *sampletime;
    netdelaysamples = (int) (netdelaysamples_f + 0.5);

    *fracDelay = -(netdelaysamples_f - netdelaysamples) * *sampletime;  // seconds
    *offset = *iteration * *fftchannels - netdelaysamples;
    if (*offset == -1) // can happen due to changing geometric delay over the subint
    {
        ++(*offset);
        *fracDelay += *sampletime;
    }
    if (*offset == maxoffset + 1) // can happen due to changing geometric delay over the subint
    {
        --(*offset);
        *fracDelay -= *sampletime;
    }
    if (*offset < 0 || *offset > maxoffset) {
        antValid[j] = 0;
        return;
    }
    antValid[j] = 1;
}