//
// Created by emichel on 01/05/24.
//

#include "kernel.h"

// Initialization function
void init_FxKernel(FxKernel *kernel, int nant, int nchan, int nfft, int nbit, double lo, double bw, u8 ** inputData) {
    kernel->numantennas = nant;
    kernel->numchannels = nchan;
    kernel->fftchannels = 2 * nchan;
    kernel->numffts = nfft;
    kernel->nbits = nbit;
    kernel->lofreq = lo;
    kernel->bandwidth = bw;
    kernel->sampletime = 1.0 / (2.0 * bw);
    kernel->iscomplex = 0; // Allow for further generalization later
    kernel->cfact = kernel->iscomplex ? 2 : 1;

    // Check for consistency and initialize lookup tables, if required
    if (kernel->nbits == 2) {
        if (kernel->fftchannels % 2) {
            fprintf(stderr, "Error: FFT length must be divisible by 2 for 2 bit data. Aborting\n");
            exit(1);
        }
        initLUT2bitReal();
    } else if (kernel->nbits != 8) {
        fprintf(stderr, "Error: Do not support %d bits. Aborting\n", kernel->nbits);
        exit(1);
    }

    // Figure out the array stride size
    kernel->stridesize = (int)sqrt(kernel->numchannels);
    if (kernel->stridesize * kernel->stridesize != kernel->numchannels) {
        fprintf(stderr, "Please choose a number of channels that is a square\n");
        exit(1);
    }
    kernel->substridesize = (2 / kernel->cfact) * kernel->stridesize;

    // Check if LO frequency has a fractional component
    kernel->fractionalLoFreq = (lo - (int)lo) > TINY;

    // Allocate memory for unpacked array
    kernel->unpacked = (cf32***)malloc(nant * sizeof(cf32**));
    if (kernel->unpacked == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        exit(1);
    }
    for (int i = 0; i < nant; i++) {
        kernel->unpacked[i] = (cf32**)malloc(2 * sizeof(cf32*));
        if (kernel->unpacked[i] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }
        for (int j = 0; j < 2; j++) {
            kernel->unpacked[i][j] = vectorAlloc_cf32(kernel->fftchannels);
            if (kernel->unpacked[i][j] == NULL) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                exit(1);
            }
        }
    }

    //allocate the arrays for holding the fringe rotation vectors
    kernel->subtoff  = vectorAlloc_f64(kernel->substridesize);
    kernel->subtval  = vectorAlloc_f64(kernel->substridesize);
    kernel->subxoff  = vectorAlloc_f64(kernel->substridesize);
    kernel->subxval  = vectorAlloc_f64(kernel->substridesize);
    kernel->subphase = vectorAlloc_f64(kernel->substridesize);
    kernel->subarg   = vectorAlloc_f32(kernel->substridesize);
    kernel->subsin   = vectorAlloc_f32(kernel->substridesize);
    kernel->subcos   = vectorAlloc_f32(kernel->substridesize);
    kernel->steptoff  = vectorAlloc_f64(kernel->stridesize);
    kernel->steptval  = vectorAlloc_f64(kernel->stridesize);
    kernel->stepxoff  = vectorAlloc_f64(kernel->stridesize);
    kernel->stepxval  = vectorAlloc_f64(kernel->stridesize);
    kernel->stepphase = vectorAlloc_f64(kernel->stridesize);
    kernel->steparg   = vectorAlloc_f32(kernel->stridesize);
    kernel->stepsin   = vectorAlloc_f32(kernel->stridesize);
    kernel->stepcos   = vectorAlloc_f32(kernel->stridesize);
    kernel->stepcplx  = vectorAlloc_cf32(kernel->stridesize);
    kernel->complexrotator = vectorAlloc_cf32(kernel->fftchannels);

    // populate the fringe rotation arrays that can be pre-populated
    for(int i=0;i<kernel->substridesize;i++)
    {
        kernel->subxoff[i] = ((double)i)/((double)kernel->fftchannels);
        kernel->subtoff[i] = i*kernel->sampletime;
    }
    for(int i=0;i<kernel->stridesize;i++)
    {
        kernel->stepxoff[i] = ((double)i* (double)kernel->stridesize)/ ((double)kernel->fftchannels);
        kernel->steptoff[i] = i*kernel->stridesize*kernel->sampletime;
    }

    // Allocate memory for FFT'ed data and initialize FFT
    int order = 0;
    while((kernel->fftchannels) >> order != 1)
    {
        order++;
    }

// Allocate memory for channelised and conjchannels arrays
    kernel->channelised = (cf32***)malloc(nant * sizeof(cf32**));
    if (kernel->channelised == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        exit(1);
    }

    kernel->conjchannels = (cf32***)malloc(nant * sizeof(cf32**));
    if (kernel->conjchannels == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        exit(1);
    }

// Allocate memory for each element of channelised and conjchannels arrays
    for(int i = 0; i < nant; i++)
    {
        kernel->channelised[i] = (cf32**)malloc(2 * sizeof(cf32*));
        if (kernel->channelised[i] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }

        kernel->conjchannels[i] = (cf32**)malloc(2 * sizeof(cf32*));
        if (kernel->conjchannels[i] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }

        for(int j = 0; j < 2; j++)
        {
            kernel->channelised[i][j] = (cf32*)vectorAlloc_cf32(kernel->fftchannels);
            if (kernel->channelised[i][j] == NULL) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                exit(1);
            }

            kernel->conjchannels[i][j] = (cf32*)vectorAlloc_cf32(kernel->numchannels);
            if (kernel->conjchannels[i][j] == NULL) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                exit(1);
            }
        }
    }

    // Get the size of, and initialise, the FFT
    int sizeFFTSpec, sizeFFTInitBuf, wbufsize;
    u8 *fftInitBuf, *fftSpecBuf;
    ippsFFTGetSize_C_32fc(order, vecFFT_NoReNorm, vecAlgHintFast, &sizeFFTSpec, &sizeFFTInitBuf, &wbufsize);
    fftSpecBuf = ippsMalloc_8u(sizeFFTSpec);
    fftInitBuf = ippsMalloc_8u(sizeFFTInitBuf);
    kernel->fftbuffer = ippsMalloc_8u(wbufsize);
    ippsFFTInit_C_32fc(&kernel->pFFTSpecC, order, vecFFT_NoReNorm, vecAlgHintFast, fftSpecBuf, fftInitBuf);
    if (fftInitBuf) ippFree(fftInitBuf);

    // Visibilities
    kernel->nbaselines = nant * (nant - 1) / 2;
    kernel->visibilities = (cf32***)malloc(kernel->nbaselines * sizeof(cf32**));
    if (kernel->visibilities == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        exit(1);
    }
    for (int i = 0; i < kernel->nbaselines; i++) {
        kernel->visibilities[i] = (cf32**)malloc(4 * sizeof(cf32*));
        if (kernel->visibilities[i] == NULL) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            exit(1);
        }
        for (int j = 0; j < 4; j++) {
            kernel->visibilities[i][j] = vectorAlloc_cf32(kernel->numchannels);
            if (kernel->visibilities[i][j] == NULL) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                exit(1);
            }
        }
    }

    // also the channel frequency arrays and the other fractional sample correction arrays
    kernel->subchannelfreqs = vectorAlloc_f32(kernel->stridesize);
    kernel->stepchannelfreqs = vectorAlloc_f32(kernel->stridesize);
    kernel->subfracsamparg = vectorAlloc_f32(kernel->stridesize);
    kernel->subfracsampsin = vectorAlloc_f32(kernel->stridesize);
    kernel->subfracsampcos = vectorAlloc_f32(kernel->stridesize);
    kernel->stepfracsamparg = vectorAlloc_f32(kernel->stridesize);
    kernel->stepfracsampsin = vectorAlloc_f32(kernel->stridesize);
    kernel->stepfracsampcos = vectorAlloc_f32(kernel->stridesize);
    kernel->stepfracsampcplx = vectorAlloc_cf32(kernel->stridesize);
    kernel->fracsamprotator = vectorAlloc_cf32(kernel->numchannels);

    // populate the channel frequency arrays
    for(int i=0;i<kernel->stridesize;i++)
    {
        kernel->subchannelfreqs[i] = (float)((TWO_PI*i*kernel->bandwidth)/(double)kernel->numchannels);
        kernel->stepchannelfreqs[i] = (float)((TWO_PI*i*kernel->stridesize*kernel->bandwidth)/(double)kernel->numchannels);
    }

    // Antenna and baseline validity/weights
    kernel->antValid = (int*)malloc(nant * sizeof(int));
    if (kernel->antValid == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        exit(1);
    }

    kernel->baselineCount = (int*)malloc(kernel->nbaselines * sizeof(int));
    if (kernel->baselineCount == NULL) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        exit(1);
    }

    // Give the fxkernel its pointer to the input data
    kernel->inputdata = inputData;
    if (kernel->inputdata == NULL) {
        fprintf(stderr, "Input data is NULL. Quitting.\n");
        exit(1);
    }
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

void setDelays(FxKernel *kernel, double ** d, double * f)
{
    kernel->delays = d;
    kernel->filestartoffsets = f;
}


void processAntennas(FxKernel *kernel, int i, int j, double *meandelay, double *fractionaldelay, double *delaya,
                     double *delayb, double *netdelaysamples_f, int *netdelaysamples, int *offset, int maxoffset) {
    // unpack
    getStationDelay(kernel, j, i, meandelay, delaya, delayb);
    *netdelaysamples_f = (*meandelay - kernel->filestartoffsets[j]) / kernel->sampletime;
    *netdelaysamples   = (int)(*netdelaysamples_f + 0.5);

    *fractionaldelay = -(*netdelaysamples_f - *netdelaysamples)*kernel->sampletime;  // seconds
    *offset = i*kernel->fftchannels - *netdelaysamples;
    if(*offset == -1) // can happen due to changing geometric delay over the subint
    {
        ++(*offset);
        *fractionaldelay += kernel->sampletime;
    }
    if(*offset == maxoffset+1) // can happen due to changing geometric delay over the subint
    {
        --(*offset);
        *fractionaldelay -= kernel->sampletime;
    }
    if(*offset < 0 || *offset>maxoffset)
    {
        kernel->antValid[j] = 0;
        return;
    }
    kernel->antValid[j] = 1;
    unpack(kernel, kernel->inputdata[j], kernel->unpacked[j], *offset);

    // fringe rotate - after this function, each unpacked array has been fringe-rotated in-place
    fringerotate(kernel, kernel->unpacked[j], *delaya, *delayb);

    // Channelise
    dofft(kernel, kernel->unpacked[j], kernel->channelised[j]);

    // If original data was real voltages, required channels will fill n/2+1 of the n channels, so move in-place

    // Fractional sample correct
    fracSampleCorrect(kernel, kernel->channelised[j], *fractionaldelay);

    // Calculate complex conjugate once, for efficency
    conjChannels(kernel, kernel->channelised[j], kernel->conjchannels[j]);
}

void process(FxKernel *kernel)
{
    // delay variables
    double meandelay; //mean delay in the middle of the FFT for a given antenna, in seconds
    double fractionaldelay; // fractional component of delay to be correction after channelisation
    double delaya, delayb; // coefficients a and b for the delay interpolation across a given FFT
    double netdelaysamples_f; // net samples delays (floating point): mean delay minus file start offset converted into units of sample time
    int netdelaysamples; // integer number of net samples delay
    int offset; //offset into the packed data vector

    // Zero visibilities
    for (int i=0;i<kernel->nbaselines; i++)
    {
        for (int j=0; j<4; j++)
        {
            vectorZero_cf32(kernel->visibilities[i][j], kernel->numchannels);
        }
    }

    int maxoffset = (kernel->numffts-1)*kernel->fftchannels;
    memset(kernel->baselineCount, 0, sizeof(int)*kernel->nbaselines); // Reset baselinecount

    // for(number of FFTs)... (parallelised via pthreads?)
    for(int i=0;i<kernel->numffts;i++)
    {
        // do station-based processing for each antenna in turn
        for(int j=0;j<kernel->numantennas;j++)
        {
            processAntennas(kernel, i, j, &meandelay, &fractionaldelay, &delaya, &delayb, &netdelaysamples_f, &netdelaysamples, &offset, maxoffset);
        }

        // then do the baseline based processing   (CJP: What about autos)
        int b = 0; // Baseline counter
        for(int j=0;j<kernel->numantennas-1;j++)
        {
            if (!kernel->antValid[j])
            {
                b += (kernel->numantennas-(j+1));
                continue;
            }
            for(int k=j+1;k<kernel->numantennas;k++)
            {
                if (!kernel->antValid[k])
                {
                    b++;
                    continue;
                }

                for(int l=0;l<2;l++)
                {
                    for(int m=0;m<2;m++)
                    {
                        // cross multiply + accumulate
                        vectorAddProduct_cf32(kernel->channelised[j][l], kernel->conjchannels[k][m], kernel->visibilities[b][2*l + m], kernel->numchannels);
                    }
                }
                kernel->baselineCount[b]++;
                b++;
            }
        }
    }

    // Normalise
    cf32 norm;
    for (int i=0; i<kernel->nbaselines; i++) {
        if (kernel->baselineCount[i]==0) continue; // Really should flag data
        norm.re = kernel->baselineCount[i];
        norm.im = 0;
        for (int j=0; j<4; j++) {
            vectorDivC_cf32_I(norm, kernel->visibilities[i][j], kernel->numchannels);
        }
    }
}

void getStationDelay(FxKernel *kernel, int antenna, int fftindex, double *meandelay, double *a, double *b) {
    double *interpolator = kernel->delays[antenna];
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

void unpack(FxKernel *kernel, u8 * inputdata, cf32 ** unpacked, int offset)
{
    if (kernel->nbits==2) {
        unpackReal2bit(inputdata, unpacked, offset, kernel->fftchannels);
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

void fringerotate(FxKernel *kernel, cf32 ** unpacked, f64 a, f64 b)
{
    int integerdelay;
    int status;

    // subtract off any integer delay (whole seconds) present (of course, this shouldn't be present).
    integerdelay = (int)b;
    b -= integerdelay;

    // Fill in the delay values, using a and b and the precomputed offsets
    status = vectorMulC_f64(kernel->subxoff, a, kernel->subxval, kernel->substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate, subval multiplication\n");
    status = vectorMulC_f64(kernel->stepxoff, a, kernel->stepxval, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate, stepval multiplication\n");
    status = vectorAddC_f64_I(b, kernel->subxval, kernel->substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate, subval addition!!!\n");

    // Turn delay into turns of phase by multiplying by the lo
    status = vectorMulC_f64(kernel->subxval, kernel->lofreq, kernel->subphase, kernel->substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate lofreq sub multiplication!!!\n");
    status = vectorMulC_f64(kernel->stepxval, kernel->lofreq, kernel->stepphase, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in linearinterpolate lofreq step multiplication!!!\n");
    if(kernel->fractionalLoFreq)
    {
        status = vectorAddC_f64_I((kernel->lofreq-(int)(kernel->lofreq))*(double)(integerdelay), kernel->subphase, kernel->substridesize);
        if(status != vecNoErr)
            fprintf(stderr, "Error in linearinterpolate lofreq non-integer freq addition!!!\n");
    }

    // Convert turns of phase into radians and bound into [0,2pi), then take sin/cos and assemble rotator vector
    for(int i=0;i<kernel->substridesize;i++)
    {
        kernel->subarg[i] = -TWO_PI*(kernel->subphase[i] - (int)(kernel->subphase[i]));
    }
    for(int i=0;i<kernel->stridesize;i++)
    {
        kernel->steparg[i] = -TWO_PI*(kernel->stepphase[i] - (int)(kernel->stepphase[i]));
    }
    status = vectorSinCos_f32(kernel->subarg, kernel->subsin, kernel->subcos, kernel->substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in sin/cos of sub rotate argument!!!\n");
    status = vectorSinCos_f32(kernel->steparg, kernel->stepsin, kernel->stepcos, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in sin/cos of step rotate argument!!!\n");
    status = vectorRealToComplex_f32(kernel->subcos, kernel->subsin, kernel->complexrotator, kernel->substridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error assembling sub into complex!!!\n");
    status = vectorRealToComplex_f32(kernel->stepcos, kernel->stepsin, kernel->stepcplx, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error assembling step into complex!!!\n");
    for(int i=1;i<kernel->stridesize;i++)
    {
        status = vectorMulC_cf32(kernel->complexrotator, kernel->stepcplx[i], &kernel->complexrotator[i*kernel->substridesize], kernel->substridesize);
        if(status != vecNoErr)
            fprintf(stderr, "Error doing the time-saving complex multiplication!!!\n");
    }

    // Actually apply the fringe rotation to each polarisation in turn
    for (int i=0;i<2;i++)
    {
        status = vectorMul_cf32_I(kernel->complexrotator, unpacked[i], kernel->fftchannels);
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

void fracSampleCorrect(FxKernel *kernel, cf32 ** channelised, f64 fracdelay)
{
    int status;

    // Create an array of phases for the fractional sample correction
    status = vectorMulC_f32(kernel->subchannelfreqs, fracdelay, kernel->subfracsamparg, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, arg generation (sub)!!! %d \n", status);

    status = vectorMulC_f32(kernel->stepchannelfreqs, fracdelay, kernel->stepfracsamparg, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, arg generation (step)!!! %d \n", status);

    //create the fractional sample correction array
    status = vectorSinCos_f32(kernel->subfracsamparg, kernel->subfracsampsin, kernel->subfracsampcos, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, sin/cos (sub)!!! %d \n", status);

    status = vectorSinCos_f32(kernel->stepfracsamparg, kernel->stepfracsampsin, kernel->stepfracsampcos, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, sin/cos (step)!!! %d \n", status);

    status = vectorRealToComplex_f32(kernel->subfracsampcos, kernel->subfracsampsin, kernel->fracsamprotator, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, real to complex (sub)!!! %d \n", status);

    status = vectorRealToComplex_f32(kernel->stepfracsampcos, kernel->stepfracsampsin, kernel->stepfracsampcplx, kernel->stridesize);
    if(status != vecNoErr)
        fprintf(stderr, "Error in frac sample correction, real to complex (step)!!! %d \n", status);
    for(int i=1;i<kernel->stridesize;i++)
    {
        status = vectorMulC_cf32(kernel->fracsamprotator, kernel->stepfracsampcplx[i], &(kernel->fracsamprotator[i*kernel->stridesize]), kernel->stridesize);
        if(status != vecNoErr)
            fprintf(stderr, "Error doing the time-saving complex multiplication in frac sample correction!!! \n");
    }

    // Apply the fractional sample correction to each polarisation in turn
    for(int i=0;i<2;i++)
    {
        status = vectorMul_cf32_I(kernel->fracsamprotator, channelised[i], kernel->numchannels);
        if(status != vecNoErr)
            fprintf(stderr, "Error in frac sample correction!!! %d \n", status);
    }
}

void conjChannels(FxKernel *kernel, cf32 ** channelised, cf32 ** conjchannels) {
    // To avoid calculating this multiple times, generate the complex conjugate of the channelised data
    vecStatus status;

    for (int i=0; i<2; i++) {
        status = vectorConj_cf32(channelised[i], conjchannels[i], kernel->numchannels); // Assumes USB and throws away 1 channel for real data
        if(status != vecNoErr) {
            fprintf(stderr, "Error calling vectorConj \n");
            exit(1);
        }
    }
}

void saveVisibilities(const char *outfile, FxKernel *kernel) {
    f32 ***amp, ***phase;

    FILE *fvis = fopen(outfile, "w");
    if (!fvis) {
        fprintf(stderr, "Error opening file for writing.\n");
        return;
    }

    amp = (f32***)malloc(kernel->nbaselines * sizeof(f32**));
    phase = (f32***)malloc(kernel->nbaselines * sizeof(f32**));
    if (!amp || !phase) {
        fprintf(stderr, "Memory allocation failed. Quitting.\n");
        return;
    }

    for (int i = 0; i < kernel->nbaselines; i++) {
        amp[i] = (f32**)malloc(4 * sizeof(f32*)); // 2 pols and crosspol
        phase[i] = (f32**)malloc(4 * sizeof(f32*)); // 2 pols and crosspol
        if (!amp[i] || !phase[i]) {
            fprintf(stderr, "Memory allocation failed. Quitting.\n");
            return;
        }
        for (int j = 0; j < 4; j++) {
            amp[i][j] = vectorAlloc_f32(kernel->numchannels);
            phase[i][j] = vectorAlloc_f32(kernel->numchannels);
            if (!amp[i][j] || !phase[i][j]) {
                fprintf(stderr, "Memory allocation failed. Quitting.\n");
                return;
            }
            vectorMagnitude_cf32(kernel->visibilities[i][j], amp[i][j], kernel->numchannels);
            vectorPhase_cf32(kernel->visibilities[i][j], phase[i][j], kernel->numchannels);
        }
    }

    for (int c = 0; c < kernel->numchannels; c++) {
        fprintf(fvis, "%5d %11.6f", c, (c + 0.5) / kernel->numchannels * kernel->bandwidth / 1e6);
        for (int i = 0; i < kernel->nbaselines; i++) {
            for (int j = 0; j < 4; j++) {
                fprintf(fvis, " %11.6f %11.6f", kernel->visibilities[i][j][c].re, kernel->visibilities[i][j][c].im);
                fprintf(fvis, " %11.6f %10.6f", amp[i][j][c], phase[i][j][c]);
            }
        }
        fprintf(fvis, "\n");
    }
    fclose(fvis);

    for (int i = 0; i < kernel->nbaselines; i++) {
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

void freeFxKernel(FxKernel *kernel) {
    if (kernel == NULL) return;

    // Free input data
    if (kernel->inputdata) {
        for (int i = 0; i < kernel->numantennas; i++) {
            free(kernel->inputdata[i]);
        }
        free(kernel->inputdata);
    }

    // Free unpacked data
    if (kernel->unpacked) {
        for (int i = 0; i < kernel->numantennas; i++) {
            for (int j = 0; j < 2; j++) {
                ippFree(kernel->unpacked[i][j]);
            }
            free(kernel->unpacked[i]);
        }
        free(kernel->unpacked);
    }

    // Free channelised data
    if (kernel->channelised) {
        for (int i = 0; i < kernel->numantennas; i++) {
            for (int j = 0; j < 2; j++) {
                ippFree(kernel->channelised[i][j]);
            }
            free(kernel->channelised[i]);
        }
        free(kernel->channelised);
    }

    // Free conjchannels data
    if (kernel->conjchannels) {
        for (int i = 0; i < kernel->numantennas; i++) {
            for (int j = 0; j < 2; j++) {
                ippFree(kernel->conjchannels[i][j]);
            }
            free(kernel->conjchannels[i]);
        }
        free(kernel->conjchannels);
    }

    // Free visibilities data
    if (kernel->visibilities) {
        for (int i = 0; i < kernel->nbaselines; i++) {
            for (int j = 0; j < 4; j++) {
                ippFree(kernel->visibilities[i][j]);
            }
            free(kernel->visibilities[i]);
        }
        free(kernel->visibilities);
    }

    // Free delays data
    if (kernel->delays) {
        for (int i = 0; i < kernel->numantennas; i++) {
            free(kernel->delays[i]);
        }
        free(kernel->delays);
    }

    // Free other internal arrays
    ippFree(kernel->subtoff);
    ippFree(kernel->subtval);
    ippFree(kernel->subxoff);
    ippFree(kernel->subxval);
    ippFree(kernel->subphase);
    ippFree(kernel->subarg);
    ippFree(kernel->subsin);
    ippFree(kernel->subcos);
    ippFree(kernel->steptoff);
    ippFree(kernel->steptval);
    ippFree(kernel->stepxoff);
    ippFree(kernel->stepxval);
    ippFree(kernel->stepphase);
    ippFree(kernel->steparg);
    ippFree(kernel->stepsin);
    ippFree(kernel->stepcos);
    ippFree(kernel->stepcplx);
    ippFree(kernel->complexrotator);
    ippFree(kernel->subfracsamparg);
    ippFree(kernel->subfracsampsin);
    ippFree(kernel->subfracsampcos);
    ippFree(kernel->subchannelfreqs);
    ippFree(kernel->stepfracsamparg);
    ippFree(kernel->stepfracsampsin);
    ippFree(kernel->stepfracsampcos);
    ippFree(kernel->stepfracsampcplx);
    ippFree(kernel->stepchannelfreqs);
    ippFree(kernel->fracsamprotator);

    // Free FFTs
    ippFree(kernel->fftbuffer);
    if (kernel->pFFTSpecC) {
        ippFree(kernel->pFFTSpecC);
    }

    // Free validity/weights
    free(kernel->antValid);
    free(kernel->baselineCount);
}

void saveLog(long long diff_ms, char *starttimestring) {
    // Open the file for writing
    FILE *ftiming = fopen("runtime.fxkernel.log", "w");
    if (ftiming == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        // Handle error as needed
    }

// Write to the file
    fprintf(ftiming, "Run time was %lld milliseconds\n", diff_ms);
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
}