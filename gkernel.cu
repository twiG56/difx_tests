//
// Created by emichel on 18/06/24.
//

#include "gkernel.cuh"

using namespace std;

__constant__ float kLevels_2bit[4];

void init_2bitLevels() {
    static const float HiMag = 3.3359;  // Optimal value
    const float lut4level[4] = {-HiMag, -1.0, 1.0, HiMag};
    gpuErrchk(cudaMemcpyToSymbol(kLevels_2bit, lut4level, 4*sizeof(float), 0, cudaMemcpyHostToDevice));
}

// Rotate inplace a complex number by theta (radians)
__device__ static __inline__ void cuRotatePhase (COMPLEX &x, float theta)
{
    float cs, sn;
    sincosf(theta, &sn, &cs);

#ifdef USEHALF
    float2 y = __half22float2(x);
  float px = y.x * cs - y.y * sn;
  float py = y.x * sn + y.y * cs;
#else
    float px = x.x * cs - x.y * sn;
    float py = x.x * sn + x.y * cs;
#endif
    x = MAKECOMPLEX(px, py);
    return;
}

__global__ void FracSampleCorrection(COMPLEX *ant, float *fractionalDelayValues,
                                     int numchannels, int fftsamples, int numffts, int subintsamples) {
    size_t ichan = threadIdx.x + blockIdx.x * blockDim.x;
    size_t ifft = blockIdx.y;
    size_t iant = blockIdx.z;

    // phase and slope for this FFT
    float dslope = fractionalDelayValues[iant*numffts + ifft];
    float theta = ichan*dslope;
    cuRotatePhase(ant[sampIdx(iant, 0, ichan+ifft*fftsamples, subintsamples)], theta);
    cuRotatePhase(ant[sampIdx(iant, 1, ichan+ifft*fftsamples, subintsamples)], theta);
}

__global__ void FracSampleCorrection1D(COMPLEX *ant, float *fractionalDelayValues,
                                     int numchannels, int fftsamples, int numffts, int subintsamples) {
    size_t ichan = threadIdx.x;
    size_t combinedIndex = blockIdx.x;
    size_t ifft = combinedIndex % numffts;
    size_t iant = combinedIndex / numffts;

    // phase and slope for this FFT
    float dslope = fractionalDelayValues[iant*numffts + ifft];
    float theta = ichan*dslope;
    cuRotatePhase(ant[sampIdx(iant, 0, ichan+ifft*fftsamples, subintsamples)], theta);
    cuRotatePhase(ant[sampIdx(iant, 1, ichan+ifft*fftsamples, subintsamples)], theta);
}

void allocDataGPU(int8_t ****packedData, COMPLEX ***unpackedData,
                  COMPLEX ***channelisedData, cuComplex ***baselineData,
                  float ***rotationPhaseInfo, float ***fractionalSampleDelays, int ***sampleShifts,
                  double ***gpuDelays, int numantenna, int subintsamples, int nbit, int nPol,
                  bool iscomplex, int nchan, int numffts, int parallelAccum, int num_streams, int num_packed_bytes) {

    unsigned long long GPUalloc = 0;

    *packedData = (int8_t ***)malloc(num_streams * sizeof(int8_t **));
    *unpackedData = (COMPLEX **)malloc(num_streams * sizeof(COMPLEX *));
    *channelisedData = (COMPLEX **)malloc(num_streams * sizeof(COMPLEX *));
    *baselineData = (cuComplex **)malloc(num_streams * sizeof(cuComplex *));
    *rotationPhaseInfo = (float **)malloc(num_streams * sizeof(float *));
    *fractionalSampleDelays = (float **)malloc(num_streams * sizeof(float *));
    *sampleShifts = (int **)malloc(num_streams * sizeof(int *));
    *gpuDelays = (double **)malloc(num_streams * sizeof(double *));

    // Unpacked data
    printf("Alloc %d complex output values per baseline\n", nchan * parallelAccum);
    for (int s=0; s<num_streams; s++) {
        (*packedData)[s] = (int8_t **)malloc(numantenna * sizeof(int8_t *));
        if ((*packedData)[s] == NULL) {
            // Handle allocation failure
            fprintf(stderr, "Memory allocation failed for packedData[%d]\n", s);
            exit(1);
        }
        for (int i=0; i<numantenna; i++) {
            gpuErrchk(cudaMalloc(&((*packedData)[s][i]), num_packed_bytes));

            GPUalloc += num_packed_bytes;
        }

        gpuErrchk(cudaMalloc(&(*unpackedData)[s], numantenna*nPol*subintsamples*sizeof(cuComplex)));
        //printf("Unpacked Data size : %d \n", numantenna*nPol*subintsamples*sizeof(cuComplex));
        GPUalloc += numantenna*nPol*subintsamples*sizeof(cuComplex);

        // FFT output
        gpuErrchk(cudaMalloc(&(*channelisedData)[s], numantenna*nPol*subintsamples*sizeof(cuComplex)));
        //printf("Channelised Data size : %d \n", numantenna*nPol*subintsamples*sizeof(cuComplex));
        GPUalloc += numantenna*nPol*subintsamples*sizeof(cuComplex);

        // Baseline visibilities
        int nbaseline = numantenna*(numantenna-1)/2;
        if (!iscomplex) subintsamples /= 2;
        gpuErrchk(cudaMalloc(&(*baselineData)[s], nbaseline*4*nchan*parallelAccum*sizeof(cuComplex)));
        //printf("Baseline Data size : %d \n", nbaseline*4*nchan*parallelAccum*sizeof(cuComplex));
        GPUalloc += nbaseline*4*nchan*parallelAccum*sizeof(cuComplex);

        // Fringe rotation vector (will contain starting phase and phase increment for every FFT of every antenna)
        gpuErrchk(cudaMalloc(&(*rotationPhaseInfo)[s], numantenna*numffts*2*sizeof(float)));
        //printf("rotationPhaseInfo size : %d \n", numantenna*numffts*2*sizeof(float));
        GPUalloc += numantenna*numffts*2*sizeof(float);

        // Fractional sample delay vector (will contain midpoint fractional sample delay [in units of radians per channel!]
        // for every FFT of every antenna)
        gpuErrchk(cudaMalloc(&(*fractionalSampleDelays)[s], numantenna*numffts*sizeof(float)));
        //printf("fractionalSampleDelays size : %d \n", numantenna*numffts*sizeof(float));
        GPUalloc += numantenna*numffts*sizeof(float);

        // Sample shifts vector (will contain the integer sample shift relative to nominal FFT start for every FFT of every antenna)
        gpuErrchk(cudaMalloc(&(*sampleShifts)[s], numantenna*numffts*sizeof(int)));
        //printf("sampleShifts size : %d \n", numantenna*numffts*sizeof(int));
        GPUalloc += numantenna*numffts*sizeof(int);

        // Delay information vectors
        gpuErrchk(cudaMalloc(&(*gpuDelays)[s], numantenna*4*sizeof(double)));
        //printf("gpuDelays size : %d \n", numantenna*4*sizeof(double));
        GPUalloc += numantenna*4*sizeof(double);
    }

    printf("Allocated %.2f Mb on GPU\n", GPUalloc / 1e6);
}

/* Calculate the starting fringe rotation phase and phase increment for each FFT of each antenna, and the fractional sample error */
__global__ void calculateDelaysAndPhases(double * gpuDelays, double lo, double sampletime, int fftsamples, int fftchannels, int samplegranularity, float * rotationPhaseInfo, int *sampleShifts, float* fractionalSampleDelays)
{
    size_t ifft = threadIdx.x + blockIdx.x * blockDim.x;
    size_t iant = blockIdx.y;
    int numffts = blockDim.x * gridDim.x;
    double meandelay, deltadelay, netdelaysamples_f, startphase;
    double d0, d1, d2, a, b;
    double * interpolator = &(gpuDelays[iant*4]);
    double filestartoffset = gpuDelays[iant*4+3];
    float fractionaldelay;
    int netdelaysamples;

    // evaluate the delay for the given FFT of the given antenna

    // calculate values at the beginning, middle, and end of this FFT
    d0 = interpolator[0]*ifft*ifft + interpolator[1]*ifft + interpolator[2];
    d1 = interpolator[0]*(ifft+0.5)*(ifft+0.5) + interpolator[1]*(ifft+0.5) + interpolator[2];
    d2 = interpolator[0]*(ifft+1.0)*(ifft+1.0) + interpolator[1]*(ifft+1.0) + interpolator[2];

    // use these to calculate a linear interpolator across the FFT, as well as a mean value
    a = d2-d0; //this is the delay gradient across this FFT
    b = d0 + (d1 - (a*0.5 + d0))/3.0; //this is the delay at the start of the FFT
    meandelay = a*0.5 + b; //this is the delay in the middle of the FFT
    deltadelay = a / fftsamples; // this is the change in delay per sample across this FFT window

    netdelaysamples_f = (meandelay - filestartoffset) / sampletime;
    netdelaysamples = __double2int_rn(netdelaysamples_f/samplegranularity) * samplegranularity;

    // Save the integer number of sample shifts
    sampleShifts[iant*numffts + ifft] = netdelaysamples;

    // Save the fractional delay
    fractionaldelay = (float)(-(netdelaysamples_f - netdelaysamples)*2*M_PI/fftsamples);  // radians per FFT channel
    fractionalSampleDelays[iant*numffts + ifft] = fractionaldelay;

    // set the fringe rotation phase for the first sample of a given FFT of a given antenna
    startphase = b*lo;
    rotationPhaseInfo[iant*numffts*2 + ifft*2] = (float)(startphase - int(startphase))*2*M_PI;
    rotationPhaseInfo[iant*numffts*2 + ifft*2 + 1] = (float)(deltadelay * lo)*2*M_PI;
}

__global__ void calculateDelaysAndPhases1D(double * gpuDelays, double lo, double sampletime, int fftsamples, int fftchannels, int samplegranularity, float * rotationPhaseInfo, int *sampleShifts, float* fractionalSampleDelays)
{
    size_t globalIdx = threadIdx.x + blockIdx.x * blockDim.x;
    size_t totalThreads = gridDim.x * blockDim.x;

    size_t numffts = totalThreads / (4 * 15);
    size_t iant = (globalIdx / numffts) % 4;
    size_t executionIndex = globalIdx / (numffts * 4);
    size_t ifft = globalIdx % numffts;

    double meandelay, deltadelay, netdelaysamples_f, startphase;
    double d0, d1, d2, a, b;
    double *interpolator = &(gpuDelays[iant * 4]);
    double filestartoffset = gpuDelays[iant * 4 + 3];
    float fractionaldelay;
    int netdelaysamples;

    // Evaluate the delay for the given FFT of the given antenna

    // Calculate values at the beginning, middle, and end of this FFT
    d0 = interpolator[0] * ifft * ifft + interpolator[1] * ifft + interpolator[2];
    d1 = interpolator[0] * (ifft + 0.5) * (ifft + 0.5) + interpolator[1] * (ifft + 0.5) + interpolator[2];
    d2 = interpolator[0] * (ifft + 1.0) * (ifft + 1.0) + interpolator[1] * (ifft + 1.0) + interpolator[2];

    // Use these to calculate a linear interpolator across the FFT, as well as a mean value
    a = d2 - d0; // This is the delay gradient across this FFT
    b = d0 + (d1 - (a * 0.5 + d0)) / 3.0; // This is the delay at the start of the FFT
    meandelay = a * 0.5 + b; // This is the delay in the middle of the FFT
    deltadelay = a / fftsamples; // This is the change in delay per sample across this FFT window

    netdelaysamples_f = (meandelay - filestartoffset) / sampletime;
    netdelaysamples = __double2int_rn(netdelaysamples_f / samplegranularity) * samplegranularity;

    // Save the integer number of sample shifts
    sampleShifts[iant * fftchannels + ifft] = netdelaysamples;

    // Save the fractional delay
    fractionaldelay = (float)(-(netdelaysamples_f - netdelaysamples) * 2 * M_PI / fftsamples); // Radians per FFT channel
    fractionalSampleDelays[iant * fftchannels + ifft] = fractionaldelay;

    // Set the fringe rotation phase for the first sample of a given FFT of a given antenna
    startphase = b * lo;
    rotationPhaseInfo[iant * fftchannels * 2 + ifft * 2] = (float)(startphase - int(startphase)) * 2 * M_PI;
    rotationPhaseInfo[iant * fftchannels * 2 + ifft * 2 + 1] = (float)(deltadelay * lo) * 2 * M_PI;
}

// Rotate a complex number by theta (radians)
__device__ static __inline__ void cuRotatePhase3 (float x, COMPLEX &y, float sinA, float cosA)
{
    y = MAKECOMPLEX(x * cosA, x * sinA);
    return;
}

__global__ void unpack2bit_2chan_rotate(COMPLEX *dest, const int8_t *src, float *rotVec, const int32_t *shifts, const int32_t fftsamples) {
    // static const float HiMag = 3.3359;  // Optimal value
    // const float levels_2bit[4] = {-HiMag, -1.0, 1.0, HiMag};
    const size_t isample = 2*(blockDim.x * blockIdx.x + threadIdx.x);
    const size_t ifft = blockIdx.y;
    const size_t osample = isample + ifft*fftsamples;
    int subintsamples = fftsamples * gridDim.y;

    // Try to Fix
    size_t idx = ((isample - shifts[ifft])/2); // FIXME: may lead to memory access outside src[] bounds, see with 'cuda-memcheck ./benchmark_gxkernel'
    // And of try to fix

    int8_t src_i = src[idx]; // Here I am just loading src into local memory to
    // reduce the number of reads from global memory

    // I have just changed the order of the writes made to dest
    // In theory this should reduce the number of write operations made
    // I have also implemented the use of constant memory for the levels_2bit
    // array

    float samp0 = kLevels_2bit[src_i&0x3];
    float samp1 = kLevels_2bit[(src_i>>4)&0x3];
    float samp2 = kLevels_2bit[(src_i>>2)&0x3];
    float samp3 = kLevels_2bit[(src_i>>6)&0x3];

    // phase and slope for this FFT
    float p0 = rotVec[ifft*2];
    float p1 = rotVec[ifft*2+1];
    float theta0 = -p0 - isample*p1;
    float theta1 = -p0 - (isample+1)*p1;

    float sinT0, cosT0, sinT1, cosT1;
    sincosf(theta0, &sinT0, &cosT0);
    sincosf(theta1, &sinT1, &cosT1);
    cuRotatePhase3(samp0, dest[(osample)], sinT0, cosT0);
    cuRotatePhase3(samp1, dest[(osample+1)], sinT1, cosT1);
    cuRotatePhase3(samp2, dest[(subintsamples + osample)], sinT0, cosT0);
    cuRotatePhase3(samp3, dest[(subintsamples + osample+1)], sinT1, cosT1);
}

// Rotate a complex number by theta (radians)
__device__ static __inline__ void cuRotatePhase4 (cuComplex x, COMPLEX &y, float sinA, float cosA)
{
    y = MAKECOMPLEX(x.x * cosA - x.y * sinA, x.x * sinA + x.y * cosA);
    return;
}

__global__ void unpack8bitcomplex_2chan_rotate(COMPLEX *dest, const int8_t *src, float *rotVec, const int32_t *shifts, const int32_t fftsamples) {
    const size_t isamp = (blockDim.x * blockIdx.x + threadIdx.x); //This can go from 0 ... fftsamples*2 (i.e., number of samples in an FFT * 2 channels)
    const size_t ifft = blockIdx.y;
    int subintsamples = fftsamples * gridDim.y;

    int ibyte = isamp*2; // 2 bytes per complex sample
    int pol = isamp % 2;
    int osamp = isamp/2 + pol*subintsamples;

    cuComplex samp = make_cuFloatComplex(src[ibyte - shifts[ifft]*4], src[ibyte - shifts[ifft]*4 + 1]);

    // phase and slope for this FFT
    float p0 = rotVec[ifft*2];
    float p1 = rotVec[ifft*2+1];
    float theta = p0 + isamp*p1;

    float sinT, cosT;
    sincosf(theta, &sinT, &cosT);
    cuRotatePhase4(samp, dest[ifft*fftsamples + osamp], sinT, cosT);
}

__global__ void CCAH3(cuComplex *accum, const COMPLEX *ants, int nant, int nfft, int nchan, int fftwidth) {
    int t = threadIdx.x+blockIdx.x*blockDim.x;
    if (t>=nchan) return;

    // Assuming nPol ==2 !!!!

    // blockIdx.y: index of first vector (antennaindex)
    // blockIdx.z: index delta to second vector, minus 1.
    int ant1 = blockIdx.y;
    int ant2 = ant1 + blockIdx.z + 1;

    if (ant2>=nant)  return;

    // index into output vector blocks: = (j-i-1) + n-1 + ... + n-i
    int b = ant1*nant-ant1*(ant1+1)/2 + -ant1 + ant2-1;

    int s = nfft*fftwidth;

    const COMPLEX* iv = ants+ant1*s*2+t;
    const COMPLEX* jv = ants+ant2*s*2+t;

    COMPLEX u1 = iv[0];
    COMPLEX v1 = jv[0];
    COMPLEX u2 = iv[s];
    COMPLEX v2 = jv[s];
    cuComplex a1;
    cuComplex a2;
    cuComplex a3;
    cuComplex a4;
    a1.x = (u1.x*v1.x + u1.y*v1.y);

    a1.y = u1.y*v1.x - u1.x*v1.y;
    a2.x = u1.x*v2.x + u1.y*v2.y;
    a2.y = u1.y*v2.x - u1.x*v2.y;
    a3.x = u2.x*v1.x + u2.y*v1.y;
    a3.y = u2.y*v1.x - u2.x*v1.y;
    a4.x = u2.x*v2.x + u2.y*v2.y;
    a4.y = u2.y*v2.x - u2.x*v2.y;

    for (int k = fftwidth; k<s; k += fftwidth) {
        u1 = iv[k];
        v1 = jv[k];
        u2 = iv[k+s];
        v2 = jv[k+s];

        a1.x += HALF2FLOAT(u1.x*v1.x + u1.y*v1.y);
        a1.y += HALF2FLOAT(u1.y*v1.x - u1.x*v1.y);
        a2.x += HALF2FLOAT(u1.x*v2.x + u1.y*v2.y);
        a2.y += HALF2FLOAT(u1.y*v2.x - u1.x*v2.y);
        a3.x += HALF2FLOAT(u2.x*v1.x + u2.y*v1.y);
        a3.y += HALF2FLOAT(u2.y*v1.x - u2.x*v1.y);
        a4.x += HALF2FLOAT(u2.x*v2.x + u2.y*v2.y);
        a4.y += HALF2FLOAT(u2.y*v2.x - u2.x*v2.y);
    }

    a1.x /= nfft;
    a1.y /= nfft;
    a2.x /= nfft;
    a2.y /= nfft;
    a3.x /= nfft;
    a3.y /= nfft;
    a4.x /= nfft;
    a4.y /= nfft;
    accum[4*b*nchan+t] = a1;
    accum[(4*b+1)*nchan+t] = a2;
    accum[(4*b+2)*nchan+t] = a3;
    accum[(4*b+3)*nchan+t] = a4;
}

// Function to compute the phase angle of a cuComplex number
float cuCargf(cuComplex z) {
    return atan2f(cuCimagf(z), cuCrealf(z));
}

void saveVisibilities(const char *outfile, cuComplex *baselines, int nbaseline, int nchan, int stride, double bandwidth) {
    cuComplex **vis;
    FILE *fvis = fopen(outfile, "w");
    if (fvis == NULL) {
        fprintf(stderr, "Error opening file %s\n", outfile);
        exit(1);
    }

    // Copy final visibilities back to CPU
    vis = (cuComplex **)malloc(nbaseline * 4 * sizeof(cuComplex *));
    for (int i = 0; i < nbaseline * 4; i++) {
        vis[i] = (cuComplex *)malloc(nchan * sizeof(cuComplex));
        gpuErrchk(cudaMemcpy(vis[i], &baselines[i * stride], nchan * sizeof(cuComplex), cudaMemcpyDeviceToHost));
    }

    printf("Test : %e \n", vis[0][0].x);
    printf("Test : %e \n", vis[0][0].y);

    for (int c = 0; c < nchan; c++) {
        fprintf(fvis, "%5d %11.6f", c, (c + 0.5) / nchan * bandwidth / 1e6);
        for (int i = 0; i < nbaseline * 4; i++) {
            fprintf(fvis, " %11.6f %11.6f %11.6f %10.6f",
                    cuCrealf(vis[i][c]), cuCimagf(vis[i][c]),
                    cuCabsf(vis[i][c]), cuCargf(vis[i][c]));
        }
        fprintf(fvis, "\n");
    }
    fclose(fvis);

    for (int i = 0; i < nbaseline * 4; i++) {
        free(vis[i]);
    }
    free(vis);
}

__global__ void CCAH2(cuComplex *accum, const COMPLEX *ants, int nant, int nfft, int nchan, int fftwidth) {
    int t = threadIdx.x+blockIdx.x*blockDim.x;
    if (t>=nchan) return;

    // blockIdx.y: index of first vector (2*antennaindex+polindex)
    // blockIdx.z: index delta to second vector, minus 1.
    int ii = blockIdx.y;
    int ij = blockIdx.z;

    ij += ii+1;

    int ai = ii/2;
    int aj = ij/2;


    if (ai>=aj || ai>=nant || aj>=nant) {
        return;
    }
    int pi = ii%2;
    int pj = ij%2;

    // index into output vector blocks: = (j-i-1) + n-1 + ... + n-i
    int b = 4*(ai*nant-ai*(ai+1)/2 + aj-ai-1)+2*pi+pj;

    int s = nfft*fftwidth;

    const COMPLEX* iv = ants+ii*s+t;
    const COMPLEX* jv = ants+ij*s+t;

    float2 u = HALF2FLOAT2(iv[0]);
    float2 v = HALF2FLOAT2(jv[0]);
    float2 a;
    a.x = u.x*v.x + u.y*v.y;
    a.y = u.y*v.x - u.x*v.y;

    for (int k = fftwidth; k<s; k += fftwidth) {
        u = HALF2FLOAT2(iv[k]);
        v = HALF2FLOAT2(jv[k]);

        a.x += u.x*v.x + u.y*v.y;
        a.y += u.y*v.x - u.x*v.y;
    }

    a.x /= nfft;
    a.y /= nfft;
    accum[b*nchan+t] = a;
}

__global__ void CrossCorrAccumHoriz(cuComplex *accum, const COMPLEX *ants, int nant, int nfft, int nchan, int fftwidth) {
    int t = threadIdx.x+blockIdx.x*blockDim.x;
    if (t>=nchan) return;

    // input vector indices in block .y and .z
    int i = blockIdx.y;
    int j = blockIdx.z;
    j += i+1;

    if (i>=nant || j>=nant) return;

    // index into output vectors: = (j-i-1) + n-1 + ... + n-i
    int b = i*nant-i*(i+1)/2 + j-i-1;

    int s = nfft*fftwidth;

    for (int pi = 0; pi<2; ++pi) {
        for (int pj = 0; pj<2; ++pj) {
            const COMPLEX* iv = &ants[antIdx(i, pi, t, s)];
            const COMPLEX* jv = &ants[antIdx(j, pj, t, s)];

            COMPLEX u = iv[0];
            COMPLEX v = jv[0];
            COMPLEX a;
            a.x = u.x*v.x + u.y*v.y;
            a.y = u.y*v.x - u.x*v.y;

            for (int k = fftwidth; k<s; k += fftwidth) {
                u = iv[k];
                v = jv[k];

                a.x += u.x*v.x + u.y*v.y;
                a.y += u.y*v.x - u.x*v.y;
            }

            a.x /= nfft;
            a.y /= nfft;
            accum[accumIdx(b, pi*2+pj, t, nchan)] = HALF2FLOAT2(a);
        }
    }
}


// TODO DELETE

void copyAndIterate(cuComplex** unpackedData, int s, int numantenna, int nPol, int subintsamples) {
    // Calculate the size of the data
    size_t dataSize = numantenna * nPol * subintsamples * sizeof(cuComplex);

    // Allocate memory on the CPU
    cuComplex* cpuData = (cuComplex*)malloc(dataSize);
    if (cpuData == NULL) {
        printf("Failed to allocate memory on the CPU \n");
        return;
    }

    // Copy data from GPU to CPU
    gpuErrchk(cudaMemcpy(cpuData, &(*unpackedData)[s], dataSize, cudaMemcpyDeviceToHost));

    // Iterate over the data on the CPU
//    for (int i = 0; i < numantenna; ++i) {
//        unpackedData[s][2 * i * subintsamples];
//        printf
//    }
    printf("test : %e\n", cpuData[2 * 0 * subintsamples].x);
    printf("test : %e\n", cpuData[2 * 0 * subintsamples + 1].x);

    printf("test : %e\n", cpuData[2 * 1 * subintsamples].x);
    printf("test : %e\n", cpuData[2 * 1 * subintsamples + 1].x);

    printf("test : %e\n", cpuData[2 * 1 * subintsamples + 3840].x);
    printf("test : %e\n", cpuData[2 * 1 * subintsamples + 1 +3840].x);

    // Free the CPU memory
    free(cpuData);
}

void copyAndIterateBaseline(cuComplex** unpackedData, int nbaseline, int nchan, int parallelAccum, int subintsamples) {
    // Calculate the size of the data
    size_t dataSize = nbaseline*4*nchan*parallelAccum*sizeof(cuComplex);

    // Allocate memory on the CPU
    cuComplex* cpuData = (cuComplex*)malloc(dataSize);
    if (cpuData == NULL) {
        printf("Failed to allocate memory on the CPU \n");
        return;
    }

    // Copy data from GPU to CPU
    gpuErrchk(cudaMemcpy(cpuData, &(*unpackedData)[0], dataSize, cudaMemcpyDeviceToHost));

    // Iterate over the data on the CPU
//    for (int i = 0; i < numantenna; ++i) {
//        unpackedData[s][2 * i * subintsamples];
//        printf
//    }
    printf("test : %e\n", cpuData[0].x);
    printf("test : %e\n", cpuData[1].x);
    printf("test : %e\n", cpuData[2].x);
    printf("test : %e\n", cpuData[3].x);


    // Free the CPU memory
    free(cpuData);
}

// TODO END DELETE