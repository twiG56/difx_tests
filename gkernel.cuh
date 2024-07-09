//
// Created by emichel on 18/06/24.
//

#ifndef TESTS_DIFX_GKERNEL_CUH
#define TESTS_DIFX_GKERNEL_CUH

#include <cuComplex.h>
#include <cufft.h>
#include "cuda_runtime_api.h"
#include <stdbool.h>
#include <stdio.h>

#define COMPLEX cuComplex
#define MAKECOMPLEX(x, y)  make_cuFloatComplex(x,y)
#define HALF2FLOAT2(x)    x
#define HALF2FLOAT(x)     x

__host__ __device__ static __inline__ int
sampIdx(int antenna, int pol, int sample, int stride) {
    const int num_pols = 2;

    return (antenna * num_pols + pol) * stride + sample;
}

__host__ __device__ static __inline__ int
antIdx(int antenna, int pol, int channel, int stride) {
    const int num_pols = 2;

    return (antenna * num_pols + pol) * stride + channel;
}

__host__ __device__ static __inline__ int
accumIdx(int baseline, int product, int channel, int stride) {
    const int num_products = 4;

    return (baseline * num_products + product) * stride + channel;
}

__global__ void FracSampleCorrection(COMPLEX *ant, float *fractionalDelayValues,
                                     int numchannels, int fftsamples, int numffts, int subintsamples);

__global__ void FracSampleCorrection1D(COMPLEX *ant, float *fractionalDelayValues,
                                       int numchannels, int fftsamples, int numffts, int subintsamples);

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        exit(code);
    }
}


#define CudaCheckError()  __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaCheckError(const char *file, const int line) {
    enum cudaError err = cudaGetLastError();
    if (cudaSuccess != err) {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                file, line - 1, cudaGetErrorString(err));
        exit(err);
    }
}

void init_2bitLevels();

void allocDataGPU(int8_t ****packedData, COMPLEX ***unpackedData,
                  COMPLEX ***channelisedData, cuComplex ***baselineData,
                  float ***rotationPhaseInfo, float ***fractionalSampleDelays, int ***sampleShifts,
                  double ***gpuDelays, int numantenna, int subintsamples, int nbit, int nPol,
                  bool iscomplex, int nchan, int numffts, int parallelAccum, int num_streams, int num_packed_bytes);

__global__ void
calculateDelaysAndPhases(double *gpuDelays, double lo, double sampletime, int fftsamples, int fftchannels,
                         int samplegranularity, float *rotationPhaseInfo, int *sampleShifts,
                         float *fractionalSampleDelays);

__global__ void
calculateDelaysAndPhases1D(double *gpuDelays, double lo, double sampletime, int fftsamples, int fftchannels,
                           int samplegranularity, float *rotationPhaseInfo, int *sampleShifts,
                           float *fractionalSampleDelays);


__global__ void unpack2bit_2chan_rotate(COMPLEX *dest, const int8_t *src, float *rotVec, const int32_t *shifts,
                                        const int32_t fftsamples);

__global__ void unpack8bitcomplex_2chan_rotate(COMPLEX *dest, const int8_t *src, float *rotVec, const int32_t *shifts,
                                               const int32_t fftsamples);

__global__ void CCAH3(cuComplex *accum, const COMPLEX *ants, int nant, int nfft, int nchan, int fftwidth);

__global__ void CCAH2(cuComplex *accum, const COMPLEX *ants, int nant, int nfft, int nchan, int fftwidth);

__global__ void CrossCorrAccumHoriz(cuComplex *accum, const COMPLEX *ants, int nant, int nfft, int nchan, int fftwidth);

void
saveVisibilities(const char *outfile, cuComplex *baselines, int nbaseline, int nchan, int stride, double bandwidth);

void copyAndIterate(cuComplex **unpackedData, int s, int numantenna, int nPol, int subintsamples);

void copyAndIterateBaseline(cuComplex **unpackedData, int nbaseline, int nchan, int parallelAccum, int subintsamples);

#endif //TESTS_DIFX_GKERNEL_CUH
