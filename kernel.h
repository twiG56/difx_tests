//
// Created by emichel on 01/05/24.
//

#ifndef TESTS_DIFX_KERNEL_H
#define TESTS_DIFX_KERNEL_H

#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <string.h>
#include <ipps.h>
#include <ippvm.h>
#include <ippcore.h>
#include "vectordefs.h"
#include <math.h>

// Structure definition
typedef struct {
    u8 ** inputdata; /**< input data array */

    cf32 *** unpacked; /**< unpacked data (fringe rotation is performed in-place here) */

    cf32 *** channelised; /**< channelised data */
    cf32 *** conjchannels; /**< conjugated values */

    // Validity/weights
    int *antValid; /**< checks if there is good data for the given antenna at a given time (if yes, it will accumulate; if not, if will leave this data out) */
    int *baselineCount; /**< counter incrementing number of baselines with successfully obtained cross correlations and accumulations */

    cf32 *** visibilities; /**< output data array */

    double ** delays; /**< delay polynomial for each antenna. Put in the time in units of FFT lengths since start of subintegration, get back delay in seconds. */

    double * filestartoffsets; /**< Offset for each antenna file from the nominal start time of the subintegration. In seconds. */

    // internal arrays
    f64 * subtoff;
    f64 * subtval;
    f64 * subxoff;
    f64 * subxval;
    f64 * subphase;
    f32 * subarg;
    f32 * subsin;
    f32 * subcos;

    f64 * steptoff;
    f64 * steptval;
    f64 * stepxoff;
    f64 * stepxval;
    f64 * stepphase;
    f32 * steparg;
    f32 * stepsin;
    f32 * stepcos;
    cf32 * stepcplx;
    cf32 * complexrotator;

    f32 * subfracsamparg;
    f32 * subfracsampsin;
    f32 * subfracsampcos;
    f32 * subchannelfreqs;

    f32 * stepfracsamparg;
    f32 * stepfracsampsin;
    f32 * stepfracsampcos;
    cf32 * stepfracsampcplx;
    f32 * stepchannelfreqs;
    cf32 * fracsamprotator;

    // FFTs
    u8 * fftbuffer;
    vecFFTSpecC_cf32 * pFFTSpecC;

    // other constants
    int numantennas; /**< number of antennas in the dataset */
    int numchannels; /**< the number of channels to create during channelisation; each channel is a unique frequency point */
    int nbits;  /**< Number of bits for voltage samples */
    int numffts;  /**< the number of FFTs computed; i.e., the length of a subint */
    int fftchannels; /**< length of an FFT; 2*nchan for real data, and nchan for complex */
    int nbaselines; /**< Number of baselines (nant*(nant-1)/2) */
    int stridesize; /**< used for the time-saving complex multiplications */
    int substridesize; /**< used for the time-saving complex multiplications. Equal to stridesize for complex data, or 2x stridesize for real data */
    double lofreq; /**< local oscillator frequency; in Hz */
    double bandwidth; /**< in Hz */
    double sampletime; /**< 1/(2*bandwidth); in seconds */
    int fractionalLoFreq; /**< if true, means we need to do an extra multiplication in fringe rotation phase calculation */
    int iscomplex;  /**< Is the original data real or complex voltages? */
    int cfact; /**< "complex factor"; either 1 (for real data) or 2 (for complex data); determines length of substridesize (twice the stridesize for real data [i.e. need 2N samples for nchan] and equal to the stridesize for complex data [need N samples for nchan]) */
} FxKernel;

static cf32 lut2bit[256][4]; /**< Look up table for two bit data, will contain the 4 time samples corresponding to a given byte */

void initLUT2bitReal();

void init_FxKernel(FxKernel *kernel, int nant, int nchan, int nfft, int nbit, double lo, double bw, u8 **inputData);

void setDelays(FxKernel *kernel, double ** d, double * f);

void getStationDelay(FxKernel *kernel, int antenna, int fftindex, double *meandelay, double *a, double *b);

void unpack(FxKernel *kernel, u8 * inputdata, cf32 ** unpacked, int offset);

void unpackReal2bit(u8 *inputdata, cf32 **unpacked, int offset, int nsamp);

void fringerotate(FxKernel *kernel, cf32 ** unpacked, f64 a, f64 b);

void dofft(FxKernel *kernel, cf32 ** unpacked, cf32 ** channelised);

void fracSampleCorrect(FxKernel *kernel, cf32 ** channelised, f64 fracdelay);

void conjChannels(FxKernel *kernel, cf32 ** channelised, cf32 ** conjchannels);

void process(FxKernel *kernel);

void saveVisibilities(const char *outfile, FxKernel *kernel);

void freeFxKernel(FxKernel *kernel);

void startTiming(struct timespec *starttime, char *starttimestring);

void saveLog(long long diff_ms, char *starttimestring);

void endTiming(struct timespec *starttime, long long *diff_ms);

#endif //TESTS_DIFX_KERNEL_H
