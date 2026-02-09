#ifndef HALALSEQ_CALIBRATE_H
#define HALALSEQ_CALIBRATE_H

/* Spike-in calibration: estimate prior parameters for d and b from
 * calibration samples with known composition.
 */

typedef struct {
    double *true_w;       /* [n_species] known weight fractions */
    double *obs_reads;    /* [n_species * n_markers] observed read counts */
    int n_species;
    int n_markers;
} calibration_sample_t;

typedef struct {
    double d_mu, d_sigma;    /* LogNormal prior on DNA yield */
    double b_mu, b_sigma;    /* LogNormal prior on PCR bias */
    int n_samples;
} calibration_result_t;

/* Estimate bias priors from spike-in calibration data */
calibration_result_t *calibrate_estimate(const calibration_sample_t *samples,
                                          int n_samples);

void calibrate_result_destroy(calibration_result_t *r);

#endif /* HALALSEQ_CALIBRATE_H */
