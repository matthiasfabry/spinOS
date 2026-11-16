#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Solve Kepler's equation for many phases:
 *   E - e * sin(E) = 2*pi*phase
 *
 * Returns 0 on success.
 * Returns non-zero on invalid inputs or solver failure.
 */
int spinos_kepler_ecc_anom_f32(float e, const float *phase, int n, float *out) {
    if (phase == NULL || out == NULL || n < 0) {
        return 1;  // Invalid input
    }
    if (!(e >= 0.0f && e < 1.0f)) {
        return 2;  // invalid eccentricity
    }
    if (n == 0) {
        return 0;  // Nothing to do
    }

    const float two_pi = (float)(2.0 * M_PI);

    // loop over all phases and solve Kepler's equation
    for (int i = 0; i < n; i++) {
        float ph = phase[i];
        if (!isfinite(ph)) {
            return 3;  // Invalid phase
        }

        ph = fmodf(ph, 1.0f);
        if (ph < 0.0f) {
            ph += 1.0f;
        }  // fold phase to [0, 1)

        float M = two_pi * ph; // Mean anomaly

        float E = M + e * sinf(M) * (1.0f + e * cosf(M));  // small e expansion guess
        float lo = 0.0f;
        float hi = two_pi;

        int converged = 0;
        // do Newton-Raphson iterations
        for (int it = 0; it < 32; it++) {
            float s = sinf(E);
            float c = cosf(E);
            float f = E - e * s - M;  // this is keplers equation
            if (fabsf(f) < 1e-6f) {
                converged = 1;
                break;
            }
            float fp = 1.0f - e * c;  // derivative of keplers equation
            if (fabsf(fp) < 1e-6f) {
                break;
            }
            float d = f / fp;
            float next_E = E - d;
            // fold E to [0, 2*pi)
            if (next_E < 0.0f) {
                next_E += two_pi;
            } else if (next_E > two_pi) {
                next_E -= two_pi;
            }
            if (!(next_E >= lo && next_E <= hi) || !isfinite(next_E)) {
                break;
            }
            E = next_E;
        }

        // If Newton-Raphson did not converge, use bisection method
        if (!converged) {
            float flo = lo - e * sinf(lo) - M;
            float fhi = hi - e * sinf(hi) - M;
            if (flo * fhi > 0.0f) {
                return 4;  // No guaranteed root in [0, 2*pi)
            }
            for (int it = 0; it < 64; it++) {
                E = 0.5f * (lo + hi);
                float fmid = E - e * sinf(E) - M;
                if (fabsf(fmid) < 1e-6f) {
                    converged = 1;
                    break;
                }
                if (flo * fmid <= 0.0f) {
                    hi = E;
                    fhi = fmid;
                } else {
                    lo = E;
                    flo = fmid;
                }
            }
            if (!converged) {
                return 5;  // Failed to converge
            }
        }

        out[i] = E;
    }

    return 0;
}
