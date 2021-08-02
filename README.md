# process-observers

MATLAB scripts for simulating process observers for online state estimation and sensor fusion.

## Contents

Observers currently included:
- [luenberger_filter.m](luenberger_filter.m) - Luenberger observer (steady-state)
- [kalman_filter_ss.m](kalman_filter_ss.m) - steady-state Kalman filter
- [kalman_filter.m](kalman_filter.m) - dynamic Kalman filter

Multi-model observers:
- [mkf_filter.m](mkf_filter.m) - general purpose multi-model Kalman filter
- [mkf_filter_RODD.m](mkf_filter_RODD.m) - multi-model Kalman filter for state estimation in the presence of randomly-occurring deterministic disturbances (RODDs) as described in Robertson et al. (1995).

## Tutorial

See the following LiveScript for an introduction on how to use these functions:

- [RODD_code_tutorial.mlx](RODD_code_tutorial.mlx)

## References

1. Kalman, R. E. (1960). A New Approach to Linear Filtering and Prediction Problems. Journal of Basic Engineering. 82: 35–45. https://doi.org/10.1115%2F1.3662552.
2. Robertson, D. G., Kesavan, P., & Lee, J. H. (1995). Detection and estimation of randomly occurring deterministic disturbances. Proceedings of 1995 American Control Conference - ACC?95, 6, 4453?4457. https://doi.org/10.1109/ACC.1995.532779
