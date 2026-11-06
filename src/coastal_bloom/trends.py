"""Exploratory seasonal linear trend with calendar-aware Newey–West errors."""
import numpy as np


def monthly_trend(dates, values, hac_lags=12):
    """Fit intercept + linear years + eleven calendar-month indicators.

    Returns slope/year and asymptotic HAC standard error, not a causal test.
    No interpolation of missing months. Bartlett lag weights use actual month
    separation. Minimum 36 finite months; inference remains exploratory for
    short records. Full analysis also needs block-bootstrap sensitivity.
    """
    dates = np.asarray(dates, dtype="datetime64[M]")
    y = np.asarray(values, float)
    if dates.ndim != 1 or dates.shape != y.shape or np.any(np.isnat(dates)):
        raise ValueError("provide aligned 1D monthly dates and values")
    if not isinstance(hac_lags, (int, np.integer)) or hac_lags < 0:
        raise ValueError("hac_lags must be a nonnegative integer")
    months = dates.astype(int)
    if len(np.unique(months)) != len(months):
        raise ValueError("duplicate months are not allowed")
    keep = np.isfinite(y)
    months, y = months[keep], y[keep]
    order = np.argsort(months)
    months, y = months[order], y[order]
    if len(y) < 36 or months[-1]-months[0] < 35:
        raise ValueError("at least 36 observed months spanning 3 years required")
    t = (months-months[0])/12
    seasonal = (months[:, None] % 12 == np.arange(1, 12)).astype(float)
    X = np.column_stack([np.ones(len(y)), t, seasonal])
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError("seasonal trend design is rank deficient")
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    residual = y-X@beta
    scores = X*residual[:, None]
    meat = scores.T@scores
    positions = {int(month): i for i, month in enumerate(months)}
    for lag in range(1, hac_lags+1):
        for i, month in enumerate(months):
            j = positions.get(int(month-lag))
            if j is not None:
                cross = np.outer(scores[i], scores[j])
                meat += (1-lag/(hac_lags+1))*(cross+cross.T)
    bread = np.linalg.inv(X.T@X)
    cov = bread@meat@bread * len(y)/(len(y)-X.shape[1])
    return {"slope_per_year": float(beta[1]),
            "hac_standard_error": float(np.sqrt(max(cov[1, 1], 0))),
            "n_months": len(y), "hac_lags": hac_lags}
