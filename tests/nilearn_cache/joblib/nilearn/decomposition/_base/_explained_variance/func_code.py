# first line: 582
def _explained_variance(X, components, per_component=True):
    """Score function based on explained variance.

    Parameters
    ----------
    X : ndarray
        Holds single subject data to be tested against components.

    components : array-like
        Represents the components estimated by the decomposition algorithm.

    per_component : bool, default=True
        Specify whether the explained variance ratio is desired for each
        map or for the global set of components_.

    Returns
    -------
    score : ndarray
        Holds the score for each subjects. score is two dimensional if
        per_component = True.

    """
    full_var = np.var(X)
    n_components = components.shape[0]
    S = np.sqrt(np.sum(components**2, axis=1))
    S[S == 0] = 1
    components = components / S[:, np.newaxis]
    projected_data = components.dot(X.T)
    if per_component:
        res_var = np.zeros(n_components)
        for i in range(n_components):
            res = X - np.outer(projected_data[i], components[i])
            res_var[i] = np.var(res)
            # Free some memory
            del res
        return np.maximum(0.0, 1.0 - res_var / full_var)
    else:
        lr = LinearRegression(fit_intercept=True)
        lr.fit(components.T, X.T)
        res = X - lr.coef_.dot(components)
        res_var = row_sum_of_squares(res).sum()
        return np.maximum(0.0, 1.0 - res_var / row_sum_of_squares(X).sum())
