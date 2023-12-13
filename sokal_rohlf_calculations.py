# statistical calculations based on "Biometry" by Robert R. SOKAL and F. James ROHLF, 3rd edition
import numpy as np
import scipy.stats as stats


def regress_confidence_sokal_rohlf(X, Y, alpha):

    q = 1-alpha
    # calculate averages
    X_bar = np.mean(X)
    Y_bar = np.mean(Y)

    #difference of datapoints from mean
    x_small = X - X_bar
    y_small = Y - Y_bar

    # sum of the square differences in x
    sum_x2 = np.nansum(np.square(x_small))

    # multiple differences from mean to get xy and sum
    xy = x_small*y_small
    sum_xy = np.nansum(xy)

    # regression coefficient (slope of the regression)
    # Check if sum_x2 is zero or NaN
    if sum_x2 == 0 or np.isnan(sum_x2):
        # Handle the case where division is not possible or meaningful
        b_yx = np.nan  # or set to some other appropriate value
        a = np.nan
        r2 = np.nan
        CI_alpha_slope = np.nan
        t_stat = np.nan
        sb = np.nan
        y_err = np.nan
    else:
        # Calculate b_yx
        b_yx = sum_xy / sum_x2

        # y intercept 
        a = Y_bar - b_yx*X_bar

        # predicted Y based on regression
        Y_hat = b_yx*X + a

        # unexplained sum of squares
        d2_yx = np.sum(np.square(Y-Y_hat))

        # explained variance:
        s2_y_hat = np.sum(np.square(Y_hat - Y_bar))/(len(Y)-1)

        # unexplained variance:
        dof = (len(Y) - 2)

        if dof==0:
            s2_yx = np.nan
        else:
            s2_yx = d2_yx/dof

        # total variance:
        s2_Y = np.sum(np.square(Y - Y_bar))/(len(Y)-1)

        # coefficient of determination (r2)
        r2 = s2_y_hat/s2_Y

        # standard error of the regression coefficient
        sb = np.sqrt(s2_yx/sum_x2);

        # regression coefficient (slope) different than zero?
        # t_s = (b_yx - 0)/sb

        # t stat
        t_stat = stats.t.ppf(1-q/2, dof)

        CI_alpha_slope = [b_yx - t_stat*sb, b_yx + t_stat*sb]

        # standard error of Y_hat for a given value of X (every value)
        sy_hat = (s2_yx*(1/len(Y) + ((X-X_bar)**2)/sum_x2))**0.5

        # upper confidence interval
        y_err = t_stat*sy_hat

    return b_yx, a, r2, CI_alpha_slope, t_stat*sb/2, y_err