import pandas as pd
from scipy import stats
import numpy as np

# from: https://towardsdatascience.com/anomaly-detection-with-generalized-extreme-studentized-deviate-in-python-f350075900e2 
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm 
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm

def test_stat(y, iteration, verboseTF):
    std_dev = np.std(y)
    avg_y = np.mean(y)
    abs_val_minus_avg = abs(y - avg_y)
    max_of_deviations = max(abs_val_minus_avg)
    max_ind = np.argmax(abs_val_minus_avg)
    cal = max_of_deviations/ std_dev
    if verboseTF:
        print('Test {}'.format(iteration))
        print("Test Statistics Value(R{}) : {}".format(iteration,cal))
    return cal, max_ind

def calculate_critical_value(size, alpha, iteration, verboseTF):
    t_dist = stats.t.ppf(1 - alpha / (2 * size), size - 2)
    numerator = (size - 1) * np.sqrt(np.square(t_dist))
    denominator = np.sqrt(size) * np.sqrt(size - 2 + np.square(t_dist))
    critical_value = numerator / denominator
    if verboseTF:
        print("Critical Value(λ{}): {}".format(iteration, critical_value))
    return critical_value

def check_values(R, C, inp, max_index, iteration):
    if R > C:
        print('{} is an outlier. R{} > λ{}: {:.4f} > {:.4f} \n'.format(inp[max_index],iteration, iteration, R, C))
    else:
        print('{} is not an outlier. R{}> λ{}: {:.4f} > {:.4f} \n'.format(inp[max_index],iteration, iteration, R, C))

def ESD_Test(input_series, alpha, max_outliers, verboseTF, NoOutputTF):
    stats = []
    critical_vals = []
    tested_values = []
    max_i = 0
    for iterations in range(1, max_outliers + 1):
        stat, max_index = test_stat(input_series, iterations, verboseTF)
        critical = calculate_critical_value(len(input_series), alpha, iterations, verboseTF)
        if verboseTF:
            check_values(stat, critical, input_series, max_index, iterations)
        tested_values.append(input_series[max_index])
        input_series = np.delete(input_series, max_index)
        critical_vals.append(critical)
        stats.append(stat)
        if stat > critical:
            max_i = iterations
    if ~NoOutputTF:
        print('H0:  there are no outliers in the data')
        print('Ha:  there are up to 10 outliers in the data')
        print('')
        print('Significance level:  α = {}'.format(alpha))
        print('Critical region:  Reject H0 if Ri > critical value')
        print('Ri: Test statistic')
        print('λi: Critical Value')
        print(' ')
    df = pd.DataFrame({'i' :range(1, max_outliers + 1), 'Ri': stats, 'λi': critical_vals , 'Vals': tested_values})
    
    def highlight_max(x):
        if x.i == max_i:
            return ['background-color: yellow']*4
        else:
            return ['background-color: white']*4
    df.index = df.index + 1
    print('Number of outliers {}'.format(max_i))
    print(max_i)
    return  df.style.apply(highlight_max, axis = 1), max_i, tested_values