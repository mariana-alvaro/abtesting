from scipy import stats
from scipy.stats import t as t_dist
from scipy.stats import chi2
import math
from abtesting_test import *

# You can comment out these lines! They are just here to help follow along to the tutorial.
# print(t_dist.cdf(-2, 20)) # should print .02963
# print(t_dist.cdf(2, 20)) # positive t-score (bad), should print .97036 (= 1 - .2963)
#
# print(chi2.cdf(23.6, 12)) # prints 0.976
# print(1 - chi2.cdf(23.6, 12)) # prints 1 - 0.976 = 0.023 (yay!)


# TODO: Fill in the following functions! Be sure to delete "pass" when you want to use/run a function!
# NOTE: You should not be using any outside libraries or functions other than the simple operators (+, **, etc)
# and the specifically mentioned functions (i.e. round, cdf functions...)

def slice_2D(list_2D, start_row, end_row, start_col, end_col):
    '''
    Splices a the 2D list via start_row:end_row and start_col:end_col
    :param list: list of list of numbers
    :param nums: start_row, end_row, start_col, end_col
    :return: the spliced 2D list (ending indices are exclsive)
    '''
    to_append = []
    for l in range(start_row, end_row):
        to_append.append(list_2D[l][start_col:end_col])

    return to_append

def get_avg(nums):
    '''
    Helper function for calculating the average of a sample.
    :param nums: list of numbers
    :return: average of list
    '''
    return sum(nums) / len(nums)


def get_stdev(nums):
    '''
    Helper function for calculating the standard deviation of a sample.
    :param nums: list of numbers
    :return: standard deviation of list
    '''
    avg = get_avg(nums)
    var = 0
    for i in nums:
        var += (i - avg)**2/ (len(nums) -1)
    return var**0.5

def get_standard_error(a, b):
    '''
    Helper function for calculating the standard error, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: standard error of a and b (see studio 6 guide for this equation!)
    '''
    term1 = get_stdev(a)**2 / len(a)
    term2 = get_stdev(b)**2 / len(b)
    se = (term1 + term2)**0.5
    return se

def get_2_sample_df(a, b):
    '''
    Calculates the combined degrees of freedom between two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: integer representing the degrees of freedom between a and b (see studio 6 guide for this equation!)
    HINT: you can use Math.round() to help you round!
    '''
    den1 = (get_stdev(a)**2) / len(a)
    den1 = den1**2
    den1 = den1 / (len(a)-1)

    den2 = (get_stdev(b)**2) / len(b)
    den2 = den2**2
    den2 = den2 / (len(b)-1)

    df = get_standard_error(a,b)**4
    df = df / (den1+den2)
    return round(df)

def get_t_score(a, b):
    '''
    Calculates the t-score, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: number representing the t-score given lists a and b (see studio 6 guide for this equation!)
    '''
    t = (get_avg(a)-get_avg(b))/ get_standard_error(a,b)
    if t > 0:
        return -1*t
    return t

def perform_2_sample_t_test(a, b):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates a p-value by performing a 2-sample t-test, given two lists of numbers.
    :param a: list of numbers
    :param b: list of numbers
    :return: calculated p-value
    HINT: the t_dist.cdf() function might come in handy!
    '''
    t_vals = get_t_score(a, b)
    df = get_2_sample_df(a, b)
    return t_dist.cdf(t_vals,df)


# [OPTIONAL] Some helper functions that might be helpful in get_expected_grid().
def row_sum(observed_grid, ele_row):
    return sum(observed_grid[ele_row])

def col_sum(observed_grid, ele_col):
    sum = 0
    for i in range(len(observed_grid)):
        sum += observed_grid[i][ele_col]
    return sum

def total_sum(observed_grid):
    total = 0
    cols = len(observed_grid[0])
    for i in range(len(observed_grid)):
        for j in range(cols):
            total += observed_grid[i][j]
    return total

def calculate_expected(row_sum, col_sum, tot_sum):
    expected = (row_sum * col_sum) / tot_sum
    return expected

def get_expected_grid(observed_grid):
    '''
    Calculates the expected counts, given the observed counts.
    ** DO NOT modify the parameter, observed_grid. **
    :param observed_grid: 2D list of observed counts
    :return: 2D list of expected counts
    HINT: To clean up this calculation, consider filling in the optional helper functions below!
    '''
    #
    cols = len(observed_grid[0])
    # expected = [[0 for i in range(cols)] for j in range(len(observed_grid))]
    expected_grid = []
    total = total_sum(observed_grid)
    for i in range(len(observed_grid)):
        expected_row = []
        for j in range(cols):
            sumr = row_sum(observed_grid, i)
            sumc = col_sum(observed_grid, j)
            expected_value = calculate_expected(sumr, sumc, total)
            expected_row.append(expected_value)
        expected_grid.append(expected_row)
    return expected_grid

def df_chi2(observed_grid):
    '''
    Calculates the degrees of freedom of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: degrees of freedom of expected counts (see studio 6 guide for this equation!)
    '''
    cols = len(observed_grid[0])
    rows = len(observed_grid)
    df = (rows -1)*(cols -1)
    return df

def chi2_value(observed_grid):
    '''
    Calculates the chi^2 value of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: associated chi^2 value of expected counts (see studio 6 guide for this equation!)
    '''
    expected = get_expected_grid(observed_grid)
    cols = len(observed_grid[0])
    chi = [[0 for i in range(cols)] for j in range(len(observed_grid))]
    for i in range(len(observed_grid)):
        for j in range(cols):
            chi[i][j] = ((expected[i][j] - observed_grid[i][j])**2) / expected[i][j]
    return total_sum(chi)

def perform_chi2_homogeneity_test(observed_grid):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates the p-value by performing a chi^2 test, given a list of observed counts
    :param observed_grid: 2D list of observed counts
    :return: calculated p-value
    HINT: the chi2.cdf() function might come in handy!
    '''
    chi2_val = chi2_value(observed_grid)
    df = df_chi2(observed_grid)
    return (1- chi2.cdf(chi2_val,df))

# These commented out lines are for testing your main functions.
# Please uncomment them when finished with your implementation and confirm you get the same values :)
def data_to_num_list(s):
  '''
    Takes a copy and pasted row/col from a spreadsheet and produces a usable list of nums.
    This will be useful when you need to run your tests on your cleaned log data!
    :param str: string holding data
    :return: the spliced list of numbers
    '''
  return list(map(float, s.split()))

print()
print("TESTING")
print()
# t_test 1:
a_t1_list = data_to_num_list(a1)
b_t1_list = data_to_num_list(b1)
print(get_t_score(a_t1_list, b_t1_list)) # this should be -129.500
print(perform_2_sample_t_test(a_t1_list, b_t1_list)) # this should be 0.0000
# why do you think this is? Take a peek at a1 and b1 in abtesting_test.py :)

# t_test 2:
a_t2_list = data_to_num_list(a2)
b_t2_list = data_to_num_list(b2)
print(get_t_score(a_t2_list, b_t2_list)) # this should be -1.48834
print(perform_2_sample_t_test(a_t2_list, b_t2_list)) # this should be .082379

# t_test 3:
a_t3_list = data_to_num_list(a3)
b_t3_list = data_to_num_list(b3)
print(get_t_score(a_t3_list, b_t3_list)) # this should be -2.88969
print(perform_2_sample_t_test(a_t3_list, b_t3_list)) # this should be .005091

# chi2_test 1:
a_c1_list = data_to_num_list(a_count_1)
b_c1_list = data_to_num_list(b_count_1)
c1_observed_grid = [a_c1_list, b_c1_list]
print(chi2_value(c1_observed_grid)) # this should be 4.103536
print(perform_chi2_homogeneity_test(c1_observed_grid)) # this should be .0427939

# chi2_test 2:
a_c2_list = data_to_num_list(a_count_2)
b_c2_list = data_to_num_list(b_count_2)
c2_observed_grid = [a_c2_list, b_c2_list]
print(chi2_value(c2_observed_grid)) # this should be 33.86444
print(perform_chi2_homogeneity_test(c2_observed_grid)) # this should be 0.0000
# Again, why do you think this is? Take a peek at a_count_2 and b_count_2 in abtesting_test.py :)

# chi2_test 3:
a_c3_list = data_to_num_list(a_count_3)
b_c3_list = data_to_num_list(b_count_3)
c3_observed_grid = [a_c3_list, b_c3_list]
print(chi2_value(c3_observed_grid)) # this should be .3119402
print(perform_chi2_homogeneity_test(c3_observed_grid)) # this should be .57649202
