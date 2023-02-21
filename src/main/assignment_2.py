import numpy as np
def main():

    x = 3.7
    x_arr = [3.6,3.8,3.9]
    y_arr = [1.675, 1.436, 1.318]
    question_1 = nevilles_method(x_arr, y_arr, x)

    x_arr2 = [7.2, 7.4, 7.5, 7.6]
    y_arr2 = [23.5492, 25.3913, 26.8224, 27.4589]
    question_2 = divided_difference_table(x_arr2, y_arr2)
    
    question_3 = get_approximate_result(question_2, x_arr2, 7.3)
    question_4 = hermite_interpolation()
    
    np.set_printoptions(precision=15, suppress=True, linewidth=1000)
    print(question_1)
    print()
    print(np.diagonal(question_2)[1:4])
    print()
    print(question_3)
    print()
    np.set_printoptions(precision=7, suppress=True, linewidth=1000)
    print(question_4)
    print()
    cubicInterpolation()

def nevilles_method(x_points, y_points, x):
    matrix = np.zeros((len(x_points), len(y_points)))
    for counter, row in enumerate(matrix):
        row[0] = y_points[counter]
    num_of_points = len(x_points)
    for i in range(1, num_of_points):
        for j in range(1, i + 1):
            first_multiplication = (x - x_points[i - j]) * matrix[i][j-1]
            second_multiplication = (x - x_points[i]) * matrix[i-1][j-1]
            denominator = x_points[i] - x_points[i-j]
            coefficient = (first_multiplication - second_multiplication) / denominator
            matrix[i][j] = coefficient
    
    return matrix[num_of_points - 1][num_of_points - 1]


def divided_difference_table(x_points, y_points):
    
    size: int = len(x_points)
    matrix: np.array = np.zeros((size, size))
    
    for index, row in enumerate(matrix):
        row[0] = y_points[index]
    
    for i in range(1, size):
        for j in range(1, size):
            numerator = matrix[i][j - 1] - matrix[i - 1][j - 1]
            denominator = x_points[i] - x_points[i - j]
            operation = numerator / denominator
            #matrix[i][j] = '{0:.7g}'.format(operation)
            matrix[i][j] = operation
    
    return matrix

def get_approximate_result(matrix, x_points, value):
    # p0 is always y0 and we use a reoccuring x to avoid having to recalculate x 
    reoccuring_x_span = 1
    reoccuring_px_result = matrix[0][0]
    
    # we only need the diagonals...and that starts at the first row...
    for index in range(1, len(x_points)):
        polynomial_coefficient = matrix[index][index]
        # we use the previous index for x_points....
        reoccuring_x_span *= (value - x_points[index - 1])
        
        # get a_of_x * the x_span
        mult_operation = polynomial_coefficient * reoccuring_x_span
        # add the reoccuring px result
        reoccuring_px_result += mult_operation
    
    # final result
    return reoccuring_px_result


def apply_div_dif(matrix: np.array):
    size = len(matrix)
    for i in range(2, size):
        for j in range(2, i+2):
            # skip if value is prefilled (we dont want to accidentally recalculate...)
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue
            
            # get left cell entry
            left: float = matrix[i][j - 1]
            # get diagonal left entry
            diagonal_left: float = matrix[i - 1][j - 1]
            # order of numerator is SPECIFIC.
            numerator: float = (left - diagonal_left)
            # denominator is current i's x_val minus the starting i's x_val....
            denominator = matrix[i][0] - matrix[i - j + 1][0]
            # something save into matrix
            operation = numerator / denominator
            matrix[i][j] = operation
    
    return matrix
def hermite_interpolation():
    x_points = [3.6, 3.8, 3.9]
    y_points = [1.675, 1.436, 1.318]
    slopes = [-1.195, -1.188, -1.182]
    # matrix size changes because of "doubling" up info for hermite 
    num_of_points = len(x_points)
    matrix = np.zeros((2 * num_of_points, 2 * num_of_points))
    # populate x values (make sure to fill every TWO rows)
    for i, x in enumerate(x_points):
        matrix[2 * i][0] = x
        matrix[2 * i + 1][0] = x

    # pre-populate y values (make sure to fill every TWO rows)
    for i, x in enumerate(y_points):
        matrix[2 * i][1] = x
        matrix[2 * i + 1][1] = x

    # pre-populate with derivatives (make sure to fill every TWO rows. starting row changes)
    for i, x in enumerate(slopes):
        matrix[2 * i + 1][2] = x
    return apply_div_dif(matrix)

def cubicInterpolation():
    x = [2, 5, 8, 10]
    y = [3, 5, 7, 9]
    
    size = len(x)
    arr1 = np.zeros(size - 1)
    arr2 = np.zeros(size)

    for i in range(size - 1):
        arr1[i] = x[i + 1] - x[i]
        arr2[i + 1] = (y[i+1] - y[i]) / arr1[i]

    arr3 = np.zeros((size, size))
    arr4 = np.zeros(size)

    for row in range(1, size - 1):
        arr3[row][row - 1] = arr1[row - 1]
        arr3[row][row] = 2 * (arr1[row - 1] + arr1[row])
        arr3[row][row + 1] = arr1[row]
        arr4[row] = 3 * (arr2[row + 1] - arr2[row])

    arr3[0][0] = 1
    arr3[size - 1][size - 1] = 1

    x = np.linalg.solve(arr3, arr4)
    print(arr3)
    print()
    print(arr4)
    print()
    print(x)
if __name__ == "__main__":
    main()