# Developed by:
# https://github.com/rebberchicken
# DO NOT REPOST OR USE FOR YOUR ACADEMICS.

# This program implements the Polynomial Regression from Exercise 6.
# Gauss-Jordan Elimination (for extracting coefficients)
GaussJordanMethod <- function(augmentedMatrix) {
  n <- nrow(augmentedMatrix)
  solution <- numeric(n)
  
  # Forward Elimination
  for (i in 1:n) {
    if (i != n) {
      # Finds the pivot row with the maximum absolute value in the current column
      pivot_row <- which.max(abs(augmentedMatrix[i:n, i])) + i - 1
      
      # Checks if the pivot element is zero, if so, it should STOP since no unique solution exists
      if (augmentedMatrix[pivot_row, i] == 0) {
        return(list(solution = NA, variables = colnames(augmentedMatrix)[1:n], augcoeffmatrix = augmentedMatrix)) # returns solution as NA
      }
      
      # Swap rows for partial pivoting
      augmentedMatrix[c(i, pivot_row),] <- augmentedMatrix[c(pivot_row, i),]
    }
    
    # Updates the pivot row
    augmentedMatrix[i,] <- augmentedMatrix[i,] / augmentedMatrix[i, i]
    
    # Eliminate other rows
    for (j in 1:n) {
      if (i == j) { next } # Continues if i is equal to j
      normalized_row <- augmentedMatrix[j, i] * augmentedMatrix[i,] # Find normalized row
      augmentedMatrix[j,] <- augmentedMatrix[j,] - normalized_row  # Subtract the normalized row to eliminate the element below the pivot
    }
  }
  
  return(list(variables = colnames(augmentedMatrix)[1:n], augcoeffmatrix = augmentedMatrix, solution = as.vector(augmentedMatrix[,n + 1])))
  # Converted the augmented matrix to a vector using the as.vector function
}

# Function for polynomial regression which accepts an integer and list (data) as inputs
PolynomialRegression <- function(order, data, x_val) {
  # Checks if the order of the polynomial is less than 1
  if (order < 1) {
    return(list(augcoeffmatrix = NA, coefficients = NA, polynomial_string = NA, polynomial_function = NA, estimate = NA)) # returns NA
  }
  
  # Check if data is not a list or if the list (data) has two vectors
  if (!is.list(data) || length(data) != 2) {
    return(list(augcoeffmatrix = NA, coefficients = NA, polynomial_string = NA, polynomial_function = NA, estimate = NA)) # returns NA
  }
  
  x <- data[[1]] # Plots x as the first vector
  y <- data[[2]] # Plots y as the second vector
  
  # Checks if the vectors are of the same length
  if (length(x) != length(y)) {
    return(list(augcoeffmatrix = NA, coefficients = NA, polynomial_string = NA, polynomial_function = NA, estimate = NA)) # returns NA
  }
  
  n <- length(x) # Initializes nth order polynomial based on the length of x, assuming that they are of the same length as y
  
  # Checks if there are enough data points for the given order
  if (n <= order) {
    return(list(augcoeffmatrix = NA, coefficients = NA, polynomial_string = NA, polynomial_function = NA, estimate = NA)) # returns NA
  }
  
  # Initialization for the augmented matrix
  augmentedMatrix <- matrix(0, nrow = order + 1, ncol = order + 2)
  
  # Appends the sum of x with respect to the vector to the matrix
  for (i in 1:(order + 1)) { # Iterates over each row
    for (j in 1:(order + 1)) { # Iterates over each column
      augmentedMatrix[i, j] <- sum(x^(i + j - 2)) # Based on y = a sub 0 + a sub 1 * x + a sub 2 * x^2 + a sub order * x^order
      # X sub i*j = summation of x^i+j-2 sub k where k = 1 // Based on the number of columns which is the order of coeff + 2 
    }
    augmentedMatrix[i, order + 2] <- sum(x^(i - 1) * y) # Updates the RHS, summation of x sub i raised to m multiplied by y sub i
  }
  
  result <- GaussJordanMethod(augmentedMatrix) # Perform Gauss Jordan Elimination to get the coefficients
  coefficients <- result$solution[1:(order + 1)] # Extracts the coefficients from the result
  
  
  # Create the polynomial strings
  polynomial_string <- paste("function(x) ", paste(coefficients[1], "+", paste(paste(coefficients[-1], "* x ^", seq(1, order), collapse = " + "), collapse = " + ")))
  # String format: function(x) a_0 + a_1 * x ^ 1 + a_2 * x ^ 2 + ...
  # paste() is responsible for appending the coefficients to another string
  # While the collapse variable allows the paste function to be part of the vector with "+" as a separator
  
  # Evaluate the polynomial function at X value
  estimate <- sum(coefficients * x_val^(0:(length(coefficients) - 1)))
  
  # Returns the labeled list
  return(list(augcoeffmatrix = augmentedMatrix, 
        coefficients = coefficients,
        polynomial_string = polynomial_string,
        polynomial_function = eval(parse(text = polynomial_string)),
        estimate = estimate))
  # eval() returns the expression, parse() converts the strings to an R expression
}
