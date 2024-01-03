# Developed by:
# https://github.com/rebberchicken
# DO NOT REPOST OR USE FOR YOUR ACADEMICS.

# This program implements the Gaussian Elimination from the previous exercise.
GaussianMethod <- function(augmentedMatrix) {
  n <- nrow(augmentedMatrix) 
  solution <- numeric(n) 
  
  # Forward Elimination
  for (i in 1:(n - 1)) {
    # Finds the pivot row with the maximum absolute value in the current column
    pivot_row <- which.max(abs(augmentedMatrix[i:n, i])) + i - 1
    
    # Checks if the pivot element is zero, if so, it should STOP since no unique solution exists
    if (augmentedMatrix[pivot_row, i] == 0) {
      return(list(solution = NA, variables = colnames(augmentedMatrix)[1:n], augcoeffmatrix = augmentedMatrix)) # returns solution as NA
    }
    
    # Swap rows for partial pivoting
    augmentedMatrix[c(i, pivot_row),] <- augmentedMatrix[c(pivot_row, i),]
    
    # Eliminate other rows
    for (j in (i + 1):n) {
      pivot_element <- augmentedMatrix[i, i] # get the pivot element
      multiplier <- augmentedMatrix[j, i] / pivot_element # calculate the multiplier for row elimination
      normalized_row <- multiplier * augmentedMatrix[i, ] # find normalized row
      augmentedMatrix[j,] <- augmentedMatrix[j,] - normalized_row # subtract the normalized row to eliminate the element below the pivot
    }
  }
  
  # Backward Substitution
  solution[n] <- augmentedMatrix[n, n + 1] / augmentedMatrix[n, n] # initializes last variable for backward substitution
  for (i in (n - 1):1) {
    solution[i] <- (augmentedMatrix[i, n + 1] - sum(augmentedMatrix[i, (i + 1):n] * solution[(i + 1):n])) / augmentedMatrix[i, i] # appends solution for each variable
  }
  
  return(list(variables = colnames(augmentedMatrix)[1:n], augcoeffmatrix = augmentedMatrix, solution = solution)) # returns list
}
