# TANDANG, Bernard Jezua R.
# 2021-09992
# CMSC 150 - B1L

# This program implements the Quadratic Spline Interpolation.
source("GaussianElimination.R") # From the previous exercise
poly.qsi <- function(data, num){ # Function that evaluates X using QSI, accepts a list and a value num to be evaluated
  
  var <- c("a","b","c") # Consistent unknown variables for the equations
  unknown <- length(var)
  unknownVar <- c() # Vector to store generated unknown variables
  x <- data[[1]]
  y <- data[[2]]
  n <- length(x) - 1 # n = Number of data points - 1
  
  # Matrix for QSI
  qsiMat <- matrix(0, nrow = 3*n, ncol = 3*n+1) 
  counter = 1
  
  # Loop to generate unknowns for interval * n unknowns equations 
  for(i in 1:n){
    for(j in 1:unknown){
      unknownVar[counter] <- paste(var[j], i, sep = "") # Creates ("a1", "b1", "c1")
      counter = counter + 1
    }
  }

  cnames <- c(unknownVar,"RHS")
  rnames <- 1:length(unknownVar)
  
  rownames(qsiMat) <- rnames
  colnames(qsiMat) <- cnames

  # Vectors to store the generated equations
  eq1 <- c()
  eq2 <- c()
  eq3 <- c()
  eq4 <- c()
  
  # Loop conditions for the equations
  # Condition 1: Internal knots
  for(i in 2:n){
    eq1[i*2-3] <- paste(x[i]*x[i], "*a",i-1," + ", x[i], "*b", i-1, " + ", 1, "*c", i-1, " + ", -1*y[i], sep = "") 
    eq1[i*2-2] <- paste(x[i]*x[i], "*a",i," + ", x[i], "*b", i, " + ", 1, "*c", i, " + ", -1*y[i], sep = "") 
  }
  
  # e.g. in test case:
  # 16.5*16.5 *a1 + 16.5 *b1 + c1 = 6.1 -> 16.5*16.5 *a1 + 16.5 *b1 + c1 + -6.1 = 0
  
  # Condition 3: Internal knots again
  for(i in 2:n){
    eq4[i-1] <- paste(2*x[i], "*a",i-1," + ", 1, "*b", i-1, " + ", -1*2*x[i], "*a", i, " + ", -1*1, "*b",i, sep = "")
  }
  
  # e.g. in test case:
  # 2*16.5 *a1 + 1*b1 = 2*16.5 *a2 + 1*b2 -> 2*16.5 *a1 + 16.5 *b1 + -1*(2*16.5) *a2 + -1*(1*b2)
  
  # Condition 2: External knots
  eq2 <- paste(x[1]*x[1], "*a",1," + ", x[1], "*b", 1, " + ", 1, "*c",1, " + ", -1*y[1], sep = "") 
  eq3 <- paste(x[n+1]*x[n+1], "*a",n," + ", x[n+1], "*b", n, " + ", 1, "*c",n, " + ", -1*y[n+1], sep = "")  
  
  # e.g. in test case:
  # 9*9 *a1 + 9 *b1 + c1 = 4 -> 9*9 *a1 + 9 *b1 + c1 + -4 = 0
  
  # Combine all equations into 1 vector
  eq5 <- c(eq1, eq2, eq3, eq4)
  count <- 2 # Counter initialization for the row of the QSI
  # Loop to populate the QSI matrix using the equations similar to third exercise
  for(j in 1:length(eq5)+1){
    equation <- strsplit(eq5[j-1], " \\+ ")[[1]] # Splits j-1th equation of the vector
    for(i in equation){ 
      if(grepl("*", i, fixed=TRUE)){ # Condition to check if i-th element is the RHS
        coeff <- strsplit(i, "\\*")[[1]] # Splits equation "3*a1" -> "3" "a1"
        if(coeff[2] == "a1"){ # Condition 4 where a1 = 0
          qsiMat[count, coeff[2]] <- 0
        } else { # Populates the coefficient side of the QSI Matrix
          qsiMat[count,coeff[2]] <- as.numeric(coeff[1]) 
        }
      } else{ # Populates RHS
        qsiMat[count,"RHS"] <- as.numeric(i)*-1
      }
    }
    count <- count + 1 # Increments counter for QSI row
  }
  
  # Remove the first row and column of QSI matrix since a1=0
  qsiMat <- qsiMat[-1,-1] # Gets rid of the first row and column of QSI matrix since a1=0
  
  # Use Gaussian Elimination to evaluate the QSI matrix and get the solution set
  qsi <- GaussianMethod(qsiMat)$solution
  #print(GaussianMethod(qsiMat)$augcoeffmatrix)
  
  # Adds Condition 4: "a1 = 0" to solution vector
  qsi <- c(0, qsi) 
  
  # Generates the Quadratic Equation for every interval
  polyMat <- matrix(qsi, n, unknown, byrow = TRUE) # Matrix to store the solution set in to 3 columns for a,b,c
  polynomial <- c() # Vector to store final equations
  
  for(i in 1:n){ # For loop to create and evaluate quadratic equations
    polynomial[i] <- paste("function (x) " , polyMat[i,1], "*x^2 + ", polyMat[i,2], "*x", " + ", polyMat[i,3], sep="" )
    polynomial[i] <- list(eval(parse(text = polynomial[i])))
    # paste() is responsible for appending the coefficients to another string
  }
   
  # For loop to evaluate num with the appropriate interval equation
  for(i in 1:n){
    if(x[i] <= num && num <= x[i+1]){
      ret <- polynomial[[i]]
      ret <- ret(num) # Estimate
      break
    } else { # Returns NA if the num is not within the intervals
      ret <- NA
    }
    
  }
  return(list(qsi.fxns = polynomial, y = ret)) # Returns a list of string for QSI functions and estimate
}
