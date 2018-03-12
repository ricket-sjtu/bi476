# Lab 1: Introductory Math For Biostatistics Using R
___


## 1. Calculus in R

### 1.1 Differentials
The R function `stats::deriv()`, `Deriv::Deriv()` can be used to infer the 
1-order derivative and second derivatives (Hessian). Here is an example:
```r
f <- function(x) sin(3*x)
g <- deriv(as.expression(body(f)), names(formals(f)),
        function.arg=TRUE, hessian=FALSE)
g
function (x) 
{
    .expr1 <- x^2
    .expr2 <- sin(x)
    .expr6 <- 2 * x
    .value <- .expr1 * .expr2 + log(.expr1)
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
    .grad[, "x"] <- .expr6 * .expr2 + .expr1 * cos(x) + .expr6/.expr1
    attr(.value, "gradient") <- .grad
    .value
}
g(5)
[1] -20.75423
attr(,"gradient")
             x
[1,] -2.097688
```

Here we have a function
$$
f(x, y) = x\sin y + \log x^2
$$
- Compute the first-order derivative w.r.t. $x, y$
- Compute the second-order derivative w.r.t. $x, y$
- Plot the function in the range of (0 ,10)
- Where will the function get the optimum (maximum/minimum) value?

### 1.2 Taylor Expansion

Let's see the Taylor expansion of $f(x) = \log(1+x)$ at $x_0=0$:
$$
f(x) = f(x_0) + \dots
$$

And in the package `pracma`, we have the `taylor()` function for computing 
the Taylor's expansion:
```r
f <- function(x) log(1+x)
p <- taylor(f, x0=0, n=4)
x <- seq(-1, 1, length.out=100)
yf <- f(x)
yp <- polyval(p, x)
plot(x, yf, type="l", col="gray", lwd=3)
lines(x, yp, col="blue")
grid()
```


### 1.3 The properties of a function
- `formals`: `formals()`
- `body`: `body`
- `environment`: 

## 2. Linear algebra in R

### 2.1 Matrix multiplication
For the following matrix:
$$
\mathbf{A} = \begin{bmatrix}
1 & 1 & 1 & 1 & 1\\\\
1 & 1 & 1 & 1 & 1\\\\
1 & 1 & 1 & 1 & 1\\\\
1 & 1 & 1 & 1 & 1\\\\
1 & 1 & 1 & 1 & 1
\end{bmatrix}, B = \mathrm{diag}(1,2,3,4,5)
$$
What happens to
- $AB$
- $BA$

### 2.2 Least squares fitting
Consider the table of pairs $(x_i, y_i)$ below:
```
x: 1.00,2.00,3.00,4.00,5.00
y: 3.70,4.20,4.90,5.70,6.00
```
Use least squres fitting to approximate the linear models:
$$
y_i = \beta_0 + \beta_1 x_i, i=1,\dots,5
$$

### 2.3 Singular value decomposition (SVD)
Consider the matrix
$$
\mathbf{A} = \begin{bmatrix}
-2 & 11\\\\
-10 & 5
\end{bmatrix}
$$
- Compute the singular values, left singular vectors, and 
    right singular vectors of $\mathbf{A}$.
- Draw a careful, labeled graph of the unit ball in 
    $\mathbb{R}^2$ and its image under $\mathbf{A}$.
- Plot the singular vectors with the coordinates of their 
    vertices marked within the above picture.
- Find $\mathbf{A}^{-1}$ via **SVD**.
- Compute the eigenvalues $\lambda_1, \lambda_2$ of $\mathbf{A}$.
- Verify that $\mathrm{det} \mathbf{A} = \lambda_1 \lambda_2$
    and $|\mathrm{det}\mathbf{A}| = \sigma_1 \sigma_2$.
- Are $\mathbf{AA^T}$ and $\mathbf{A^TA}$ positive-definite? 
    positive semi-definite? Use the eigenvalues to answer the question.

## 3. Optimization in R

Optimization is a big and essential part of the 
statistical modeling, which can be used to solve many 
modeling problems, e.g., maximum likelihood.

Here we will cover 5 popular optimization techniques.

### 3.1 Golden-section search
Golden section search is a line-search method for global 
optimization in 1-dimensional context. It is a direct-search 
technique as it samples the function to approximate a 
derivative rather than computing it directly.

Here is an example for Golden section search method in R for 
solving a 1d nonlinear unconstrained optimization function:

```r
# define a 1d basin function, optima at f(0)=0
basin <- function(x) x^2 - 3

# locate the minimum of the function using golden section line search
result <- optimize(
    basin, # the function to be minimized
    c(-5, 5), # the bounds of the parameters
    maximum=FALSE, # function minima
    tol=1e-8) # the size of the final bracketing

# display the result
print(result$minimum)
print(result$objective)

# plot the function
x <- seq(-5, 5, length.out=100)
y <- basin(expand.grid(x))
plot(x, y, xlab='x', ylab='f(x)', type='l')

# plot the solution as a point
points(result$minimum, result$objective, 
        col='red', pch=19)
rect(result$minimum-0.3, result$objective-0.7,
    result$minimum+0.3, result$objective+0.7,
    lwd=2) 
```

### 3.2 Nelder-Mead

The Nelder-Mead method is an optimization method for 
**multi-dimensional nonlinear unconstrained function**.

- It is also a direct search method in that it does not 
require a function gradient during the procedure.
- It is a pattern search that it uses a geometric pattern 
to explore the problem space.

Here is an example to illustrate the usage of Nelder-Mead 
algorithm for solving a two-dimensional nonlinear optimization 
problem.

```r
# definition of the 2D Rosenbrock function, optima at c(1,1)
rosenbrock <- function(v) {
    (1-v[1])^2 + 100*(v[2]-v[1]*v[1])^2
}

# find the minimum using the Nelder-Mead method
result <- optim(
    c(runif(1,-3,3), runif(1,-3,3)), # start position
    rosenbrock, # the function to be minized
    NULL, # no function gradient
    method="Nelder-Mead", # use the Nelder-Mead
    control=c(  # configuration of Nelder-Mead
        maxit=100, # maximum iterations
        reltol=1e-8, # response tolerance over-one step
        alpha=1.0, # reflection factor
        beta=0.5, # contraction factor
        gamma=2.0)) # expansion factor

# summarize the results
# coordinate of the minimum
print(result$par)
# function response of the minimum
print(result$value)
# the number of function calls performed
print(result$counts)

# display the function as a contour plot
x <- seq(-3, 3, length.out=100)
y <- seq(-3, 3, length.out=100)
z <- rosenbrock(expand.grid(x, y))
contour(x, y, matrix(log10(z), length(x)), 
        xlab="x", ylab="y")
# draw the optima
points(result$par[1], result$par[2], col="red", pch=19)
rect(result$par[1]-0.2, result$par[2]-0.2,
    result$par[1]+0.2, result$par[2]+0.2,
    lwd=2) 
```

### 3.3 Gradient Descent
Gradient descent is the optimization method for __unconstrained 
nonlinear function__ using the *first-order derivative*. When 
it is used for function mamximization, it may be referred to as 
gradient ascent.

- Steepest descent is an extension that performs a line search on 
the line of gradient to locate the optimal neighboring point (
determining the optimal step or steepest step).
- Batch gradient descent (BGD) is an extension where the cost 
function and its derivative are computed as the summed error on 
a collection of training examples.
- For stochastic gradient descent (SGD, or Online gradient descent), 
the cost function and derivative are computed for each training 
example.

The following example illustrates the usage of gradient descent 
algorithm in R for solving a two-dimensional nonlinear 
optimization function.

```r
# define a 2D basin function, optima at c(0,0)
basin <- function(x) x[1]^2 + x[2]^2

# define the derivative function for the basin function
deriv <- function(x) {
    c(2*x[1], 2*x[2])
}

# definition of the gradient descent
gradient_descent <- function(func, derv, start, step=.05, tol=1e-8) {
    pt1 <- start
    grdnt <- derv(pt1)
    pt2 <- c(pt1[1] - step*grdnt[1], pt1[2] - step*grdnt[2])
    while (abs(func(pt1) - func(pt2)) > tol) {
        pt1 <- pt2
        grdnt <- derv(pt1)
        pt2 <- c(pt1[1] - step*grdnt[1], pt1[2] - step*grdnt[2])
        print(func(pt2)) # print progress
    }
    pt2 # return the last point
}

# find the minimum of the basin
result <- gradient_descent(
    basin, # function to be minimized
    deriv, # the gradient
    c(runif(1,-3,3), runif(1,-3,3)), # start point
    step=0.05, # step size
    tol=1e-8) # relative tolerance for one-step

# display the summary of the result
print(result) # coordinate of the function minimum
print(basin(result)) # the function minimum value

# display the function as a contour plot
x <- seq(-3, 3, length.out=100)
y <- seq(-3, 3, length.out=100)
z <- basin(expand.grid(x, y))
contour(x, y, matrix(z, length(x)), 
        xlab="x", ylab="y")
# draw the optimal point
points(result[1], result[2], col="red", pch=19)
rect(result[1]-.2, result[2]-.2, 
    result[1]+.2, result[2]+.2,
    lwd=2)
```

### 3.4 Conjugate Gradient

Conjugate Gradient Method is a first-order derivative optimization 
method for multidimensional nonlinear unconstrained functions. It 
is related to other first-order derivative optimization algorithms 
such as Gradient Descent and Steepest Descent.

The information processing objective of the technique is to locate 
the extremum of a function. From a starting position, the method first 
computes the gradient to locate the direction of steepest descent, then 
performs a line search to locate the optimum step size (alpha). The 
method then repeats the process of computing the steepest direction, 
computes direction of the search, and performing a line search to 
locate the optimum step size. A parameter beta defines the direction 
update rule based on the gradient and can be computed using one of a 
number of methods.

The difference between Conjugate Gradient and Steepest Descent is that 
it uses conjugate directions rather than local gradients to move 
downhill towards the function minimum, which can be very efficient.

The example provides a code listing of the Conjugate Gradient method 
in R solving a two-dimensional nonlinear optimization function.

```r
# definition of the 2D Rosenbrock function, optima is at (1,1)
rosenbrock <- function(v) { 
(1 - v[1])^2 + 100 * (v[2] - v[1]*v[1])^2
}

# definition of the gradient of the 2D Rosenbrock function
derivative <- function(v) {
    c(-400 * v[1] * (v[2] - v[1]*v[1]) - 2 * (1 - v[1]), 
        200 * (v[2] - v[1]*v[1]))
}

# locate the minimum of the function using the Conjugate Gradient method
result <- optim(
    c(runif(1,-3,3), runif(1,-3,3)), # start at a random position
    rosenbrock, # the function to minimize
    derivative, # no function gradient 
    method="CG", # use the Conjugate Gradient method
    control=c( # configure Conjugate Gradient
        maxit=100, # maximum iterations of 100
        reltol=1e-8, # response tolerance over-one step
        type=2)) # use the Polak-Ribiere update method

# summarize results
print(result$par) # the coordinate of the minimum
print(result$value) # the function response of the minimum
print(result$counts) # the number of function calls performed

# display the function as a contour plot
x <- seq(-3, 3, length.out=100)
y <- seq(-3, 3, length.out=100)
z <- rosenbrock(expand.grid(x, y))
contour(x, y, matrix(log10(z), length(x)), xlab="x", ylab="y")
# draw the optima as a point
points(result$par[1], result$par[2], col="red", pch=19)
rect(result$par[1]-0.2, result$par[2]-0.2, 
    result$par[1]+0.2, result$par[2]+0.2, 
    lwd=2)
```

### 3.5 BFGS

BFGS is an optimization method for multidimensional nonlinear unconstrained 
functions.

BFGS belongs to the family of quasi-Newton (Variable Metric) optimization 
methods that make use of both first-derivative (gradient) and 
second-derivative (Hessian matrix) based information of the function being 
optimized. More specifically, it is a quasi-Newton method which means that 
it approximates the second-order derivative rather than compute it directly. 
It is related to other quasi-Newton methods such as the DFP Method, 
Broyden’s method and the SR1 Method.

Two popular extension of BFGS is L-BFGS (Limited Memory BFGS) which has 
lower memory resource requirements and L-BFGS-B (Limited Memory Boxed BFGS) 
which extends L-BFGS and imposes box constraints on the method.

The information processing objective of the BFGS Method is to locate the 
extremum of a function.

This is achieved by iteratively building up a good approximation of the 
inverse Hessian matrix. Given an initial starting position, it prepares 
an approximation of the Hessian matrix (square matrix of second-order 
partial derivatives). It then repeats the process of computing the search 
direction using the approximated Hessian, computes an optimum step size 
using a Line Search then, updates the position, and updates the 
approximation of the Hessian. The method for updating the Hessian each 
iteration is called the BFGS rule which insures the updated matrix 
is positive definite.

The example provides a code listing of the BFGS method in R solving a 
two-dimensional nonlinear optimization function.

```r
# definition of the 2D Rosenbrock function, optima is at (1,1)
rosenbrock <- function(v) { 
    (1 - v[1])^2 + 100 * (v[2] - v[1]*v[1])^2
}

# definition of the gradient of the 2D Rosenbrock function
derivative <- function(v) {
    c(-400 * v[1] * (v[2] - v[1]*v[1]) - 2 * (1 - v[1]), 
        200 * (v[2] - v[1]*v[1]))
}

# locate the minimum of the function using the BFGS method
result <- optim(
    c(runif(1,-3,3), runif(1,-3,3)), # start at a random position
    rosenbrock, # the function to minimize
    derivative, # no function gradient 
    method="BFGS", # use the BFGS method
    control=c( # configure BFGS
    maxit=100, # maximum iterations of 100
    reltol=1e-8)) # response tolerance over-one step

# summarize results
print(result$par) # the coordinate of the minimum
print(result$value) # the function response of the minimum
print(result$counts) # the number of function calls performed

# display the function as a contour plot
x <- seq(-3, 3, length.out=100)
y <- seq(-3, 3, length.out=100)
z <- rosenbrock(expand.grid(x, y))
contour(x, y, matrix(log10(z), length(x)), xlab="x", ylab="y")
# draw the optima as a point
points(result$par[1], result$par[2], col="red", pch=19)
rect(result$par[1]-0.2, result$par[2]-0.2, 
    result$par[1]+0.2, result$par[2]+0.2, 
    lwd=2)
```
