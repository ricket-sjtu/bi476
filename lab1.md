# Lab 1: Introductory Math For Biostatistics Using R
___

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
