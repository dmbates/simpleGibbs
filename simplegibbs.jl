## Gibbs sampler for function:

## f(x,y) = x x^2 \exp(-xy^2 - y^2 + 2y - 4x)

## using conditional distributions:

## x|y \sim Gamma(3, y^2 +4)
## y|x \sim Normal(\frac{1}{1+x}, \frac{1}{2(1+x)})

## Direct Julia translations of the R sampler
## Rgibbs <- function(N,thin) {
##     mat <- matrix(0,ncol=2,nrow=N)
##     x <- 0
##     y <- 0
##     for (i in 1:N) {
##         for (j in 1:thin) {
##             x <- rgamma(1,3,y*y+4)
##             y <- rnorm(1,1/(x+1),1/sqrt(2*(x+1)))
##         }
##         mat[i,] <- c(x,y)
##     }
##     mat
## }


## First using the Julia functions that call libRmath
require("Rmath.jl")
using Distributions

function JGibbs1(N::Integer, thin::Integer)
    mat = Array(Float64, (N, 2))
    x = 0.; y = 0.
    for i in 1:N
        for j in 1:thin
            x = rgamma(1,3,1/(y*y + 4))[1] # rgamma uses rate not scale
            y = rnorm(1, 1/(x+1),1/sqrt(2(x + 1)))[1]
        end
        mat[i,:] = [x,y]
    end
    mat
end

## With direct calls to the underlying C functions in libRmath
function JGibbs2(N::Int, thin::Int)
    mat = Array(Float64, (N, 2))
    x = 0.; y = 0.
    for i in 1:N
        for j in 1:thin
            x = ccall((:rgamma, :libRmath), Float64,
                      (Float64, Float64), 3., 1/(y*y + 4))
            y = ccall((:rnorm, :libRmath), Float64,
                      (Float64, Float64), 1/(x+1), 1/sqrt(2(x + 1)))
        end
        mat[i,:] = [x,y]
    end
    mat
end

## Calling the Julia random number generators
function JGibbs3(N::Int, thin::Int)
    mat = Array(Float64, (N, 2))
    x = 0.; y = 0.
    G = Gamma(3.)
    for i in 1:N
        for j in 1:thin
            x = rand(G) / (y*y + 4.)
            y = randn()/sqrt(2.x+2.) + 1./(x+1.)
        end
        mat[i,:] = [x,y]
    end
    mat
end

## Sampling vectors
function JGibbs4(N::Int, thin::Int)
    X = Array(Float64,(N,)); Y = Array(Float64, (N,)) 
    x = 0.; y = 0.
    gg = rand(Gamma(3.), N*thin)
    nn = randn(N*thin)
    ii = 0
    for i in 1:N
        for j in 1:thin
            ii += 1
            x = gg[ii] / (y*y + 4.)
            y = nn[ii]/sqrt(2.x+2.) + 1./(x+1.)
        end
        X[i] = x; Y[i] = y
    end
    hcat(X,Y)
end

## Distributed versions - keeping the results as a distributed array
dJGibbs3a(N::Int, thin::Int) = DArray(I->JGibbs3(length(I[1]),thin), (N, 2))
## Converting the results to an array controlled by the parent process
dJGibbs3b(N::Int, thin::Int) = convert(Array{Float64,2}, dJGibbs3a(N, thin))

## Timings
#println("JGibbs1: $([@elapsed JGibbs1(20000, 200) for i=1:10])")
#println("JGibbs2: $([@elapsed JGibbs2(20000, 200) for i=1:10])")
#println("JGibbs3: $([@elapsed JGibbs3(20000, 200) for i=1:10])")
#println("dJGibbs3a: $([@elapsed dJGibbs3a(20000, 200) for i=1:10])")
#println("dJGibbs3b: $([@elapsed dJGibbs3b(20000, 200) for i=1:10]')")
