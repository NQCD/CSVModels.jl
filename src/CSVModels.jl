module CSVModels
using BSplineKit
using NQCModels
using Parameters: Parameters
using NQCModels

#==
Since CSV files come in many shapes and sizes, leave importing them into a matrix to the user. 

CSV model needs an input matrix of shape (2,n), where the first column of the matrix denotes the position, the second half the potential value in Hartree. 
The model assumes that the positions are evenly spaced, forming a uniform grid, for faster fitting. 
==#

Parameters.@with_kw struct CSVModel_1D <: NQCModels.AdiabaticModels.AdiabaticModel
    # Input matrix
    potential_matrix::Matrix{Float64}
    # build interpolation and derivative object once, then evaluate with functions
    potential_function=interpolate(potential_matrix[:,1],potential_matrix[:,2],BSplineOrder(3))
    derivative=Derivative(1)*potential_function
end



NQCModels.ndofs(model::CSVModel_1D)=Int(size(model.potential_matrix)[2]//2)

function NQCModels.potential(model::CSVModel_1D, R::AbstractArray)
    return(model.potential_function.(R))
end

function NQCModels.derivative!(model::CSVModel_1D, D::AbstractArray, R::AbstractArray)
    D=model.derivative.(R)
    return(D)
end

CSVModel_1D(x)=CSVModel_1D(potential_matrix=x) # Definition shortcut. 

export CSVModel_1D

end
