# Values for cohesive zone elements
abstract type CohesiveCellValues{dim,T,refshape} <: CellValues{dim,T,refshape} end

struct CohesiveVectorValues{dim,T<:Real,refshape<:Ferrite.AbstractRefShape,M} <: CohesiveCellValues{dim,T,refshape}
    N::Matrix{Vec{dim,T}}
    dNdξ::Matrix{Tensor{2,dim,T,M}}
    dNdx::Matrix{Tensor{2,dim,T,M}}
    detJdA::Vector{T}
    M::Matrix{T}
    dMdξ::Matrix{Vec{dim,T}}
    # deliberately abstractly typed, weights can be inferred correctly
    qr::QuadratureRule{<:Any,refshape,T}
    R::Vector{Tensor{2,dim,T,M}}
    # The following fields are deliberately abstract -- they are never used in
    # performace critical code, just stored here for convenience.
    func_interp::Interpolation{dim,refshape}
    geo_interp::Interpolation{dim,refshape}
end

Ferrite.FieldTrait(::Type{<:CohesiveVectorValues}) = Ferrite.VectorValued()

function CohesiveVectorValues(
    qr::QuadratureRule,
    func_interpol::SurfaceInterpolation,
    geom_interpol::SurfaceInterpolation = MidPlaneInterpolation(func_interpol.ip_base)
)
    return CohesiveVectorValues(Float64, qr, func_interpol, geom_interpol)
end

function CohesiveVectorValues(
    ::Type{T},
    quad_rule::QuadratureRule{dim_base,shape},
    func_interpol::SurfaceInterpolation{dim,shape},
    geom_interpol::SurfaceInterpolation{dim,shape} = MidPlaneInterpolation(func_interpol.ip_base)
) where {T, dim_base, dim, shape} 

    @assert Ferrite.getdim(func_interpol) == Ferrite.getdim(geom_interpol)
    @assert Ferrite.getrefshape(func_interpol) == Ferrite.getrefshape(geom_interpol) == shape
    n_qpoints = length(getweights(quad_rule))


    # Function interpolation
    n_func_basefuncs = getnbasefunctions(func_interpol) * dim
    N    = fill(zero(Vec{dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)
    dNdξ = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)
    dNdx = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)

    covar_base = fill(zero(Tensor{2,dim,T}) * T(NaN), n_qpoints)

    # Geometry interpolation
    n_geom_basefuncs = getnbasefunctions(geom_interpol)
    M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints)
    dMdξ = fill(zero(Vec{dim,T}) * T(NaN), n_geom_basefuncs, n_qpoints)

    for (qp, ξ) in enumerate(quad_rule.points)
        basefunc_count = 1
        for basefunc in 1:getnbasefunctions(func_interpol)
            dNdξ_temp, N_temp = Ferrite.gradient(ξ -> Ferrite.value(func_interpol, basefunc, ξ), ξ, :all)
            for comp in 1:dim
                N_comp = zeros(T, dim)
                N_comp[comp] = N_temp
                N[basefunc_count, qp] = Vec{dim,T}((N_comp...,))

                dN_comp = zeros(T, dim, dim)
                dN_comp[comp, 1:dim-1] = dNdξ_temp
                dNdξ[basefunc_count, qp] = Tensor{2,dim,T}((dN_comp...,))
                basefunc_count += 1
            end
        end
        for basefunc in 1:n_geom_basefuncs
            dMdξ_temp, M[basefunc, qp] = Ferrite.gradient(ξ -> Ferrite.value(geom_interpol, basefunc, ξ), ξ, :all)
            dM_comp = zeros(T, dim)
            dM_comp[1:dim-1] = dMdξ_temp
            dMdξ[basefunc, qp] = Vec{dim,T}((dM_comp...,))
        end
    end

    detJdA = fill(T(NaN), n_qpoints)

    MM = Tensors.n_components(Tensors.get_base(eltype(dNdx)))

    CohesiveVectorValues{dim,T,shape,MM}(N, dNdξ, dNdx, detJdA, M, dMdξ, quad_rule, covar_base, func_interpol, geom_interpol)
end

@inline getdetJdA(cv::CohesiveCellValues, q_point::Int) = cv.detJdA[q_point]
@inline Ferrite.getdetJdV(cv::CohesiveCellValues, q_point::Int) = getdetJdA(cv, q_point)
@inline get_rotation(cv::CohesiveCellValues, qp::Int) = cv.R[qp]

function Ferrite.reinit!(
    cv::CohesiveCellValues{dim},
    x::AbstractVector{Vec{dim,T}},
) where {dim,T}

    n_geom_basefuncs = Ferrite.getngeobasefunctions(cv)
    n_func_basefuncs = getnbasefunctions(cv)
    @assert length(x) == n_geom_basefuncs

    @inbounds for qp in 1:length(cv.qr.weights)
        w = cv.qr.weights[qp]
        fecv_J = zero(Tensor{2,dim})
        for j in 1:n_geom_basefuncs
            fecv_J += x[j] ⊗ cv.dMdξ[j, qp]
        end
        n = mid_plane_normal(fecv_J) # normal vector of the cohesive element
        e = basevec(n, dim)
        fecv_J = fecv_J + n/norm(n) ⊗ e 
        detJ = det(fecv_J)
        detJ > 0.0 || throw(ArgumentError("det(J) is not positive: det(J) = $(detJ)"))
        cv.detJdA[qp] = detJ * w
        # compute rotation matrix
        cv.R[qp] = rotation_matrix(fecv_J)
        # compute dNdx
        Jinv = inv(fecv_J)
        for i in 1:n_func_basefuncs
            cv.dNdx[i, qp] = cv.dNdξ[i, qp] ⋅ Jinv
        end
    end
end

mid_plane_normal(dMdx::Tensor{2,2,T}) where T = Vec{2,T}((-dMdx[2,1], dMdx[1,1]))
mid_plane_normal(dMdx::Tensor{2,3,T}) where T = dMdx[:,1] × dMdx[:,2]

function rotation_matrix(dMdx::Tensor{2,dim,T}) where {dim,T}
    R = dMdx
    for d = 1:dim
        v = R[:, d]
        N = Tensor{2,dim,T}((i,j)->i==j ? (i==d ? inv(norm(v)) : one(T) ) : zero(T))
        R = R ⋅ N
    end
    return R
end



