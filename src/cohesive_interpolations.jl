abstract type SurfaceInterpolation{dim,shape,order,ip} <: Interpolation{dim,shape,order} end

Ferrite.getnbasefunctions(ip::SurfaceInterpolation{2, RefCube, order, Lagrange{1, RefCube, order}}) where {order} = 2*Ferrite.getnbasefunctions(ip.ip_base)
Ferrite.vertexdof_indices(ip::SurfaceInterpolation{2, RefCube, order, Lagrange{1, RefCube, order}}) where {order} = ((1,),(2,),(order-1+3,),(order-1+4,)) # TODO we should be able to do it via calling `vertexdof_indices` of `ip_base` + arithmetic
Ferrite.facedof_interior_indices(ip::SurfaceInterpolation{2, RefCube, order, Lagrange{1, RefCube, order}}) where {order} = (ntuple(i->2+i,order-1), ntuple(i->(order-1+4)+i,order-1)) # TODO we should be able to do it via calling `facedof_indices` of `ip_base` + arithmetic 

struct JumpInterpolation{dim,shape,order,ip} <: SurfaceInterpolation{dim,shape,order,ip}
    ip_base::ip
end

nvertices(ip::Interpolation) = length(Ferrite.vertexdof_indices(ip))
nvertexdofs(ip::Interpolation) = sum(length.(Ferrite.vertexdof_indices(ip)); init=0)
nfacedofs(ip::Interpolation) = sum(length.(Ferrite.facedof_indices(ip)); init=0)
ncelldofs(ip::Interpolation) = length(Ferrite.celldof_interior_indices(ip))

JumpInterpolation(ip::Interpolation{dim,shape,order}) where {dim,shape,order} = JumpInterpolation{dim+1,shape,order,typeof(ip)}(ip)

function Ferrite.value(ip::JumpInterpolation, i::Int, ξ::Vec{dim_base}) where dim_base
    (; ip_base) = ip
    @assert Ferrite.getdim(ip_base) == dim_base
    0 < i <= Ferrite.getnbasefunctions(ip) || throw(ArgumentError("no shape function $i for interpolation $ip"))

    ndofs_base = Ferrite.getnbasefunctions(ip_base)
    if i <= ndofs_base
        return -Ferrite.value(ip_base, i, ξ)
    else 
        return Ferrite.value(ip_base, i-ndofs_base, ξ)
    end
end

struct MidPlaneInterpolation{dim,shape,order,ip} <: SurfaceInterpolation{dim,shape,order,ip}
    ip_base::ip
end

MidPlaneInterpolation(ip::Interpolation{dim,shape,order}) where {dim,shape,order} = MidPlaneInterpolation{dim+1,shape,order,typeof(ip)}(ip)

# no need to handle edge dofs, as underlying interpolation is <= 2D
function Ferrite.value(ip::MidPlaneInterpolation, i::Int, ξ::Vec{dim_base}) where dim_base
    (; ip_base) = ip
    @assert Ferrite.getdim(ip_base) == dim_base
    0 < i <= getnbasefunctions(ip) || throw(ArgumentError("no shape function $i for interpolation $ip"))
    ndofs_base = Ferrite.getnbasefunctions(ip_base)
    if i <= ndofs_base
        return 0.5Ferrite.value(ip_base, i, ξ)
    else 
        return 0.5Ferrite.value(ip_base, i-ndofs_base, ξ)
    end
end
