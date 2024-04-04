abstract type SurfaceInterpolation{refshape,order,ip} <: ScalarInterpolation{refshape,order} end

Ferrite.getnbasefunctions(ip::SurfaceInterpolation) = Ferrite.getnbasefunctions(ip.ip_base)*2
#=
Ferrite.nvertexdofs(ip::SurfaceInterpolation) = Ferrite.nvertexdofs(ip.ip_base)
Ferrite.nedgedofs(ip::SurfaceInterpolation) = Ferrite.nfacedofs(ip.ip_base)
Ferrite.nfacedofs(ip::SurfaceInterpolation) = Ferrite.ncelldofs(ip.ip_base) 
Ferrite.ncelldofs(ip::SurfaceInterpolation) = 0 # cohesive elements never have dofs inside the cell

# stuff for handling higher order elements correctly
nvertices(::Interpolation{1,RefCube}) = 2
nvertices(::Interpolation{2,RefCube}) = 4
nvertices(::Interpolation{2,RefTetrahedron}) = 3
nvertices(::Interpolation{3,RefCube}) = 8
nvertices(::Interpolation{3,RefTetrahedron}) = 4
=#

##########################################
# adapting to Ferrite master
_ndofs(dof_indices::Tuple) = sum(idxs->length(idxs), dof_indices)
nvertexdofs(ip_base) = _ndofs(Ferrite.vertexdof_indices(ip_base))
nedgedofs(ip_base) = _ndofs(Ferrite.edgedof_interior_indices(ip_base))
nfacedofs(ip_base) = _ndofs(Ferrite.facedof_interior_indices(ip_base))
ncelldofs(ip_base) = _ndofs(Ferrite.celldof_interior_indices(ip_base))


function Ferrite.vertexdof_indices(ip::SurfaceInterpolation)
    (; ip_base) = ip
    vdofs = Ferrite.vertexdof_indices(ip_base)
    vdofs2 = map(v->map(d->d+nvertexdofs(ip_base), v), vdofs)
    return (vdofs..., vdofs2...)
end

function Ferrite.edgedof_interior_indices(ip::SurfaceInterpolation{3})
    (; ip_base) = ip
    edofs = Ferrite.facedof_interior_indices(ip_base)
    offset = nvertexdofs(ip_base)
    edofs1 = map(v->map(d->d + offset, v), edofs)
    edofs2 = map(v->map(d->d + offset + nfacedofs(ip_base), v), edofs)
    return (edofs1..., edofs2...)
end

function Ferrite.facedof_interior_indices(ip::SurfaceInterpolation{dim}) where dim
    (; ip_base) = ip
    fdofs = Ferrite.celldof_interior_indices(ip_base)
    offset = nvertexdofs(ip_base) + (dim==3 ? nedgedofs(ip_base) : 0)
    fdofs1 = map(v->map(d->d+offset, v), fdofs)
    fdofs2 = map(v->map(d->d+offset+ncelldofs(ip_base), v), fdofs)
    return (fdofs1, fdofs2)
end

Ferrite.adjust_dofs_during_distribution(ip::SurfaceInterpolation) = Ferrite.adjust_dofs_during_distribution(ip.ip_base)
##########################################
_promote_baserefshape(::Type{RefLine}) = RefQuadrilateral
_promote_baserefshape(::Type{RefTetrahedron}) = RefPrism
_promote_baserefshape(::Type{RefQuadrilateral}) = RefHexahedron

struct JumpInterpolation{refshape,order,ip} <: SurfaceInterpolation{refshape,order,ip}
    ip_base::ip
end

JumpInterpolation(ip::Interpolation{refshape,order}) where {refshape,order} = JumpInterpolation{_promote_baserefshape(refshape),order,typeof(ip)}(ip)
# no need to handle edge dofs, as underlying interpolation is <= 2D
function Ferrite.value(ip::JumpInterpolation, i::Int, ξ::Vec{dim_base}) where dim_base
    (; ip_base) = ip
    @assert Ferrite.getdim(ip_base) == dim_base
    0 < i <= getnbasefunctions(ip) || throw(ArgumentError("no shape function $i for interpolation $ip"))

    vertex_dofs = nvertexdofs(ip_base) # nvertices(ip_base)*Ferrite.nvertexdofs(ip_base)
    face_dofs = nfacedofs(ip_base) # length(Ferrite.faces(ip_base))*Ferrite.nfacedofs(ip_base)
    cell_dofs = ncelldofs(ip_base) # Ferrite.ncelldofs(ip_base)
    if i <= 2vertex_dofs
        if i <= vertex_dofs
            return -Ferrite.value(ip_base, i, ξ)
        else
            return Ferrite.value(ip_base, i-vertex_dofs, ξ)
        end
    elseif i - 2vertex_dofs <= 2face_dofs 
        if i-2vertex_dofs <= face_dofs 
            return -Ferrite.value(ip_base, i-vertex_dofs, ξ)
        else
            return Ferrite.value(ip_base, i-vertex_dofs-face_dofs, ξ)
        end
    elseif i - 2vertex_dofs - 2face_dofs <= 2cell_dofs 
        if i - 2vertex_dofs - 2face_dofs <= cell_dofs 
            return -Ferrite.value(ip_base, i-vertex_dofs - face_dofs, ξ)
        else
            return Ferrite.value(ip_base, i-vertex_dofs-face_dofs-cell_dofs, ξ)
        end
    end
end

struct MidPlaneInterpolation{refshape,order,ip} <: SurfaceInterpolation{refshape,order,ip}
    ip_base::ip
end

MidPlaneInterpolation(ip::Interpolation{refshape,order}) where {refshape,order} = MidPlaneInterpolation{_promote_baserefshape(refshape),order,typeof(ip)}(ip)

# no need to handle edge dofs, as underlying interpolation is <= 2D
function Ferrite.value(ip::MidPlaneInterpolation, i::Int, ξ::Vec{dim_base}) where dim_base
    (; ip_base) = ip
    @assert Ferrite.getdim(ip_base) == dim_base
    0 < i <= getnbasefunctions(ip) || throw(ArgumentError("no shape function $i for interpolation $ip"))

    vertex_dofs = nvertexdofs(ip_base) # nvertices(ip_base)*Ferrite.nvertexdofs(ip_base)
    face_dofs = nfacedofs(ip_base) # length(Ferrite.faces(ip_base))*Ferrite.nfacedofs(ip_base)
    cell_dofs = ncelldofs(ip_base) # Ferrite.ncelldofs(ip_base)
    if i <= 2vertex_dofs
        if i <= vertex_dofs
            return 0.5Ferrite.value(ip_base, i, ξ)
        else
            return 0.5Ferrite.value(ip_base, i-vertex_dofs, ξ)
        end
    elseif i - 2vertex_dofs <= 2face_dofs 
        if i-2vertex_dofs <= face_dofs 
            return 0.5Ferrite.value(ip_base, i-vertex_dofs, ξ)
        else
            return 0.5Ferrite.value(ip_base, i-vertex_dofs-face_dofs, ξ)
        end
    elseif i - 2vertex_dofs - 2face_dofs <= 2cell_dofs 
        if i - 2vertex_dofs - 2face_dofs <= cell_dofs 
            return 0.5Ferrite.value(ip_base, i-vertex_dofs - face_dofs, ξ)
        else
            return 0.5Ferrite.value(ip_base, i-vertex_dofs-face_dofs-cell_dofs, ξ)
        end
    end
end
