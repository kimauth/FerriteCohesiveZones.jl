abstract type SurfaceInterpolation{dim,shape,order,ip} <: Interpolation{dim,shape,order} end

Ferrite.getnbasefunctions(ip::SurfaceInterpolation) = Ferrite.getnbasefunctions(ip.ip_base)*2
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

struct JumpInterpolation{dim,shape,order,ip} <: SurfaceInterpolation{dim,shape,order,ip}
    ip_base::ip
end

JumpInterpolation(ip::Interpolation{dim,shape,order}) where {dim,shape,order} = JumpInterpolation{dim+1,shape,order,typeof(ip)}(ip)
# no need to handle edge dofs, as underlying interpolation is <= 2D
function Ferrite.value(ip::JumpInterpolation, i::Int, ξ::Vec{dim_base}) where dim_base
    (; ip_base) = ip
    @assert Ferrite.getdim(ip_base) == dim_base
    0 < i <= getnbasefunctions(ip) || throw(ArgumentError("no shape function $i for interpolation $ip"))

    vertex_dofs = nvertices(ip_base)*Ferrite.nvertexdofs(ip_base)
    face_dofs = length(Ferrite.faces(ip_base))*Ferrite.nfacedofs(ip_base)
    cell_dofs = Ferrite.ncelldofs(ip_base)
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

struct MidPlaneInterpolation{dim,shape,order,ip} <: SurfaceInterpolation{dim,shape,order,ip}
    ip_base::ip
end

MidPlaneInterpolation(ip::Interpolation{dim,shape,order}) where {dim,shape,order} = MidPlaneInterpolation{dim+1,shape,order,typeof(ip)}(ip)

# no need to handle edge dofs, as underlying interpolation is <= 2D
function Ferrite.value(ip::MidPlaneInterpolation, i::Int, ξ::Vec{dim_base}) where dim_base
    (; ip_base) = ip
    @assert Ferrite.getdim(ip_base) == dim_base
    0 < i <= getnbasefunctions(ip) || throw(ArgumentError("no shape function $i for interpolation $ip"))

    vertex_dofs = nvertices(ip_base)*Ferrite.nvertexdofs(ip_base)
    face_dofs = length(Ferrite.faces(ip_base))*Ferrite.nfacedofs(ip_base)
    cell_dofs = Ferrite.ncelldofs(ip_base)
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