
struct CohesiveCell{dim,N,M} <: Ferrite.AbstractCell{dim,N,M}
    nodes::NTuple{N,Int}
end

const CohesiveQuadrilateral = CohesiveCell{2,4,2}
const CohesiveQuadraticQuadrilateral = CohesiveCell{2,6,2}

const CohesiveTetrahedron = CohesiveCell{3,6,2}
const CohesiveQuadraticTetrahedron = CohesiveCell{3,12,2}

const CohesiveHexahedron = CohesiveCell{3,8,2}

Ferrite.vertices(c::Union{CohesiveQuadrilateral, CohesiveQuadraticQuadrilateral}) = (c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4])
Ferrite.faces(c::Union{CohesiveQuadrilateral, CohesiveQuadraticQuadrilateral}) = ((c.nodes[1], c.nodes[2]), (c.nodes[3],c.nodes[4]))
# TODO: will FaceValues normals be correct with these definitions?

Ferrite.default_interpolation(::Type{CohesiveQuadrilateral}) = JumpInterpolation(Lagrange{1,RefCube,1}())
Ferrite.default_interpolation(::Type{CohesiveQuadraticQuadrilateral}) = JumpInterpolation(Lagrange{1,RefCube,2}())
Ferrite.default_interpolation(::Type{CohesiveTetrahedron}) = JumpInterpolation(Lagrange{2,RefTetrahedron,1}())
Ferrite.default_interpolation(::Type{CohesiveQuadraticTetrahedron}) = JumpInterpolation(Lagrange{2,RefTetrahedron,2}())
Ferrite.default_interpolation(::Type{CohesiveHexahedron}) = JumpInterpolation(Lagrange{2,RefCube,1}())

Ferrite.cell_to_vtkcell(::Type{CohesiveQuadrilateral}) = VTKCellTypes.VTK_QUAD
Ferrite.cell_to_vtkcell(::Type{CohesiveQuadraticQuadrilateral}) = VTKCellTypes.VTK_BIQUADRATIC_QUAD
Ferrite.cell_to_vtkcell(::Type{CohesiveTetrahedron}) = VTKCellTypes.VTK_WEDGE
Ferrite.cell_to_vtkcell(::Type{CohesiveQuadraticTetrahedron}) = VTKCellTypes.VTK_QUADRATIC_WEDGE
Ferrite.cell_to_vtkcell(::Type{CohesiveHexahedron}) = VTKCellTypes.VTK_HEXAHEDRON

Ferrite.nodes_to_vtkorder(cell::CohesiveQuadrilateral) = [cell.nodes[i] for i in [1,2,4,3]]

cohesive_celltypes = Dict{DataType, String}(CohesiveCell{2,4,2}  => "CohesiveQuadrilateral",
                                            CohesiveCell{2,6,2}  => "CohesiveQuadraticQuadrilateral",
                                            CohesiveCell{3,6,2}  => "CohesiveTetrahedron",
                                            CohesiveCell{3,12,2} => "CohesiveQuadraticTetrahedron",
                                            CohesiveCell{3,8,2}  => "CohesiveHexahedron")
