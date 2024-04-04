abstract type CohesiveCell{refshape} <: Ferrite.AbstractCell{refshape} end

struct CohesiveQuadrilateral <: CohesiveCell{RefQuadrilateral} nodes::NTuple{4, Int} end
struct CohesiveQuadraticQuadrilateral <: CohesiveCell{RefQuadrilateral} nodes::NTuple{6, Int} end

struct CohesiveWedge <: CohesiveCell{RefPrism} nodes::NTuple{6, Int} end
struct CohesiveQuadraticWedge <: CohesiveCell{RefPrism} nodes::NTuple{12, Int} end

struct CohesiveHexahedron <: CohesiveCell{RefHexahedron} nodes::NTuple{8, Int} end
struct CohesiveQuadraticHexahedron <: CohesiveCell{RefHexahedron} nodes::NTuple{18, Int} end

Ferrite.vertices(c::Union{CohesiveQuadrilateral, CohesiveQuadraticQuadrilateral}) = (c.nodes[1], c.nodes[2], c.nodes[3], c.nodes[4])
Ferrite.faces(c::Union{CohesiveQuadrilateral, CohesiveQuadraticQuadrilateral}) = ((c.nodes[1], c.nodes[2]), (c.nodes[3],c.nodes[4]))
# TODO: will FaceValues normals be correct with these definitions?

Ferrite.default_interpolation(::Type{CohesiveQuadrilateral}) = JumpInterpolation(Lagrange{1,RefCube,1}())
Ferrite.default_interpolation(::Type{CohesiveQuadraticQuadrilateral}) = JumpInterpolation(Lagrange{1,RefCube,2}())
Ferrite.default_interpolation(::Type{CohesiveWedge}) = JumpInterpolation(Lagrange{2,RefTetrahedron,1}())
Ferrite.default_interpolation(::Type{CohesiveQuadraticWedge}) = JumpInterpolation(Lagrange{2,RefTetrahedron,2}())
Ferrite.default_interpolation(::Type{CohesiveHexahedron}) = JumpInterpolation(Lagrange{2,RefCube,1}())

Ferrite.cell_to_vtkcell(::Type{CohesiveQuadrilateral}) = VTKCellTypes.VTK_QUAD
Ferrite.cell_to_vtkcell(::Type{CohesiveQuadraticQuadrilateral}) = VTKCellTypes.VTK_BIQUADRATIC_QUAD
Ferrite.cell_to_vtkcell(::Type{CohesiveWedge}) = VTKCellTypes.VTK_WEDGE
Ferrite.cell_to_vtkcell(::Type{CohesiveQuadraticWedge}) = VTKCellTypes.VTK_QUADRATIC_WEDGE
Ferrite.cell_to_vtkcell(::Type{CohesiveHexahedron}) = VTKCellTypes.VTK_HEXAHEDRON

Ferrite.nodes_to_vtkorder(cell::CohesiveQuadrilateral) = [cell.nodes[i] for i in [1,2,4,3]]

#=
cohesive_celltypes = Dict{DataType, String}(CohesiveCell{2,4,2}  => "CohesiveQuadrilateral",
                                            CohesiveCell{2,6,2}  => "CohesiveQuadraticQuadrilateral",
                                            CohesiveCell{3,6,2}  => "CohesiveTetrahedron",
                                            CohesiveCell{3,12,2} => "CohesiveQuadraticTetrahedron",
                                            CohesiveCell{3,8,2}  => "CohesiveHexahedron")
=#
