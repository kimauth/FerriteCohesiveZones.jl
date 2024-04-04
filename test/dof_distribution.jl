dim = 2
grid = generate_grid(Quadrilateral, (3,1))

# Ferrite default grid:
# 5 ____ 6 ____ 7 ____ 8
# |      |      |      |
# |      |      |      |
# 1 ____ 2 ____ 3 ____ 4

# make middle cell zero-thickness:
grid.nodes[2:3] = [Node(mean(n->n.x, grid.nodes[2:3])) for _ in 1:2]
grid.nodes[6:7] = [Node(mean(n->n.x, grid.nodes[6:7])) for _ in 1:2]

cohesivecell = CohesiveQuadrilateral((3,7,2,6))

cohesive_grid = Grid([grid.cells[1], grid.cells[3], cohesivecell], grid.nodes)

ip_quads = Lagrange{RefQuadrilateral, 1}()^dim
ip_base = Lagrange{RefLine, 1}()
ip_jump = JumpInterpolation(ip_base)^dim
ip_mid = MidPlaneInterpolation(ip_base)^dim

dh = DofHandler(cohesive_grid)
sdh_quads = SubDofHandler(dh, Set(1:2))
add!(sdh_quads, :u, ip_quads)
sdh_coh = SubDofHandler(dh, Set(3))
add!(sdh_coh, :u, ip_jump)
close!(dh)

cc = CellCache(dh)

reinit!(cc, 1)
@test celldofs(cc) == collect(1:8)
reinit!(cc, 2)
@test celldofs(cc) == collect(9:16)
reinit!(cc, 3) # cohesive cell
@test celldofs(cc) == [9,10,15,16,3,4,5,6]

