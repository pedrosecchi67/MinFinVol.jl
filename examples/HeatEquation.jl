using MinFinVol
using Plots

# defining a mesh
Lx=1.0
Ly=1.0

nx=20
ny=20

xs=LinRange(0.0, Lx, nx)
ys=LinRange(0.0, Ly, ny)

# mesh generation
msh=FinVolMesh{2}()

ptinds=zeros(Int64, nx, ny)
for (i, x) in enumerate(xs)
    for (j, y) in enumerate(ys)
        isedge=(i==1 || i==nx || j==1 || j==ny)
        w=((isedge & i==1) ? 1.1 : 1.0)

        ptinds[i, j]=add_point!(msh, [x, y]; fixed=isedge, weight=w)
    end
end

left_face=Vector{Int64}()
lower_face=Vector{Int64}()
upper_face=Vector{Int64}()
right_face=Vector{Int64}()

volume_inds=zeros(Int64, nx-1, ny-1)
for i=1:(nx-1)
    for j=1:(ny-1)
        newfcs=[
            add_face!(msh, [ptinds[i, j], ptinds[i, j+1]]; univocal=true),
            add_face!(msh, [ptinds[i, j+1], ptinds[i+1, j+1]]; univocal=true),
            add_face!(msh, [ptinds[i+1, j+1], ptinds[i+1, j]]; univocal=true),
            add_face!(msh, [ptinds[i+1, j], ptinds[i, j]]; univocal=true)
        ]

        if i==1
            push!(left_face, newfcs[1])
        elseif i==(nx-1)
            push!(right_face, newfcs[3])
        end

        if j==1
            push!(lower_face, newfcs[4])
        elseif j==(ny-1)
            push!(upper_face, newfcs[2])
        end

        volume_inds[i, j]=add_volume!(msh, newfcs)
    end
end

process_geometry!(msh)
# mesh_smoothing!(msh)

Dirichlet_left=1.0
Dirichlet_right=0.0
Dirichlet_lower=0.0
# Dirichlet_upper=0.0
x_Neumann_upper=zeros(Float64, length(upper_face))
y_Neumann_upper=ones(Float64, length(upper_face))
Neumann_upper=[collect(t) for t in zip(x_Neumann_upper, y_Neumann_upper)]

# pre-allocations
_at_faces=zeros(Float64, nx*(ny-1)+ny*(nx-1))
# _at_volumes=zeros(Float64, (nx-1)*(ny-1))

_qx=zeros(Float64, size(_at_faces)...)
_qy=zeros(Float64, size(_at_faces)...)

function getres!(r, u)
    # variable values at faces
    face_average!(msh, u, _at_faces)

    # applying Dirichlet conditions
    _at_faces[left_face].=Dirichlet_left
    _at_faces[right_face].=Dirichlet_right
    _at_faces[lower_face].=Dirichlet_lower
    # _at_faces[upper_face].=Dirichlet_upper
    Neumann_conditions_impose!(msh, u, _at_faces, (upper_face, Neumann_upper))

    # obtain derivatives
    _qx.=_at_faces
    face_gradient!(msh, u, _qx, 1)

    _qy.=_at_faces
    face_gradient!(msh, u, _qy, 2)

    flux_integral!(msh, r, _qx, _qy)
end

u=zeros(Float64, (nx-1)*(ny-1))
r=zeros(Float64, size(u)...)

dt=1e-4

for i=1:2000
    getres!(r, u)

    u.-=(dt.*r)
end

plotlyjs()

x_plot=[msh.volume_centers[(nx-1)*(i-1)+j][1] for i=1:(nx-1), j=1:(ny-1)]
y_plot=[msh.volume_centers[(nx-1)*(i-1)+j][2] for i=1:(nx-1), j=1:(ny-1)]

p=surface(x_plot, y_plot, u[volume_inds])
display(p)
