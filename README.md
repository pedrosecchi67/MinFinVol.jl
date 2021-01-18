# MinFinVol.jl

MinFinVol.jl is a minimalist package for assisting the programming of finite volume formulations in Julia.

To install it, use:

```
]add https://github.com/pedrosecchi67/MinFinVol.jl
```

## Generating meshes

Mutable struct to define a mesh:

* `ndim_space`, `ndim_mesh`: numbers of dimensions for the space in which the mesh is inserted and for the mesh
    itself, respectively. A 2D, surface mesh defined as the surface of a sphere, for example, would yield `ndim_mesh=2, ndim_space=3`
* `points`: array of arrays of floats for corner points
* `point_connectivity`: array of arrays of ints, indicating points with which each point shares
    a face. Used for mesh smoothing
* `point_fixed`: array of bools. Indicates whether or not each point in the mesh should be fixed when
    performing mesh smoothing
* `point_weight`: weights for each point when performing elliptic mesh smoothing. Array of floats
* `faces`: array of arrays of ints for faces (point indices)
* `face_centers`: array of arrays of floats for face centers
* `face_normals`: array of arrays of floats indicating face normals (with norm equal to the face's area)
* `face_normal_unit`: the same as `face_normals`, but with unit norm
* `face_connectivity`: array of arrays of ints indicating, respectively, the volume behind the face and the volume ahead of it
    (or more specifically their indices)
* `face_dists`: distances between the face and its adjacent volumes, in the same format as `face_connectivity`
* `volumes`: array of arrays of ints for faces (face indices)
* `volume_centers`: array of arrays of floats with positions for volume centers
* `volume_volumes`: volume of each finite volume (array of floats)
```
mutable struct FinVolMesh{N}
    ndim_space::Int64
    ndim_mesh::Int64
    points::Vector{Vector{Float64}}
    point_connectivity::Vector{Vector{Int64}}
    point_fixed::Vector{Bool}
    point_weight::Vector{Float64}
    faces::Vector{Vector{Int64}}
    face_centers::Vector{Vector{Float64}}
    face_normals::Vector{Vector{Float64}}
    face_normal_unit::Vector{Vector{Float64}}
    face_connectivity::Vector{Vector{Int64}}
    face_dists::Vector{Vector{Float64}}
    volumes::Vector{Vector{Int64}}
    volume_centers::Vector{Vector{Float64}}
    volume_volumes::Vector{Float64} 
end
``` 

### Instantiating a Mesh

Constructor for an empty mesh:

* `ndim_mesh`: number of dimensions for the mesh, if differnet from space dimensionality.
   A 2D, surface mesh defined as the surface of a sphere, for example, would yield `ndim_mesh=2, ndim_space=3`
```
function FinVolMesh{N}(; ndim_mesh::Union{Int64, Nothing}=nothing) where N 
```

### Adding Points, Faces and Cells

Add a point to a mesh:

* `msh`: the mesh
* `pt`: a vector with the point at hand
* `fixed`: a boolean identifying whether the point should be fixed during mesh smoothing
* `weight`: a weight for mesh smoothing
* `geps`: if set to a non-zero value, is set as a minimum threshold for merging with other pre-existing points

* return: the new point's index
```
function add_point!(msh::FinVolMesh{N}, pt::Vector{Float64}; 
    fixed::Bool=false, weight::Float64=1.0, geps::Float64=0.0) where N
```

Add a face to a mesh:

* `msh`: the mesh at hand
* `pts`: a vector of point indices
* `univocal`: if true (default), allows for merging with equal, previously set faces

* return: the new face's index

```
function add_face!(msh::FinVolMesh{N}, pts::Vector{Int64}; univocal::Bool=true) where N
```

Add volume to a mesh:

* `msh`: the mesh
* `fcs`: vector of face indices

* return: an index for the volume
```
function add_volume!(msh::FinVolMesh{N}, fcs::Vector{Int64}) where N
```

Function to pre-process the geometry of a finite volume mesh and obtain pertinent geometrical properties, after smoothing and geometry alterations are performed:

* `msh`: the mesh at hand
```
function process_geometry!(msh::FinVolMesh{N}) where N
```

## Mesh Smoothing

Function to apply elliptic equation mesh smoothing. Must be used after `process_geometry!`, which is also automatically ran after changes are made

* `msh`: the mesh at hand
* `niter`: number of iterations (defaults to 100)
* `relaxation_parameter`: relaxation parameter for iteration
```
function mesh_smoothing!(msh::FinVolMesh{N}; niter::Int64=100, relaxation_parameter::Float64=1.0) where N
```

## Solving Equations

### Juggling Cell and Face Values

Linearly interpolate values from volumes to faces:

* `msh`: the mesh at hand
* `at_volumes`: array of values at volume centers (input)
* `at_faces`: array of values at face centers (output)
```
function linear_interpolation!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}) where N
```

Average volume values at face centers:

* `msh`: the mesh at hand
* `at_volumes`: array of values at volume centers (input)
* `at_faces`: array of values at face centers (output)
```
function face_average!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}) where N 
```

Fetch values at faces from upstream volumes:

* `msh`: the mesh at hand
* `at_volumes`: values at volume centers (input)
* `at_faces`: values at face centers (output)
* `u`, `v`, ... and any other set of velocities **at faces** with which to determine
    flow direction
```
function upwind_values!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}, vels...) where N
```

Obtain face values from the imposition of Neumann conditions:

* `msh`: the mesh at hand
* `at_volumes`: values at volume centers, from which to apply variations
* `at_faces`: vector with values at all face centers
* `Neumann_conditions`: tuple with an array of face indexes
    and an array of float arrays with the imposed gradients at those faces

    Example:

    ```
    Neumann_conditions=(
        [1, 3, 4], # face indexes
        [ # values of gradients at faces
            [1.0, 0.0], # ∇f=(1, 0)
            [0.0, 1.0], # ∇f=(0, 1)
            [0.5, 0.5] # ∇f=(1/2, 1/2)
        ]
    )
    ```
```
function Neumann_conditions_impose!(msh::FinVolMesh{N}, 
    at_volumes::Vector{Float64}, at_faces::Vector{Float64}, 
    Neumann_conditions::Tuple{Vector{Int64}, Vector{Vector{Float64}}}) where N
```

### Differentiation Through the Mesh

Obtain Green-Gauss gradient at cell centroids:

* `msh`: the mesh at hand
* `at_volumes`: gradients at cell centroids (output)
* `at_faces`: values at faces (input)
* `dim`: the dimension along which to obtain the gradient
```
function volume_gradient!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}, dim::Int64) where N
```

Obtain gradients at faces (averaged from finite differences with adjacent volumes):

* `msh`: the mesh at hand
* `at_volumes`: array of values at volume centers (input)
* `at_faces`: array of values at face centers (input and output)
* `dim`: spatial dimension for gradient
```
function face_gradient!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}, dim::Int64) where N
```

### Applying the Divergence Theorem

Accumulate fluxes from faces into variations at cell centroids:

Example:

Given `f1, f2, f3` in a three-dimensional mesh,
`du/dt=-1/v∫∇⋅fdΩ` is calculated at cell centroids and stored at the given vectors

* `msh`: the mesh at hand
* `at_volumes`: flux integral according to Gauss's theorem (output)
* `u, v, w...` and remaining arguments: fluxes at each spatial dimension, as vectors for face values (input) 
```
function flux_integral!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, fs...) where N
```

## Example: Forward Euler Scheme for Laplace Equation

(For the complete code, see `examples/HeatEquation.jl`)

```
_qx=zeros(Float64, size(_at_faces)...)
_qy=zeros(Float64, size(_at_faces)...)

function getres!(r, u)
    # variable values at faces
    face_average!(msh, u, _at_faces)

    # applying Dirichlet conditions
    _at_faces[left_face].=Dirichlet_left
    _at_faces[right_face].=Dirichlet_right
    _at_faces[lower_face].=Dirichlet_lower

    # applying Neumann conditions
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
```
