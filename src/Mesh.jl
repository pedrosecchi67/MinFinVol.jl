export FinVolMesh, add_point!, add_face!, add_volume!, process_geometry!, mesh_smoothing!

using LinearAlgebra

"""
Mutable struct to define a mesh

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
"""
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
    volume_volumes::Vector{Float64} # duh
end

"""
Jacobian space volume
"""
function Jac_volume(Jac::Matrix{Float64})
    v=1.0
    a=1.0
    d=1.0
    
    for i=1:size(Jac, 2)
        for j=1:(i-1)
            Jac[:, i].-=dot(Jac[:, i], Jac[:, j]).*Jac[:, j]
        end

        locn=norm(Jac[:, i])
        Jac[:, i]./=(locn+1e-18)

        if i==size(Jac, 2)
            a=v
            d=locn
        end

        v*=locn/i
    end

    Jac[:, end]*a, v, d 
end

"""
Function to process the geometry of a finite volume mesh and obtain pertinent geometrical properties 

* `msh`: the mesh at hand
"""
function process_geometry!(msh::FinVolMesh{N}) where N
    for i=1:length(msh.points)
        msh.point_connectivity[i]=[]
    end

    for i=1:length(msh.faces)
        msh.face_centers[i].=0.0
        msh.face_normals[i].=0.0
        msh.face_normal_unit[i].=0.0
        msh.face_dists[i].=0.0

        msh.face_connectivity[i].=0
    end

    msh.volume_volumes.=0.0

    for i=1:length(msh.faces)
        for ptind in msh.faces[i]
            msh.face_centers[i].+=msh.points[ptind]
        end

        msh.face_centers[i]./=length(msh.faces[i])

        for j=1:length(msh.faces[i])
            p1=msh.faces[i][j]
            p2=msh.faces[i][j%length(msh.faces[i])+1]

            if !(p1 in msh.point_connectivity[p2])
                push!(msh.point_connectivity[p2], p1)
            end

            if !(p2 in msh.point_connectivity[p1])
                push!(msh.point_connectivity[p1], p2)
            end
        end
    end

    for i=1:length(msh.volumes)
        for fcind in msh.volumes[i]
            msh.volume_centers[i].+=msh.face_centers[fcind]
        end

        msh.volume_centers[i]./=length(msh.volumes[i])
    end

    Jac=Matrix{Float64}(undef, N, msh.ndim_mesh)

    for i=1:length(msh.volumes)
        for fc in msh.volumes[i]
            fcpts=msh.faces[fc]
            nfc=length(fcpts)

            freeside=1
            if msh.face_connectivity[fc][freeside]!=0
                freeside=2
            end

            if msh.face_connectivity[fc][freeside]!=0
                throw(
                    error(
                        "process_geometry!:ERROR:face $(fc) connected to more than two volumes"
                    )
                )
            end

            for j=1:(nfc-msh.ndim_mesh+1)
                for k=1:(msh.ndim_mesh-1)
                    kind=k+j

                    Jac[:, k].=msh.points[fcpts[kind]].-msh.points[fcpts[1]]
                end

                if freeside==1
                    Jac[:, msh.ndim_mesh].=msh.face_centers[fc].-msh.volume_centers[i]
                else
                    Jac[:, msh.ndim_mesh].=msh.volume_centers[i].-msh.face_centers[fc]
                end

                shat, v, d=Jac_volume(Jac)

                msh.face_normals[fc].+=shat
                msh.face_dists[fc][freeside]=d
                msh.face_connectivity[fc][freeside]=i

                msh.volume_volumes[i]+=v
            end

            if freeside==2
                msh.face_normals[fc]./=2
            end
        end
    end

    for fc=1:length(msh.faces)
        msh.face_normal_unit[fc].=msh.face_normals[fc]/norm(msh.face_normals[fc])
    end
end

"""
Constructor for an empty mesh

* `ndim_mesh`: number of dimensions for the mesh, if differnet from space dimensionality.
   A 2D, surface mesh defined as the surface of a sphere, for example, would yield `ndim_mesh=2, ndim_space=3`
"""
function FinVolMesh{N}(; ndim_mesh::Union{Int64, Nothing}=nothing) where N
    FinVolMesh{N}(
        N,
        something(ndim_mesh, N),
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        []
    )
end

"""
Add a point to a mesh

* `msh`: the mesh
* `pt`: a vector with the point at hand
* `fixed`: a boolean identifying whether the point should be fixed during mesh smoothing
* `weight`: a weight for mesh smoothing
* `geps`: if set to a non-zero value, is set as a minimum threshold for merging with other pre-existing points

* return: the new point's index
"""
function add_point!(msh::FinVolMesh{N}, pt::Vector{Float64}; fixed::Bool=false, weight::Float64=1.0, geps::Float64=0.0) where N
    @assert N==length(pt)

    if geps>0.0
        for (i, ept) in enumerate(msh.points)
            if norm(ept.-pt)<geps
                return i
            end
        end
    end

    push!(msh.points, pt)

    push!(msh.point_fixed, fixed)
    push!(msh.point_weight, weight)
    push!(msh.point_connectivity, [])

    length(msh.points)
end

function _is_face_eq(f1, f2)
    if length(f1)!=length(f2)
        return false
    end

    for p1 in f1
        if !(p1 in f2)
            return false
        end
    end

    true
end

"""
Add a face to a mesh

* `msh`: the mesh at hand
* `pts`: a vector of point indices
* `univocal`: if true (default), allows for merging with equal, previously set faces

* return: the new face's index
"""
function add_face!(msh::FinVolMesh{N}, pts::Vector{Int64}; univocal::Bool=true) where N
    for (i, fc) in enumerate(msh.faces)
        if _is_face_eq(fc, pts)
            return i
        end
    end

    push!(msh.faces, pts)

    push!(msh.face_centers, zeros(Float64, N))
    push!(msh.face_normals, zeros(Float64, N))
    push!(msh.face_normal_unit, zeros(Float64, N))
    push!(msh.face_connectivity, zeros(Int64, 2))
    push!(msh.face_dists, zeros(Float64, 2))

    length(msh.faces)
end

"""
Add volume to a mesh

* `msh`: the mesh
* `fcs`: vector of face indices

* return: an index for the volume
"""
function add_volume!(msh::FinVolMesh{N}, fcs::Vector{Int64}) where N
    push!(msh.volumes, fcs)

    push!(msh.volume_centers, zeros(Float64, N))
    push!(msh.volume_volumes, 0.0)

    length(msh.volumes)
end

"""
Function to apply elliptic equation mesh smoothing. Must be used after `process_geometry!`, which is
also automatically re-ran after changes are made

* `msh`: the mesh at hand
* `niter`: number of iterations (defaults to 100)
* `relaxation_parameter`: relaxation parameter for iteration
"""
function mesh_smoothing!(msh::FinVolMesh{N}; niter::Int64=100, relaxation_parameter::Float64=1.0) where N
    pos=zeros(Float64, N)
    w=0.0

    for nit=1:niter
        for (i, (pt, cons)) in enumerate(zip(msh.points, msh.point_connectivity))
            if !msh.point_fixed[i]
                w=msh.point_weight[i]
                pos.=pt.*w

                for c in cons
                    pos.+=msh.points[c].*msh.point_weight[c]
                    w+=msh.point_weight[c]
                end

                pos./=w

                msh.points[i].=relaxation_parameter.*pos.+(1.0-relaxation_parameter).*msh.points[i]
                
                w/=(length(cons)+1)
                msh.point_weight[i]=(1.0-relaxation_parameter)*msh.point_weight[i]+relaxation_parameter*w
            end
        end
    end

    process_geometry!(msh)
end
