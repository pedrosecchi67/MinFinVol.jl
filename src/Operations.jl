export linear_interpolation!, face_average!, upwind_values!, Neumann_conditions_impose!, volume_gradient!, face_gradient!, flux_integral!

"""
Linearly interpolate values from volumes to faces

* `msh`: the mesh at hand
* `at_volumes`: array of values at volume centers (input)
* `at_faces`: array of values at face centers (output)
"""
function linear_interpolation!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}) where N
    v1=0.0
    v2=0.0

    for (i, (fc_adj, fc_dists)) in enumerate(zip(msh.face_connectivity, msh.face_dists))
        if fc_adj[1]==fc_adj[2] && fc_adj[1]==0
            throw(
                error(
                    "linear_interpolation!:ERROR:unconnected face detected"
                )
            )
        end

        v1=(fc_adj[1]==0 ? 0.0 : at_volumes[fc_adj[1]])
        v2=(fc_adj[2]==0 ? 0.0 : at_volumes[fc_adj[2]])

        at_faces[i]=(v1*fc_dists[2]+v2*fc_dists[1])/(fc_dists[1]+fc_dists[2])
    end
end

"""
Average volume values at face centers.
When a face is connected to only one volume, it incorporates the value of that volume

* `msh`: the mesh at hand
* `at_volumes`: array of values at volume centers (input)
* `at_faces`: array of values at face centers (output)
"""
function face_average!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}) where N
    for (fc, cons) in enumerate(msh.face_connectivity)
        if (cons[1]==0 || cons[2]==0)
            if (cons[1]==cons[2])
                throw(
                    error(
                        "face_average!:ERROR:face $(fc) not connected to any volumes"
                    )
                )
            end

            if cons[1]==0
                at_faces[fc]=at_volumes[cons[2]]
            else
                at_faces[fc]=at_volumes[cons[1]]
            end
        else
            at_faces[fc]=(at_volumes[cons[1]]+at_volumes[cons[2]])/2
        end
    end
end

"""
Fetch values at faces from upstream volumes

* `msh`: the mesh at hand
* `at_volumes`: values at volume centers (input)
* `at_faces`: values at face centers (output)
* `u`, `v`, ... and any other set of velocities **at faces** with which to determine
    flow direction
"""
function upwind_values!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}, vels...) where N
    nv=0.0

    for (fc, fccons) in enumerate(msh.face_connectivity)
        if fccons[1]==0
            if fccons[2]==0
                throw(
                    error(
                        "upwind_values!:ERROR:unconnected face detected"
                    )
                )
            end

            at_faces[fc]=at_volumes[fccons[2]]
        elseif fccons[2]==0
            at_faces[fc]=at_volumes[fccons[1]]
        else
            nv=0.0

            for (n, v) in zip(msh.face_normals[fc], vels)
                nv+=n*v[fc]
            end

            if nv>=0.0
                at_faces[fc]=at_volumes[fccons[1]]
            else
                at_faces[fc]=at_volumes[fccons[2]]
            end
        end
    end
end

"""
Impose Neumann conditions for certain faces

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
"""
function Neumann_conditions_impose!(msh::FinVolMesh{N}, 
    at_volumes::Vector{Float64}, at_faces::Vector{Float64}, 
    Neumann_conditions::Tuple{Vector{Int64}, Vector{Vector{Float64}}}) where N
    ccell=0

    for (fc, nc) in zip(Neumann_conditions...)
        ccell=msh.face_connectivity[fc][1]

        if ccell==0
            ccell=msh.face_connectivity[fc][2]
        end

        if ccell==0
            throw(
                error(
                    "Neumann_conditions_impose!:ERROR:face disconnected from any cells"
                )
            )
        end

        at_faces[fc]=at_volumes[ccell]

        for (nd, d) in enumerate(nc)
            at_faces[fc]+=d*(msh.face_centers[fc][nd]-msh.volume_centers[ccell][nd])
        end
    end
end

"""
Green-Gauss gradient at cell centroids

* `msh`: the mesh at hand
* `at_volumes`: gradients at cell centroids (output)
* `at_faces`: values at faces (input)
* `dim`: the dimension along which to obtain the gradient
"""
function volume_gradient!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}, dim::Int64) where N
    at_volumes.=0.0

    for (i, (cons, shat)) in enumerate(zip(msh.face_connectivity, msh.face_normals))
        if cons[1]!=0
            at_volumes[cons[1]]+=shat[dim]*at_faces[i]
        end

        if cons[2]!=0
            at_volumes[cons[2]]-=shat[dim]*at_faces[i]
        end
    end

    for (i, v) in enumerate(msh.volume_volumes)
        at_volumes[i]/=v
    end
end

"""
Gradients at faces (averaged from finite differences with adjacent volumes)

* `msh`: the mesh at hand
* `at_volumes`: array of values at volume centers (input)
* `at_faces`: array of values at face centers (input and output)
* `dim`: spatial dimension for gradient
"""
function face_gradient!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, at_faces::Vector{Float64}, dim::Int64) where N
    for (fc, (cons, shat, dists)) in enumerate(zip(msh.face_connectivity, msh.face_normal_unit, msh.face_dists))
        if (cons[1]==0 || cons[2]==0)
            if (cons[1]==cons[2])
                throw(
                    error(
                        "face_average!:ERROR:face $(fc) not connected to any volumes"
                    )
                )
            end

            if cons[1]==0
                at_faces[fc]=(at_volumes[cons[2]]-at_faces[fc])*shat[dim]/dists[2]
            else
                at_faces[fc]=(at_faces[fc]-at_volumes[cons[1]])*shat[dim]/dists[1]
            end
        else
            at_faces[fc]=((at_volumes[cons[2]]-at_faces[fc])*shat[dim]/dists[2]+
                (at_faces[fc]-at_volumes[cons[1]])*shat[dim]/dists[1])/2
        end
    end
end

"""
Accumulate fluxes from faces into variations at cell centroids

Example:

Given `f1, f2, f3` in a three-dimensional mesh,
`du/dt=-1/v∫∇⋅fdΩ` is calculated at cell centroids and stored at the given vectors

* `msh`: the mesh at hand
* `at_volumes`: flux integral according to Gauss's theorem (output)
* `u, v, w...` and remaining arguments: fluxes at each spatial dimension, as vectors for face values (input) 
"""
function flux_integral!(msh::FinVolMesh{N}, at_volumes::Vector{Float64}, fs...) where N
    at_volumes.=0.0

    for (i, (cons, n)) in enumerate(zip(msh.face_connectivity, msh.face_normals))
        if cons[1]!=0
            for (j, f) in enumerate(fs)
                at_volumes[cons[1]]-=n[j]*f[i]
            end
        end

        if cons[2]!=0
            for (j, f) in enumerate(fs)
                at_volumes[cons[2]]+=n[j]*f[i]
            end
        end
    end

    for (i, v) in enumerate(msh.volume_volumes)
        at_volumes[i]/=v
    end
end
