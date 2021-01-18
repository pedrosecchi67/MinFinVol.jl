@assert begin
    msh=FinVolMesh{2}()

    newptind1=add_point!(msh, [0.0, 0.0])
    newptind2=add_point!(msh, [1.0, 0.0])
    newptind3=add_point!(msh, [0.0, 1.0])

    @assert newptind1==1
    @assert newptind2==2
    @assert newptind3==3

    newfaceind1=add_face!(msh, [newptind1, newptind2])
    newfaceind2=add_face!(msh, [newptind2, newptind3])
    newfaceind3=add_face!(msh, [newptind3, newptind1])

    @assert newfaceind1==1
    @assert newfaceind2==2
    @assert newfaceind3==3

    newvol=add_volume!(msh, [newfaceind1, newfaceind2, newfaceind3])

    @assert newvol==1

    process_geometry!(msh)
    @assert isapprox(msh.volume_volumes[1], 0.5)

    mesh_smoothing!(msh)
    @assert isapprox(msh.volume_volumes[1], 0.0)

    true
end
