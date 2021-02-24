@assert begin
    msh=FinVolMesh{1}()

    nnds=10

    pts=[add_point!(msh, [x]) for x in LinRange(0.0, 1.0, nnds)]

    fcs=[
        add_face!(msh, [pt]) for pt in pts
    ]

    vols=[
        add_volume!(msh, [fcs[i], fcs[i+1]]) for i=1:(length(fcs)-1)
    ]

    process_geometry!(msh)

    u_nodes=collect(LinRange(0.0, 1.0, nnds))
    u_cells=(u_nodes[1:(end-1)].+u_nodes[2:end])./2

    u_faces=Vector{Float64}(undef, nnds)

    face_average!(msh, u_cells, u_faces)
    
    eqs=isapprox.(u_faces, u_nodes)

    @assert all(eqs[2:(end-1)])

    u_faces[1]=0.0
    u_faces[end]=1.0

    grads=zeros(Float64, nnds-1)
    volume_gradient!(msh, grads, u_faces, 1)

    @assert isapprox(grads, ones(Float64, nnds-1); atol=1e-3)

    flux_ints=zeros(Float64, nnds-1)

    flux_integral!(msh, flux_ints, u_faces)
    @assert isapprox(flux_ints, -1.0*ones(Float64, nnds-1); atol=1e-3)

	spm = empty_spmatrix(msh)

	for i = 1:(nnds - 1)
		if i < (nnds - 1)
			@assert spm[i, i + 1] == spm[i + 1, i]
			@assert spm[i, i + 1] == 1.0
		end

		@assert spm[i, i] == 1.0
	end

    true
end
