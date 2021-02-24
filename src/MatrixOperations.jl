using SparseArrays

"""
Function to generate an empty sparse matrix with the pre-allocated
structure as derivable from the adjacencies between volumes in a mesh
"""
function empty_spmatrix(msh::FinVolMesh{N}) where N

	n_non_border = 0

	for cons in msh.face_connectivity
		if cons[1] != 0 && cons[2] != 0
			n_non_border += 1
		end
	end

	n_entries = 2 * n_non_border + length(msh.volumes)

	lins = zeros(Int64, n_entries)
	cols = zeros(Int64, n_entries)
	data = ones(Float64, n_entries)

	entry_counter = 1

	for (i, cons) in enumerate(msh.face_connectivity)
		if cons[1] != 0 && cons[2] != 0
			lins[entry_counter] = cons[1]
			cols[entry_counter] = cons[2]

			entry_counter += 1

			lins[entry_counter] = cons[2]
			cols[entry_counter] = cons[1]

			entry_counter += 1
		end
	end

	for i = 1:length(msh.volumes)
		lins[entry_counter] = i
		cols[entry_counter] = i

		entry_counter += 1
	end

	return sparse(lins, cols, data)

end

export empty_spmatrix

