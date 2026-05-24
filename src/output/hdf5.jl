# hdf5.jl
#
# Thin bridge between HDF5.jl's group API and StructArrayTables. HDF5.jl
# provides per-dataset reads and (low-level) creation-order iteration; we add
# a `_group_to_namedtuple` helper that preserves the C++ side's column order,
# and expose `h5read(StructArray, file, group)` as the public ad-hoc reader.

function _creation_order_keys(g::HDF5.Group)
    n = length(g)
    [HDF5.API.h5l_get_name_by_idx(g, ".", HDF5.API.H5_INDEX_CRT_ORDER,
        HDF5.API.H5_ITER_INC, i, HDF5.API.H5P_DEFAULT) for i in 0:(n-1)]
end

"""
    _group_to_namedtuple(g::HDF5.Group) -> NamedTuple

Read every dataset in `g` (in creation order) into a NamedTuple keyed by the
raw HDF5 dataset names. Output is suitable for `StructArray(::NamedTuple)`.
"""
function _group_to_namedtuple(g::HDF5.Group)
    ks = _creation_order_keys(g)
    NamedTuple{Tuple(Symbol.(ks))}(Tuple(read(g[k]) for k in ks))
end

"""
    h5read(::Type{StructArray}, filename, groupname) -> StructArray

Read an HDF5 group into a `StructArray` with creation-order column preservation.
Column names are returned verbatim from HDF5 — see `load(outdir)` for the
high-level reader with unit attachment and suffix cleanup.
"""
function HDF5.h5read(::Type{StructArray}, filename::AbstractString, groupname::AbstractString)
    h5open(filename) do f
        StructArray(_group_to_namedtuple(f[groupname]))
    end
end

"""
    h5display(filename)

Print the HDF5 file's structure tree to stdout.
"""
function h5display(filename)
    h5open(filename) do f
        display(f)
    end
end
