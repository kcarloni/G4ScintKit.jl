
function read_sipm_voltage_trace( filename )
    d = h5read( filename, "sipm_voltage_trace" )
    V = d["voltages_V"] * 1u"V" # 512 × n_events in Julia (column-major)
    n = size(V, 2)
    StructArray(;
        event_id       = d["event_id"],
        t_min_ns       = d["t_min_ns"],
        t_max_ns       = d["t_max_ns"],
        time_bin_width_ns = d["time_bin_width_ns"],
        voltages       = [@view(V[:, i]) for i in 1:n],
    )
end
