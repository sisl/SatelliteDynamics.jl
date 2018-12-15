using SatelliteDynamics

load_eop(:C04_14)

# println(collect(keys(EOP.data)))

println(UT1_UTC(37666))