# Define test parameters
fcycle = 1.0 / 365.0
t = 1:365  # One year of daily data
overtones = 2

# Test fit_mean function
E_mean, F_mean = mean_matrices(t)
@test size(E_mean) == (365, 1)
@test all(E_mean[:, 1] .== 1.0)

# Test seasonal_matrices with fit_mean = true
E1, F1 = seasonal_matrices(fcycle, t, overtones, true)
@test size(E1) == (365, (overtones + 1) * 2 + 1)
@test all(E1[:, 1] .== 1.0)

# Test seasonal_matrices with fit_mean = false
E2, F2 = seasonal_matrices(fcycle, t, overtones, false)
@test size(E2) == (365, (overtones + 1) * 2)

# Check that the sine and cosine values are correctly computed
ω = 2π * fcycle
for i = 1:overtones+1
    @test all(E1[:, i+1] .≈ sin.(i * ω * t))
    @test all(E1[:, i+overtones+2] .≈ cos.(i * ω * t))
    @test all(E2[:, i] .≈ sin.(i * ω * t))
    @test all(E2[:, i+overtones+1] .≈ cos.(i * ω * t))
end