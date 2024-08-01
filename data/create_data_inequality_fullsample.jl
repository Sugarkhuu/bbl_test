#------------------------------------------------------------------------------
# Creates final data set for model calibration/estimation using raw data
#------------------------------------------------------------------------------

using XLSX, Statistics, DataFrames, CSV

cd("bayer_rayfair_comment/BBL_Inequality_Replication/data")

TL = range(1954.75, stop=2019.75, step=0.25)

OBSMAT = fill!(zeros(maximum(size(TL)), 41), NaN)

# FRED data
data_raw_fred = XLSX.readdata("data_raw/data_raw.xlsx", "Data!B6:M293")
data = convert(Matrix{Float64}, data_raw_fred)
timeline = range(1948, stop=2019.75, step=0.25)

nom_cons_d = data[:, 1]
nom_cons_s = data[:, 2]
nom_cons_nd = data[:, 3]
nom_inv = data[:, 4]
nom_gov_spend = data[:, 5]
gdp_defl = data[:, 6]
nom_gdp = data[:, 7]
pop_temp = data[:, 8]
tot_hours = data[:, 9]
nom_wage_comp = data[:, 10]
cons_defl = data[:, 11]
inv_defl = data[:, 12]

# smooth hp trend of population to solve "best levels problem"
# see section 3.1.3 in https://sites.google.com/site/pfeiferecon/Pfeifer_2013_Observation_Equations.pdf
function hp_filter(y, lambda)
    T = size(y, 1)
    FilterMat = zeros(T, T)

    FilterMat[1, 1:3] = [1 + lambda, -2 * lambda, lambda]
    FilterMat[2, 1:4] = [-2 * lambda, 1 + 5 * lambda, -4 * lambda, lambda]

    for i = 3 : T - 2
        FilterMat[i, i-2 : i+2] = [lambda, -4*lambda, 1 + 6 * lambda, -4 * lambda, lambda]
    end

    FilterMat[T-1, T-3:T] = [lambda, -4 * lambda, 1 + 5 * lambda, -2 * lambda]
    FilterMat[T, T-2:T] = [lambda, -2 * lambda, 1 + lambda]

    trend = FilterMat \ y
    cycle = y - trend

    return trend, cycle
end

pop = hp_filter(pop_temp, 10000)[1]

real_cons_pc = (nom_cons_nd + nom_cons_d + nom_cons_s) ./ (cons_defl .* pop)
real_inv_pc = nom_inv ./ (inv_defl .* pop)
real_gdp_sum_pc = (nom_cons_nd + nom_cons_d + nom_cons_s + nom_inv + nom_gov_spend) ./ (gdp_defl .* pop)
real_wage_comp = nom_wage_comp ./ gdp_defl

hours_pc = tot_hours ./ pop
inflation = [NaN; diff(log.(gdp_defl))]

# Tax progressivity
# uses intermediate data file created by /tax_progressivity/tax_prog.m
toptaxraw_fn = CSV.read("data_final/tax_progressivity.csv", DataFrame, header=0)[!, 2]
timeline_toptax_fn = range(1946.75, stop=2017.75, step=0.25)
timeline_toptax_annual_fn = range(1946.75, stop=2017.75, step=1)

toptax_fn = fill!(zeros(size(timeline_toptax_fn)), NaN)
toptax_fn[1:4:end] .= toptaxraw_fn

# Shadow federal funds rate
WuXiaShadowRate = XLSX.readdata("data_raw/WuXiaShadowRate.xlsx", "Data!C590:C673")[:]
shadow_rate_monthly = convert(Vector{Float64}, WuXiaShadowRate)
FedFundsRate = XLSX.readdata("data_raw/FedFundsRate.xlsx", "Data!B12:B797")[:]
fedfunds_rate_monthly = convert(Vector{Float64}, FedFundsRate)
# replace ZLB period by shadow rate
fedfunds_rate_monthly[655:738] .= shadow_rate_monthly

# take quarterly average
fedfunds_rate_quarterly = map(mean, Iterators.partition(fedfunds_rate_monthly, 3))

# Convert to quarterly gross rate
shadow_rate_quarterly = 1 .+ (fedfunds_rate_quarterly ./ 400)

# Timeline
timeline_shadow_rate = range(1954.5, step=0.25, stop=2019.75)

# Uncertainty
income_uncertainty = XLSX.readdata("data_raw/BayerEtAl2019_income_uncertainty.xlsx", "Data!B6:B126")
timeline_uncertainty = range(1983, step=0.25, stop=2013)

# Income and Wealth shares
inequality_raw = XLSX.readdata("data_raw/WID_Data_16022023-171829.xlsx", "Data!C2:D74")
timeline_inequality = range(1947.75, step=0.25, stop=2019.75)
timeline_inequality_annual = range(1947.75, step=1, stop=2019.75)

WStop10 = fill!(similar(timeline_inequality), NaN)
IStop10 = fill!(similar(timeline_inequality), NaN)

WStop10[1:4:end] .= inequality_raw[:, 2]
IStop10[1:4:end] .= inequality_raw[:, 1]

# additional data for targeted moments
# fixed_assets_annual = XLSX.readdata("data_raw/FixedAssets1947.xlsx", "Sheet1!A2:A74")
fixed_assets_annual = XLSX.readdata("data_raw/K1TTOTL1ES000.xlsx", "Data!B34:B106")
# depreciation_annual = XLSX.readdata("data_raw/M1TTOTL1ES000.xlsx", "Data!B34:B106")

# depr_rate = mean(depreciation_annual[end-30:end]./fixed_assets_annual[end-30:end])

timeline_fixedassets_quarterly = range(1947.75, stop=2019.75, step=0.25)
timeline_fixedassets_annual = range(1947.75, stop=2019.75, step=1)

fixed_assets = fill!(zeros(size(timeline_fixedassets_quarterly)), NaN)
fixed_assets[1:4:end] .= fixed_assets_annual
KYratio = fixed_assets[2:end] ./ (nom_cons_nd + nom_cons_d + nom_cons_s + nom_inv + nom_gov_spend) # [2:end] to start in 1948Q1

debt_share_annual = XLSX.readdata("data_raw/FYPUGDA188S.xlsx", "Data!B20:B92")
timeline_debtshare_quarterly = range(1947.75, stop=2019.75, step=0.25)
timeline_debtshare_annual = range(1947.75, stop=2019.75, step=1)

debt_share = fill!(zeros(size(timeline_debtshare_quarterly)), NaN)
debt_share[1:4:end] .= debt_share_annual

Gshare = nom_gov_spend ./ (nom_cons_nd + nom_cons_d + nom_cons_s + nom_inv + nom_gov_spend)

# Compute targeted moments
nanmean(x) = mean(filter(!isnan,x))
TargetedMoments = fill!(zeros(5), NaN)

TargetedMoments[1]= nanmean(KYratio[(timeline_fixedassets_quarterly[2:end] .>= TL[1]) .& (timeline_fixedassets_quarterly[2:end] .<= TL[end])])*4.0

TargetedMoments[2] = nanmean(debt_share[(timeline_debtshare_quarterly .>= TL[1]) .& (timeline_debtshare_quarterly .<= TL[end])])*4.0/100.0
# TargetedMoments[2]= nanmean(debt_share[(timeline_debt_MV .>= TL[1]) .& (timeline_debt_MV .<= TL[end])])

TargetedMoments[3] = nanmean(toptax_fn[(timeline_toptax_fn .>= TL[1]) .& (timeline_toptax_fn .<= TL[end])])

TargetedMoments[4] = nanmean(WStop10[(timeline_inequality .>= TL[1]) .& (timeline_inequality .<= TL[end])])*100.0

TargetedMoments[5] = nanmean(Gshare[(timeline .>= TL[1]) .& (timeline .<= TL[end])])*100.0

dfMoments = DataFrame(Target=["K-to-Y share", "Debt share", "Mean tax prog.", "Top10 wealth share", "G-to-Y share"], Value=TargetedMoments)
dfMoments[!, :Value] = round.(dfMoments[:, :Value], digits=2)
CSV.write("data_final/data_moments.csv", dfMoments)
# println(TargetedMoments)

# Adjust timing of R, today's observation is t+1 in model
timeline_shadow_rate = range(timeline_shadow_rate[1] + 0.25, step=0.25, stop=timeline_shadow_rate[end] + 0.25)

# In Log Differences
C_obs = diff(log.(real_cons_pc[(timeline .>= TL[1]-0.25) .& (timeline .<= TL[end])])) .- mean(diff(log.(real_cons_pc[(timeline .>= TL[1]-0.25) .& (timeline .<= TL[end])])))
I_obs = diff(log.(real_inv_pc[(timeline .>= TL[1]-0.25) .& (timeline .<= TL[end])])) .- mean(diff(log.(real_inv_pc[(timeline .>= TL[1]-0.25) .& (timeline .<= TL[end])])))
Y_sum_obs = diff(log.(real_gdp_sum_pc[(timeline .>= TL[1]-0.25) .& (timeline .<= TL[end])])) .- mean(diff(log.(real_gdp_sum_pc[(timeline .>= TL[1]-0.25) .& (timeline .<= TL[end])])))
W_comp_obs = diff(log.(real_wage_comp[(timeline .>= TL[1]-0.25) .& (timeline .<= TL[end])])) .- mean(diff(log.(real_wage_comp[(timeline .>= TL[1]-0.25) .& (timeline .<= TL[end])]))) 

# In Log Deviations
TOPT_fn_obs = (log.(toptax_fn[(timeline_toptax_fn .>= TL[1]) .& (timeline_toptax_fn .<= TL[end])])) .- nanmean(log.(toptax_fn[(timeline_toptax_fn .>= TL[1]) .& (timeline_toptax_fn .<= TL[end])]))
W10_obs = (log.(WStop10[(timeline_inequality .>= TL[1]) .& (timeline_inequality .<= TL[end])])) .- nanmean(log.(WStop10[(timeline_inequality .>= TL[1]) .& (timeline_inequality .<= TL[end])]))
I10_obs = (log.(IStop10[(timeline_inequality .>= TL[1]) .& (timeline_inequality .<= TL[end])]) .- nanmean(log.(IStop10[(timeline_inequality .>= TL[1]) .& (timeline_inequality .<= TL[end])])))
sigma_obs = log.(income_uncertainty[(timeline_uncertainty .>= TL[1]) .& (timeline_uncertainty .<= TL[end])]) .- mean(log.(income_uncertainty[(timeline_uncertainty .>= TL[1]) .& (timeline_uncertainty .<= TL[end])]))
R_obs = log.(shadow_rate_quarterly[(timeline_shadow_rate .>= TL[1]) .& (timeline_shadow_rate .<= TL[end])]) .- mean(log.(shadow_rate_quarterly[(timeline_shadow_rate .>= TL[1]) .& (timeline_shadow_rate .<= TL[end])]))
Pi_obs = inflation[(timeline .>= TL[1]) .& (timeline .<= TL[end])] .- mean(inflation[(timeline .>= TL[1]) .& (timeline .<= TL[end])])
N_obs = (log.(hours_pc[(timeline .>= TL[1]) .& (timeline .<= TL[end])]) .- mean(log.(hours_pc[(timeline .>= TL[1]) .& (timeline .<= TL[end])])))

# Save in csv file
OBSMAT[(TL .>= timeline[1]) .& (TL .<= timeline[end]), 1] = Y_sum_obs
OBSMAT[(TL .>= timeline[1]) .& (TL .<= timeline[end]), 2] = I_obs
OBSMAT[(TL .>= timeline[1]) .& (TL .<= timeline[end]), 3] = C_obs
OBSMAT[(TL .>= timeline[1]) .& (TL .<= timeline[end]), 4] = N_obs
OBSMAT[(TL .>= timeline[1]) .& (TL .<= timeline[end]), 5] = W_comp_obs
OBSMAT[(TL .>= timeline_shadow_rate[1]) .& (TL .<= timeline_shadow_rate[end]), 6] = R_obs
OBSMAT[(TL .>= timeline[1]) .& (TL .<= timeline[end]), 7] = Pi_obs
OBSMAT[(TL .>= timeline_uncertainty[1]) .& (TL .<= timeline_uncertainty[end]), 8] = sigma_obs
OBSMAT[(TL .>= timeline_inequality[1]) .& (TL .<= timeline_inequality[end]), 9] = W10_obs
OBSMAT[(TL .>= timeline_inequality[1]) .& (TL .<= timeline_inequality[end]), 10] = I10_obs
OBSMAT[(TL .>= timeline_toptax_fn[1]) .& (TL .<= timeline_toptax_fn[end]), 11] = TOPT_fn_obs

T = DataFrame(OBSMAT[:,[1,2,3,4,5,6,7,8,9,10,11]], 
              [:Ygrowth, :Igrowth, :Cgrowth, :N, :wgrowth, :RB, :pi, :sigma2, :TOP10Wshare, :TOP10Ishare, :tauprog])

CSV.write("data_final/bbl_data_inequality_fullsample.csv", T)

exit()