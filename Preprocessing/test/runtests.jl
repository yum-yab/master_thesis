using Base: return_types
using Test
using Preprocessing

using NetCDF
using Statistics
using NumericalIntegration

function get_hus_sampleset()
  return ncread("./testdata/hus_sampleset.nc", "hus")
end

function get_plevs_sampleset()
  return ncread("./testdata/hus_sampleset.nc", "plev")
end

function get_iwv_sampleset()
  return ncread("./testdata/prw_sampleset.nc", "prw")
end

function myiwv(pressure_levels, hus_data)
  

  function get_avg_humidity(i)
    if i == length(pressure_levels)
      return pressure_levels[end]/2, hus_data[end]/2
    else 
      return (pressure_levels[i] + pressure_levels[i+1])/2, (hus_data[i] + hus_data[i+1])/2
    end
  end

  result = 0.

  layer_number = length(pressure_levels) - 1

  for i in 1:layer_number
    layer_thickness = pressure_levels[i] - pressure_levels[i + 1]
    hus_avg = (hus_data[i] + hus_data[i + 1])/2
    result += layer_thickness/9.806 * hus_avg
  end 

  return result

end

function calculate_stats_on_difference(ivt_fun, iwv_data, hus_data, pressure_levels, stat_funs...)

  (lonsize, latsize, levsize, timesize) = size(hus_data)
  wind_vector = fill(1., levsize)

  diffs_percentage = [ 
    abs(ivt_fun(pressure_levels, hus_data[lo, la, :, t], wind_vector, wind_vector).eastward_integral - iwv_data[lo, la, t])/iwv_data[lo, la, t] 
    for lo in 1:lonsize 
    for la in 1:latsize 
    for t in 1:timesize
  ]

  results = zeros(length(stat_funs))

  for (i, stat_fun) in enumerate(stat_funs)
    results[i] = stat_fun(diffs_percentage)
  end

  return  results
end


@testset "IVT function test" begin
  
  # load the actual data from CMIP6 simulation
  hus_sample_data = get_hus_sampleset()
  iwv_sample_data = get_iwv_sampleset()
  pressure_levels = get_plevs_sampleset()
  
  # create a artificial wind vector 
  wind_vector = fill(1., size(hus_sample_data, 3))

  # calculate IVT 
  ivt_result = IVT.ivt_of_column(pressure_levels, hus_sample_data[1, 1, :, 1], wind_vector, wind_vector)
  
  (median_diff, mean_diff, maximum_diff) = calculate_stats_on_difference(IVT.ivt_of_column, iwv_sample_data, hus_sample_data, pressure_levels, median, mean, maximum) 
  println("Statistics of differnences of simulation iwv and my own integral of the same humidity data")
  println("Median diff: $median_diff\tMean Diff: $mean_diff\tMaximum diff: $maximum_diff")
  @test median_diff < .1
  @test mean_diff < .2

  # @test ivt_result.eastward_integral == iwv_sample_data[1, 1, 1]
end
