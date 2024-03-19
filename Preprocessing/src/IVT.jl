module IVT
  
  import NumericalIntegration
  export ivt_of_column, VerticalColumnData, IVTResult

  """
  Struct for one vertical column of atmospheric data. Containing the vertical column of specific humidity (hus), eastward_wind (ua), northward wind (va) and the pressure levels the data is available at. Additionally the surface pressure at the given grid point is required as a boundary
  """
  struct VerticalColumnData
    specific_humidity::Vector{Union{Float32, Missing}}
    eastward_wind::Vector{Union{Float32, Missing}}
    northward_wind::Vector{Union{Float32, Missing}}
    pressure_levels::Vector{Float64}
    surface_pressure::Float32
  end

  """
  Struct for the result of IVT calculation. Contains the eastward and northward integral spanning the vector.
  """
  struct IVTResult{T <: Real}
    eastward_integral::T
    northward_integral::T
  end

  function ph(i::Int, pressure_levels::AbstractArray{<:AbstractFloat, 1}, surface_pressure::Float64)
    if i == 1
      return surface_pressure
    elseif i == length(pressure_levels) + 1
      return 0
    else
      return (pressure_levels[i] + pressure_levels[i-1])/2
    end
    
  end
  
  function calculate_layer_values(i::Int, pressure_levels::AbstractArray{<:AbstractFloat, 1}, surface_pressure::Float64, specific_humidity::AbstractArray{<:AbstractFloat, 1}, northward_wind::AbstractArray{<:AbstractFloat, 1}, eastward_wind::AbstractArray{<:AbstractFloat, 1}, g::Float64 = 9.806)::Tuple{Float32, Float32}
      dp = ph(i, pressure_levels, surface_pressure) - ph(i + 1, pressure_levels, surface_pressure)
      # println("PH Values Index $i:\t$(ph(i, pressure_levels, surface_pressure))\t$(ph(i + 1, pressure_levels, surface_pressure)) = $dp")
      dm = dp/g
      # now here we strive from the description and multiply it also with the wind component
      
      va_dq = specific_humidity[i] * dm * northward_wind[i]
      ua_dq = specific_humidity[i] * dm * eastward_wind[i]
      
      # println("Multiplied Vals Index $i:\t$ua_dq")
      
      return ua_dq, va_dq
    
  end
  
  function ivt_of_column(column_data::VerticalColumnData)::IVTResult{Float32}

    nmax = size(column_data.pressure_levels, 1) + 1
    
    sum_va_hus = 0.
    sum_ua_hus = 0.
  
  
    for i in 1:nmax-1
      (ua_layer_value, va_layer_value) = calculate_layer_values(i, column_data.pressure_levels, column_data.surface_pressure, column_data.specific_humidity, column_data.northward_wind, column_data.eastward_wind)
      sum_va_hus += va_layer_value
      sum_ua_hus += ua_layer_value
    end
    
    return IVTResult(sum_ua_hus, sum_va_hus)
  end

  # function ivt_of_column(ps::Float32, plevs::Vector{Float64}, hus_data::Vector{Float32}, ua_data::Vector{Float32}, va_data::Vector{Float32})::IVTResult{Float32}
  #
  #   sum_va_hus = 0.
  #   sum_ua_hus = 0.
  # 
  # 
  #   for i ∈ eachindex(plevs)
  #     (ua_layer_value, va_layer_value) = calculate_layer_values(i, plevs, ps, hus_data, va_data, ua_data)
  #     sum_va_hus += va_layer_value
  #     sum_ua_hus += ua_layer_value
  #   end
  #   
  #   return IVTResult(sum_ua_hus, sum_va_hus)
  # end
  #
  # function ivt_of_column(ps::Float64, data::Matrix{Union{Float32, Missing}})
  #
  #   # corresponds in order with VerticalColumnData
  #   nmax = size(data[4, :], 1) + 1
  #   
  #   g = 9.806
  # 
  #   function ph(i::Int)::Float32
  #     if i == 1
  #       return 0 
  #     elseif i == nmax
  #       return ps
  #     else
  #     return (data[4, :] + data[i-1, :])/2
  #     end
  #   end
  # 
  #   function calculate_layer_values(i::Int)::Float32
  #     dp = ph(i + 1) - ph(i)
  #     dm = dp/g
  #     # now here we strive from the description and multiply it also with the wind component
  #     
  #     va_dq = data[1, i] * dm * data[3, i]
  #     ua_dq = data[1, i] * dm * data[2, i]
  #     
  #     return ua_dq, va_dq
  #   end
  # 
  #   sum_va_hus = 0.
  #   sum_ua_hus = 0.
  # 
  # 
  #   for i in 1:nmax-1
  #     (ua_layer_value, va_layer_value) = calculate_layer_values(i)
  #     sum_va_hus += va_layer_value
  #     sum_ua_hus += ua_layer_value
  #   end
  #   
  #   return sqrt(sum_ua_hus^2 + sum_va_hus^2)
  # end
  
function ivt_of_column(plevs::Vector{<: AbstractFloat}, hus_data::Vector{<: AbstractFloat}, ua_data::Vector{<: AbstractFloat}, va_data::Vector{<: AbstractFloat}, g::Float32 = Float32(9.806))::IVTResult{<: AbstractFloat}
    
    if plevs[1] < plevs[end]
      northward_integral = NumericalIntegration.integrate(plevs, hus_data .* va_data)/g
      eastward_integral = NumericalIntegration.integrate(plevs, hus_data .* ua_data)/g
    else
      northward_integral = NumericalIntegration.integrate(reverse(plevs), reverse(hus_data .* va_data))/g
      eastward_integral = NumericalIntegration.integrate(reverse(plevs), reverse(hus_data .* ua_data))/g
    end 
    return IVTResult(eastward_integral, northward_integral)
  end 


  function ivt_of_column_old(plevs::Vector{<: AbstractFloat}, hus_data::Vector{<: AbstractFloat}, ua_data::Vector{<: AbstractFloat}, va_data::Vector{<: AbstractFloat})::IVTResult{<: AbstractFloat}
    
    pressure_levels = view(plevs, 2:length(plevs))
    
    ps = plevs[1]
    
    sum_va_hus = 0.
    sum_ua_hus = 0.


    for i ∈ eachindex(pressure_levels)
      (ua_layer_value, va_layer_value) = calculate_layer_values(i, pressure_levels, ps, hus_data, va_data, ua_data)
      sum_va_hus += va_layer_value
      sum_ua_hus += ua_layer_value
    end
    
    return IVTResult(sum_ua_hus, sum_va_hus)
  end
end


