module IVT

  export ivt_of_column, VerticalColumnData

  struct VerticalColumnData
  """Struct for one vertical column of atmospheric data. Containing the vertical column of specific humidity (hus), """
    specific_humidity::Vector{Union{Float32, Missing}}
    eastward_wind::Vector{Union{Float32, Missing}}
    southward_wind::Vector{Union{Float32, Missing}}
    surface_pressure::Float32
  end
  
  function ivt_of_column(column_data::VerticalColumnData, pressure_levels::Vector{Float64})::Union{Float32, Missing}

    nmax = size(pressure_levels, 1) + 1
    
    g = 9.806

    ps = column_data.surface_pressure
  
    function ph(i::Int)::Float32
      if i == 1
        return 0 
      elseif i == nmax
        return ps
      else
      return ((pressure_levels[i] * ps) + (pressure_levels[i-1] * ps))/2
      end
    end
  
    function calculate_layer_values(i::Int)::Tuple{Union{Float32, Missing}, Union{Float32, Missing}}
      dp = ph(i + 1) - ph(i)
      dm = dp/g
      # now here we strive from the description and multiply it also with the wind component
      
      va_dq = column_data.specific_humidity[i] * dm * column_data.southward_wind[i]
      ua_dq = column_data.specific_humidity[i] * dm * column_data.eastward_wind[i]
  
      return va_dq, ua_dq
    end
  
    sum_va_hus = 0
    sum_ua_hus = 0
  
  
    for i in 1:nmax-1
      (va_layer_value, ua_layer_value) = calculate_layer_values(i)
      sum_va_hus += va_layer_value
      sum_ua_hus += ua_layer_value
    end
    
    return sqrt(sum_ua_hus^2 + sum_va_hus^2)
  end
  
end


