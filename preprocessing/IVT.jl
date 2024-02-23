module IVT

  export ivt_of_column, VerticalColumnData

  """
  Struct for one vertical column of atmospheric data. Containing the vertical column of specific humidity (hus), eastward_wind (ua), northward wind (va) and the pressure levels the data is available at. Additionally the surface pressure at the given grid point is required as a boundary
  """
  struct VerticalColumnData
    specific_humidity::Vector{Union{Float32, Missing}}
    eastward_wind::Vector{Union{Float32, Missing}}
    southward_wind::Vector{Union{Float32, Missing}}
    pressure_levels::Vector{Float64}
    surface_pressure::Float32
  end
  
  function ivt_of_column(column_data::VerticalColumnData)::Union{Float64, Missing}

    nmax = size(column_data.pressure_levels, 1) + 1
    
    g = 9.806

    ps = column_data.surface_pressure
  
    function ph(i::Int)
      if i == 1
        return 0 
      elseif i == nmax
        return ps
      else
      return (column_data.pressure_levels[i] + column_data.pressure_levels[i-1])/2
      end
    end
  
    function calculate_layer_values(i::Int)
      dp = ph(i + 1) - ph(i)
      dm = dp/g
      # now here we strive from the description and multiply it also with the wind component
      
      va_dq = column_data.specific_humidity[i] * dm * column_data.southward_wind[i]
      ua_dq = column_data.specific_humidity[i] * dm * column_data.eastward_wind[i]
      
      return ua_dq, va_dq
    end
  
    sum_va_hus = 0.
    sum_ua_hus = 0.
  
  
    for i in 1:nmax-1
      (ua_layer_value, va_layer_value) = calculate_layer_values(i)
      sum_va_hus += va_layer_value
      sum_ua_hus += ua_layer_value
    end
    
    return sqrt(sum_ua_hus^2 + sum_va_hus^2)
  end
  
end


