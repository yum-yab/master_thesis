module IVT

  export ivt_of_column, PressureLevelData

  struct PressureLevelData
    hus::Union{Float32, Missing}
    ua::Union{Float32, Missing}
    va::Union{Float32, Missing}
    p::Float64
  end
  
  function ivt_of_column(column_data::Vector{PressureLevelData})::Float32

    nmax = size(column_data, 1) + 1
    
    ps = maximum(map(pld -> pld.p, column_data))
    g = 9.806
  
    function ph(i::Int)::Float32
      if i == 0
        return 0 
      elseif i == nmax
        return ps
      else
        return (column_data[i].p + column_data[i-1])/2
      end
    end
  
    function calculate_layer_values(i::Int)::Tuple{Float32, Float32}
      dp = ph(i + 1) - ph(i)
      dm = dp/g
      # now here we strive from the description and multiply it also with the wind component
      
      va_dq = column_data[i].hus * dm * column_data[i].va 
      ua_dq = column_data[i].hus * dm * column_data[i].ua 
  
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


