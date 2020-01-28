#Gravity model of commuting in Kenya

#Calculate interaction force for gravity model --- assuming Urban pop
function TupleDist(x1::Tuple{Float64,Float64},x2::Tuple{Float64,Float64})
    return sqrt((x1[1]-x2[1])^2 + (x1[2]-x2[2])^2 )
end

transport_matrix = zeros(n,n)
for i = 1:n
    for j = 1:n
        if i != j
        transport_matrix[i,j] = KenyaTbl[:Urban][i]*KenyaTbl[:Urban][j]/(TupleDist(KenyaTbl[:Loc][i],KenyaTbl[:Loc][j])^2)
    end
    end
end
#Normalise the columns
for j = 1:n
    transport_matrix[:,j] = transport_matrix[:,j]/sum(transport_matrix[:,j])
end
