#Gravity model of commuting in Kenya

#Calculate interaction force for gravity model --- assuming Urban pop
function TupleDist(x1::Tuple{Float64,Float64},x2::Tuple{Float64,Float64})
    return sqrt((x1[1]-x2[1])^2 + (x1[2]-x2[2])^2 )
end

F = zeros(n,n)
for i = 1:n
    for j = 1:n
        if i != j
        F[i,j] = KenyaTbl[:Urban][i]*KenyaTbl[:Urban][j]/(TupleDist(KenyaTbl[:Loc][i],KenyaTbl[:Loc][j])^2)
    end
    end
end
#Normalise the columns
for j = 1:n
    F[:,j] = F[:,j]/sum(F[:,j])
end
