# Read in the data for the model


function get_Kenyadata(filename)
    KenyaTbl = readtable(filename)
    KenyaTbl = KenyaTbl[1:(end-1),:]
    n, = size(KenyaTbl)
    #Equator conversion 1 deg lat = 110.57km, 1 deg long = 111.32 km
    #Parse the string into an x,y tuple (in kms)
    KenyaTbl[:Loc] = [(0.,0.) for i = 1:n]
    KenyaTbl[:UrbanIndex] = ones(n)
    KenyaTbl[:Phase] = zeros(n)
    for i = 1:n
        x = join(split(KenyaTbl[:Location_1][i],"("))
        x = join(split(x,")"))
        y = split(x,",")
        z = parse(Float64,y[2])*111.32,parse(Float64,y[1])*110.57
        KenyaTbl[:Loc][i] = z
        if ~ismissing(KenyaTbl[:Rural][i])
            KenyaTbl[:UrbanIndex][i] = KenyaTbl[:Urban][i]/(KenyaTbl[:Urban][i]+KenyaTbl[:Rural][i])
        end
        KenyaTbl[:Phase] = 30.
    end
    return KenyaTbl,n
end
