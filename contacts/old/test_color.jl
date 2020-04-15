using Plots, Images#, ImageView
img = load("./contacts/kenya.jpg")
size(img)[1]


plot(img,aspect_ratio=1,axis=false,xlims = (0,size(img)[2]), ylims=(0,size(img)[1]))
#size_factor=.8
#plot(img,aspect_ratio=1,axis=false,xlims = (0,size(img)[1]*size_factor), ylims=(0,size(img)[2]*size_factor))
annotate!(1.,1.,text("test",5))
circle!(10.,10.)

#gui=imshow(img)


using Luxor
img = readpng("C:/Users/rabia/Documents/GitHub/KenyaCoV/contacts/kenya.png")
w = img.width
h = img.height
placeimage(img, -w/2, -h/2, .5)
placeimage(img, 0, 0; centered=false)



using Plots, Images, ImageDraw
img = load("./contacts/kenya.jpg")
draw!(img, CirclePointRadius(Point(floor(628/2),floor(735/2)),500))

c=ColorGradient([:red,:yellow,:blue])
CList = reshape( range(colorant"red", stop=colorant"blue",length=100), 1, 100 )
CList[90]
