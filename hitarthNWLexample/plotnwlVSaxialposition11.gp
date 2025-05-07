# plot nwl vs axial position in mirror geometry
#
clear;reset
#
set term wxt enhanced font "arial,20"
#
set title 'Realta Anvil System'
set xlabel "Axial position in Mirror (cm) \n"
set ylabel "NWL (MW/m2)"
#
set key top right box
#
set xrange [-600:600]
set mxtics 4
#
set yrange [0:.4]
#set ytics .5
set mytics 5
#
# axial boundaries and regions
# -597.5
# lower cyl
# -325
# lower cone
# -127.97
# central cyl
# 127.97
# upper cone
# 325
# upper cyl
# 597.5
#
# draw vertical lines
set arrow from -597.5,0 to -597.5,0.4 nohead front lt -0 lw 1 # lower end region
set arrow from -325,0 to -325,0.4 nohead front lt -0 lw 1 # lower end region
set label font "arial,12" "Lower End Region (cylinder)" at -550,.33 # label for lower end region
set label font "arial,12" "Lower Conical Region \n (cone)" at -300,.35 # label for lower cone region
#
set arrow from -127.9675,0 to -127.9675,0.4 nohead front lt -0 lw 1 # mirror central region
set arrow from 127.9675,0 to 127.9675,0.36 nohead front lt -0 lw 1 # mirror central region
set label font "arial,12" "Central Region (cylinder)" at -100,.33 # label for central region
#
set label font "arial,12" "Upper Conical Region \n (cone)" at 150,.33 # label for Upper cone region
set arrow from 325,0 to 325,0.36 nohead front lt -0 lw 1 # upper end region
set arrow from 597.5,0 to 597.5,0.4 nohead front lt -0 lw 1 # lower end region
set label font "arial,12" "Upper End Region (cylinder)" at 350,.33 # label for upper end region
#
set label font "arial,10" "(Number of Src. Particles=1e6)" at 350,.355 #
#
#
# plot the data (src_strength=1.6386e18 n/s, surfsrcSample=1e6, scoring region height= 5 cm or 0.05 m) for the Anvil 2-D source
#                                                             Samples     convert MeV to MW
plot 'testCountsBinnedRadius.txt' using ($1+2.5):(($2*1.6386e18/1e6*14.1*1.602e-13/1e6)/(0.05*2*3.14159*$4/100.)) title 'Anvil 2-D Src Front of FWarmor' with linespoints lw 3 pt 7 lc rgbcolor "black"
#
#
# plot the data (src_strength=1.6386e18 n/s, surfsrcSample=1e6, scoring region height= 5 cm or 0.05 m) for Anvil 2-D source and for a uniform cylindrical source
plot 'testCountsBinnedRadius.txt' using ($1+2.5):(($2*1.6386e18/1e6*14.1*1.602e-13/1e6)/(0.05*2*3.14159*$4/100.)) title 'Anvil 2-D Source' with linespoints lw 3 pt 7 lc rgbcolor "red",'../testCountsBinnedRadius.txt' using ($1+2.5):(($2*1.6386e18/1e6*14.1*1.602e-13/1e6)/(0.05*2*3.14159*$4/100.)) title 'Cyl. Source (r=0-20, z=-450 to 450cm)' with linespoints lw 3 pt 9 lc rgbcolor "black"
