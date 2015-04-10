set terminal png
set output "demoPowSpecSamples2d.png"
set datafile separator ","
set xyplane at 0
set palette rgbformulae 33,13,10
splot "plot-powSpecSamples2D.dat" w pm3d
