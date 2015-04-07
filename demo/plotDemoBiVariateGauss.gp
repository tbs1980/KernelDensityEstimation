set terminal png
set output "demoBiVariateGauss.png"
set datafile separator ","
set xyplane at 0
set palette rgbformulae 33,13,10
splot "plot-bi-variate-Gauss.dat" w pm3d
