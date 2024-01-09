
set fit quiet
set fit logfile "/dev/null"
set key right bottom

ERRFILE = "error.dat"

set logscale xy
set format x "%2.0t.10^{%L}"
set format y "%2.0t.10^{%L}"

set title "Convergence Analysis"

# Function for linear fit in log-log scale
f(x) = a*x + b

# Fit for density
fit f(x) ERRFILE u (log($1)):(log($2)) via a, b
set label 1 sprintf("Density slope: %.2f", a) at graph 0.02, 0.95
plot exp(f(log(x))) lc rgb 'gray' dt 4 t sprintf("Density fit")
replot ERRFILE u 1:2 w lp t "Density"

# Fit for velocity_x
fit f(x) ERRFILE u (log($1)):(log($3)) via a, b
set label 2 sprintf("Velocity_x slope: %.2f", a) at graph 0.02, 0.90
replot exp(f(log(x))) lc rgb 'blue' dt 4 t sprintf("Velocity_x fit")
replot ERRFILE u 1:3 w lp t "Velocity_x"

# Fit for velocity_y
fit f(x) ERRFILE u (log($1)):(log($4)) via a, b
set label 3 sprintf("Velocity_y slope: %.2f", a) at graph 0.02, 0.85
replot exp(f(log(x))) lc rgb 'green' dt 4 t sprintf("Velocity_y fit")
replot ERRFILE u 1:4 w lp t "Velocity_y"

# Fit for energy
fit f(x) ERRFILE u (log($1)):(log($5)) via a, b
set label 4 sprintf("Energy slope: %.2f", a) at graph 0.02, 0.80
replot exp(f(log(x))) lc rgb 'red' dt 4 t sprintf("Energy fit")
replot ERRFILE u 1:5 w lp t "Energy"

pause -1 
