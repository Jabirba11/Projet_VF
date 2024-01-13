
set fit quiet
set fit logfile "/dev/null"
set key right bottom

set xlabel "delta_x" font ",12"
set ylabel "Erreur" font ",12"

ERRFILE = "WENO3_500x3.dat"

set logscale xy
set format x "%2.0t.10^{%L}"
set format y "%2.0t.10^{%L}"

set title "Analyse de convergence du schéma de WENO3" font ",14"

# Function for linear fit in log-log scale
f(x) = a*x + b

# Fit for density
fit f(x) ERRFILE u (log($1)):(log($2)) via a, b
set label 1 sprintf("Densite pente: %.2f", a) at graph 0.02, 0.95
plot exp(f(log(x))) lc rgb 'gray' dt 4 t sprintf("Densit=e fit")
replot ERRFILE u 1:2 w lp t "Densite"

# Fit for velocity_x
fit f(x) ERRFILE u (log($1)):(log($3)) via a, b
set label 2 sprintf("Vitesse_x pente: %.2f", a) at graph 0.02, 0.90
replot exp(f(log(x))) lc rgb 'blue' dt 4 t sprintf("Vitesse_x fit")
replot ERRFILE u 1:3 w lp t "Vitesse_x"

# Fit for velocity_y
fit f(x) ERRFILE u (log($1)):(log($4)) via a, b
set label 3 sprintf("Vitesse_y pente: %.2f", a) at graph 0.02, 0.85
replot exp(f(log(x))) lc rgb 'green' dt 4 t sprintf("Vitesse_y fit")
replot ERRFILE u 1:4 w lp t "Vitesse_y"

# Fit for energy
fit f(x) ERRFILE u (log($1)):(log($5)) via a, b
set label 4 sprintf("Energie pente: %.2f", a) at graph 0.02, 0.80
replot exp(f(log(x))) lc rgb 'red' dt 4 t sprintf("Energie fit")
replot ERRFILE u 1:5 w lp t "Energie"

pause -1 
