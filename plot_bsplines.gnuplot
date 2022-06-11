reset

unset key; 

max=200

plot for [i=2:max:2] 'bsplines.dat' u 1:(column(i)) w l lc rgb 'black' t ''.i, \
     for [i=3:max:2] '' u 1:(column(i)) w l lc rgb 'red' t ''.i


