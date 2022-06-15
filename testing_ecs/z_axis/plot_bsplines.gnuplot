reset

unset key; 

min=0
max=500

plot for [i=2:(2*max):2] 'basis.txt' u 1:(column(i)) w l lc rgb 'black' t ''.i, \
     for [i=3:(2*max):2] '' u 1:(column(i)) w l lc rgb 'red' t ''.i


