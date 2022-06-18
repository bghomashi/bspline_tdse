reset

unset key; 

max=300

plot [0:50] \
     for [i=2:max:2] 'basis.txt' u 1:(column(i)) w l lc rgb 'black' t ''.i, \
     for [i=3:max:2] '' u 1:(column(i)) w l lc rgb 'red' t ''.i


