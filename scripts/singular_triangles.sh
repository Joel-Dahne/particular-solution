#!/bin/bash
# Run computations for all singular triangles (two or more singular corners)

# Parameters
beg="4" # Starting value for N
step="2" # Steps in N
end="16" # End value for N
tol="1e-5" # Relative tolerance to use in the minimization
prec="64" # Precision to use
type="1" # Type of output

args="-N_beg $beg -N_step $step -N_end $end -tol $tol -prec $prec -o $type -time"
dir="data/enclosures_"$beg"_"$step"_"$end"_"$tol"_"$prec
prog="build/examples/triangles"

echo $args

mkdir -p $dir

echo "beg = "$beg   > $dir/parameters
echo "step = "$step >> $dir/parameters
echo "end = "$end   >> $dir/parameters
echo "tol = "$tol   >> $dir/parameters
echo "prec = "$prec >> $dir/parameters
echo "type = "$type >> $dir/parameters

echo "nohup $prog $args -i 6 &> $dir/2_3_3_4_3_4 &"
nohup $prog $args       -i 6 &> $dir/2_3_3_4_3_4 &
echo "nohup $prog $args -i 7 &> $dir/2_3_2_3_2_3 &"
nohup $prog $args       -i 7 &> $dir/2_3_2_3_2_3 &
echo "nohup $prog $args -i 8 &> $dir/1_2_2_3_3_4 &"
nohup $prog $args       -i 8 &> $dir/1_2_2_3_3_4 &
echo "nohup $prog $args -i 9 &> $dir/1_2_2_3_2_3 &"
nohup $prog $args       -i 9 &> $dir/1_2_2_3_2_3 &
