#!/bin/bash
# Run computations for all regular triangles (at most one singular corner)

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

echo "nohup $prog $args -i 0 &> $dir/3_4_1_3_1_2 &"
nohup $prog $args       -i 0 &> $dir/3_4_1_3_1_2 &
echo "nohup $prog $args -i 1 &> $dir/2_3_1_3_1_2 &"
nohup $prog $args       -i 1 &> $dir/2_3_1_3_1_2 &
echo "nohup $prog $args -i 2 &> $dir/2_3_1_4_1_2 &"
nohup $prog $args       -i 2 &> $dir/2_3_1_4_1_2 &
echo "nohup $prog $args -i 3 &> $dir/2_3_1_3_1_3 &"
nohup $prog $args       -i 3 &> $dir/2_3_1_3_1_3 &
echo "nohup $prog $args -i 4 &> $dir/3_4_1_4_1_3 &"
nohup $prog $args       -i 4 &> $dir/3_4_1_4_1_3 &
echo "nohup $prog $args -i 5 &> $dir/2_3_1_4_1_4 &"
nohup $prog $args       -i 5 &> $dir/2_3_1_4_1_4 &
