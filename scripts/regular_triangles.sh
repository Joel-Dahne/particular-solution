#!/bin/bash
# Run computations for all regular triangles (at most one singular corner)

# Parameters
beg="4" # Starting value for N
step="4" # Steps in N
end="8" # End value for N
tol="1e-10" # Tolerance to use in the minimization
wid="1e-4" # Width of the interval around the initial guess
prec="53" # Precision to use
type="1" # Type of output

args="-b $beg -s $step -e $end -t $tol -w $wid -p $prec -o $type"
dir="data/nus_"$beg"_"$step"_"$end"_"$tol"_"$wid"_"$prec
prog="build/particular-solution"

echo $args

mkdir -p $dir

echo "beg = "$beg   > $dir/parameters
echo "step = "$step >> $dir/parameters
echo "end = "$end   >> $dir/parameters
echo "tol = "$tol   >> $dir/parameters
echo "wid = "$wid   >> $dir/parameters
echo "prec = "$prec >> $dir/parameters
echo "type = "$type >> $dir/parameters

echo "nohup $prog $args -n 3.056691018 -- 3 4 1 3 1 2 &> $dir/3_4_1_3_1_2 &"
nohup $prog $args -n 3.056691018 -- 3 4 1 3 1 2 &> $dir/3_4_1_3_1_2 &
echo "nohup $prog $args -n 3.240902298 -- 2 3 1 3 1 2 &> $dir/2_3_1_3_1_2 &"
nohup $prog $args -n 3.240902298 -- 2 3 1 3 1 2 &> $dir/2_3_1_3_1_2 &
echo "nohup $prog $args -n 4.063109028 -- 2 3 1 4 1 2 &> $dir/2_3_1_4_1_2 &"
nohup $prog $args -n 4.063109028 -- 2 3 1 4 1 2 &> $dir/2_3_1_4_1_2 &
echo "nohup $prog $args -h -i -n 4.143210850 -- 2 3 1 3 1 3 &> $dir/2_3_1_3_1_3 &"
nohup $prog $args -h -i -n 4.143210850 -- 2 3 1 3 1 3 &> $dir/2_3_1_3_1_3 &
echo "nohup $prog $args -n 4.470604591 -- 3 4 1 4 1 3 &> $dir/3_4_1_4_1_3 &"
nohup $prog $args -n 4.470604591 -- 3 4 1 4 1 3 &> $dir/3_4_1_4_1_3 &
echo "nohup $prog $args -h -i -n 6.525663100 -- 2 3 1 4 1 4 &> $dir/2_3_1_4_1_4 &"
nohup $prog $args -h -i -n 6.525663100 -- 2 3 1 4 1 4 &> $dir/2_3_1_4_1_4 &
