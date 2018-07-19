#!/bin/bash
# Run computations for all singular triangles (two or more singular corners)

# Parameters
beg="4" # Starting value for N
step="2" # Steps in N
end="16" # End value for N
tol="1e-5" # Relative tolerance to use in the minimization
wid="1e-4" # Width of the interval around the initial guess
prec="53" # Precision to use
type="4" # Type of output

args="-b $beg -s $step -e $end -t $tol -w $wid -p $prec -o $type"
dir="data/enclsures_"$beg"_"$step"_"$end"_"$tol"_"$wid"_"$prec
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

echo "nohup $prog $args -n 1.624084509 -- 2 3 3 4 3 4 &> $dir/2_3_3_4_3_4 &"
nohup $prog $args -n 1.624084509 -- 2 3 3 4 3 4 &> $dir/2_3_3_4_3_4 &
echo "nohup $prog $args -n 1.825757081 -- 2 3 2 3 2 3 &> $dir/2_3_2_3_2_3 &"
nohup $prog $args -n 1.825757081 -- 2 3 2 3 2 3 &> $dir/2_3_2_3_2_3 &
echo "nohup $prog $args -n 2.047890892 -- 2 3 3 4 1 2 &> $dir/2_3_3_4_1_2 &"
nohup $prog $args -n 2.047890892 -- 2 3 3 4 1 2&> $dir/2_3_3_4_1_2 &
echo "nohup $prog $args -n 2.150869291 -- 2 3 2 3 1 2 &> $dir/2_3_2_3_1_2 &"
nohup $prog $args -n 2.150869291 -- 2 3 2 3 1 2 &> $dir/2_3_2_3_1_2 &
