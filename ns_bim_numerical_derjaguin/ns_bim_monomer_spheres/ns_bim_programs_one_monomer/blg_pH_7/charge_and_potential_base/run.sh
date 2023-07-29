export OMP_NUM_THREADS=30
./run_DH > myjob.log 2>&1& # Run the given command line in the background.
pid=$!
echo $pid > myjob_pid.log
