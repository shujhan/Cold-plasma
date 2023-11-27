    for OMP_NUM_THREADS in 1 2 3 4 5 6 7 8 12 16 18
	do 
	    export OMP_NUM_THREADS=$OMP_NUM_THREADS
	    echo "For $OMP_NUM_THREADS threads: "
	    ./omp.exe
	done
