#source /opt/intel/mkl/bin/mklvars.sh intel64
#export LDFLAGS="-L/usr/local/opt/openblas/lib"
#export CPPFLAGS="-I/usr/local/opt/openblas/include"

# export LD_LIBRARY_PATH=/usr/local/lib/
# export LD_LIBRARY_PATH=/usr/local/opt/lapack/lib
# export LD_LIBRARY_PATH=/usr/local/opt/openblas/lib

# For compilers to find openblas you may need to set:
#  export LDFLAGS="-L/usr/local/opt/openblas/lib"
#  export CPPFLAGS="-I/usr/local/opt/openblas/include"

# For compilers to find lapack you may need to set:
#    export LDFLAGS="-L/usr/local/opt/lapack/lib"
#    export CPPFLAGS="-I/usr/local/opt/lapack/include"

make -f mymakefile
rm *.o
rm *.mod
