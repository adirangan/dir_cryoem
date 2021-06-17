%make -f test_ver18.make clean;
cd /data/rangan/dir_cryoem/dir_rangan_playground/;
make -f test_ver18.make;
./test_ver18.out 2 48 1 0 0 0 0 1 0 0
./test_ver18.out 2 48 1 0 0 0 0 1 0 1
./test_ver18.out 2 48 1 0 0 0 0 1 3 0
./test_ver18.out 2 48 1 0 0 0 0 1 3 1
./test_ver18.out 2 48 1 0 0 0 0 1 5 0
./test_ver18.out 2 48 1 0 0 0 0 1 5 1

./test_ver18.out 2 64 1 0 0 0 0 1 3 1
