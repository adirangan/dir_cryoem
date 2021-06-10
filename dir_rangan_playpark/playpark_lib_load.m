system('make -f playpark_OptiPlex.make clean; make -f playpark_OptiPlex.make playpark_lib;');
if (libisloaded('playpark_lib')); unloadlibrary('playpark_lib'); end;
loadlibrary('playpark_lib','playpark_lib.h','includepath','/home/rangan/dir_finufft/include/');
libfunctions('playpark_lib');
