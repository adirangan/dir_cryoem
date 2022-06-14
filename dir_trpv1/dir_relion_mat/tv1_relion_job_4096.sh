/mnt/home/rangan/relion/build/bin/relion_refine \
--o job_4096/run \
--sgd_ini_iter 50 \
--sgd_inbetween_iter 200 \
--sgd_fin_iter 50 \
--sgd_write_iter 20 \
--sgd_ini_resol 35 \
--sgd_fin_resol 15 \
--sgd_ini_subset 100 \
--sgd_fin_subset 500 \
--sgd \
--denovo_3dref \
--i tv1_relion_job_4096.star \
--ctf \
--K 1 \
--sym C1 \
--flatten_solvent \
--zero_mask \
--dont_combine_weights_via_disc \
--pool 3 \
--pad 1 \
--skip_gridding \
--particle_diameter 400 \
--oversampling 1 \
--healpix_order 1 \
--offset_range 8 \
--offset_step 2 \
--j 6 \
--pipeline_control job_4096 \
;

