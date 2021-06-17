mv ti8_excerpt_define_Zstore_0a.f ti8_excerpt_define_Zstore_part_a.f
mv ti8_excerpt_define_Zstore_0b.f ti8_excerpt_define_Zstore_part_b.f
mv ti8_excerpt_define_Zstore_0c.f ti8_excerpt_define_Zstore_part_c.f
mv ti8_excerpt_define_Zstore_0x.f ti8_excerpt_define_Zstore_part_x.f
mv ti8_excerpt_define_Zstore_0d.f ti8_excerpt_define_Zstore_part_d.f
mv ti8_excerpt_define_Zstore_0e.f ti8_excerpt_define_Zstore_part_e.f
mv ti8_excerpt_define_Zstore_1.f ti8_excerpt_define_Zstore_full.f
mv ti8l_excerpt_define_Zstore_1.f ti8l_excerpt_define_Zstore_full.f

find . -type f -name "*.f" -exec sed -i 's/ti8_excerpt_define_Zstore_0a.f/ti8_excerpt_define_Zstore_part_a.f/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_excerpt_define_Zstore_0b.f/ti8_excerpt_define_Zstore_part_b.f/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_excerpt_define_Zstore_0c.f/ti8_excerpt_define_Zstore_part_c.f/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_excerpt_define_Zstore_0x.f/ti8_excerpt_define_Zstore_part_x.f/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_excerpt_define_Zstore_0d.f/ti8_excerpt_define_Zstore_part_d.f/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_excerpt_define_Zstore_0e.f/ti8_excerpt_define_Zstore_part_e.f/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_excerpt_define_Zstore_1.f/ti8_excerpt_define_Zstore_full.f/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8l_excerpt_define_Zstore_1.f/ti8l_excerpt_define_Zstore_full.f/g' {} \;

mv ti8_Zstore_3a.f ti8_build_Zstore_part_a.f
mv ti8_Zstore_3b.f ti8_build_Zstore_part_b.f
mv ti8_Zstore_3c.f ti8_build_Zstore_part_c.f
mv ti8_Zstore_3d.f ti8_build_Zstore_part_d.f
mv ti8_Zstore_3e.f ti8_build_Zstore_part_e.f
mv ti8_Zstore_3x.f ti8_build_Zstore_part_x.f
mv ti8_Zstore_3y.f ti8_build_Zstore_part_y.f
mv ti8_Zstore_4.f ti8_build_Zstore_full_0.f

find . -type f -name "*.f" -exec sed -i 's/ti8_Zstore_3a/ti8_build_Zstore_part_a/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_Zstore_3b/ti8_build_Zstore_part_b/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_Zstore_3c/ti8_build_Zstore_part_c/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_Zstore_3d/ti8_build_Zstore_part_d/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_Zstore_3e/ti8_build_Zstore_part_e/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_Zstore_3x/ti8_build_Zstore_part_x/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_Zstore_3y/ti8_build_Zstore_part_y/g' {} \;
find . -type f -name "*.f" -exec sed -i 's/ti8_Zstore_4/ti8_build_Zstore_full_0/g' {} \;

find . -type f -name "*.make" -exec sed -i 's/ti8_Zstore_3a/ti8_build_Zstore_part_a/g' {} \;
find . -type f -name "*.make" -exec sed -i 's/ti8_Zstore_3b/ti8_build_Zstore_part_b/g' {} \;
find . -type f -name "*.make" -exec sed -i 's/ti8_Zstore_3c/ti8_build_Zstore_part_c/g' {} \;
find . -type f -name "*.make" -exec sed -i 's/ti8_Zstore_3d/ti8_build_Zstore_part_d/g' {} \;
find . -type f -name "*.make" -exec sed -i 's/ti8_Zstore_3e/ti8_build_Zstore_part_e/g' {} \;
find . -type f -name "*.make" -exec sed -i 's/ti8_Zstore_3x/ti8_build_Zstore_part_x/g' {} \;
find . -type f -name "*.make" -exec sed -i 's/ti8_Zstore_3y/ti8_build_Zstore_part_y/g' {} \;
find . -type f -name "*.make" -exec sed -i 's/ti8_Zstore_4/ti8_build_Zstore_full_0/g' {} \;

