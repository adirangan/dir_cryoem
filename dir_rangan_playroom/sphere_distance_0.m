function distance = sphere_distance_0(polar_a_0,azimu_b_0,polar_a_1,azimu_b_1);
k_0_0 = sin(polar_a_0).*cos(azimu_b_0);
k_1_0 = sin(polar_a_0).*sin(azimu_b_0);
k_2_0 = cos(polar_a_0);
k_0_1 = sin(polar_a_1).*cos(azimu_b_1);
k_1_1 = sin(polar_a_1).*sin(azimu_b_1);
k_2_1 = cos(polar_a_1);
distance = sqrt( (k_0_0-k_0_1).^2 + (k_1_0-k_1_1).^2 + (k_2_0-k_2_1).^2 );
