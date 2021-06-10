c     CTF parameters

      real *8 CTF_Spherical_Aberration !spherical aberation of the lens ;
      real *8, allocatable :: CTF_Spherical_Aberration_(:)
      real *8 CTF_Voltage_kV,CTF_Voltage_1V !voltage ;
      real *8,allocatable :: CTF_Voltage_kV_(:)
      real *8 CTF_Amplitude_Contrast !amplitude contrast ;
      real *8, allocatable :: CTF_Amplitude_Contrast_(:)
      real *8 CTF_Magnification !magnification of microscope ;
      real *8, allocatable :: CTF_Magnification_(:)
      real *8 CTF_Detector_Pixel_Size !pixel size of scanner in microns ;
      real *8 CTF_B_factor
      real *8, allocatable :: CTF_Detector_Pixel_Size_(:)
      real *8 CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle !defocus (in Angstroms) and angle of astigmatism ;
      real *8, allocatable :: CTF_Defocus_U_(:)
      real *8, allocatable :: CTF_Defocus_V_(:)
      real *8, allocatable :: CTF_Defocus_Angle_(:)
      real *8 tmp_w1,tmp_w2,tmp_k_c_1,tmp_k_c_2,tmp_ctf_value,tmp_theta
      real *8 CTF_lambda,CTF_Object_Pixel_Size,CTF_lambda_per_box
      real *8 CTF_Defocus_max, CTF_Defocus_min
      real *8 CTF_Astigmatism_max, CTF_Astigmatism_min
      real *8 CTF_Astigmatism_value

      allocate(CTF_Voltage_kV_(0:1+n_ctf-1))
      call cs1_r8(n_ctf,CTF_Voltage_kV_);
      allocate(CTF_Defocus_U_(0:1+n_ctf-1))
      call cs1_r8(n_ctf,CTF_Defocus_U_);
      allocate(CTF_Defocus_V_(0:1+n_ctf-1))
      call cs1_r8(n_ctf,CTF_Defocus_V_);
      allocate(CTF_Defocus_Angle_(0:1+n_ctf-1))
      call cs1_r8(n_ctf,CTF_Defocus_Angle_);
      allocate(CTF_Spherical_Aberration_(0:1+n_ctf-1))
      call cs1_r8(n_ctf,CTF_Spherical_Aberration_);
      allocate(CTF_Detector_Pixel_Size_(0:1+n_ctf-1))
      call cs1_r8(n_ctf,CTF_Detector_Pixel_Size_);
      allocate(CTF_Magnification_(0:1+n_ctf-1))
      call cs1_r8(n_ctf,CTF_Magnification_);
      allocate(CTF_Amplitude_Contrast_(0:1+n_ctf-1))
      call cs1_r8(n_ctf,CTF_Amplitude_Contrast_);
c$$$  the CTF parameters are taken from the rib80s dataset. ;
      CTF_Voltage_kV = 300.0d0
      CTF_Spherical_Aberration = 2.0d0
      CTF_Detector_Pixel_Size = 3.7687499999999998d0
      CTF_Amplitude_Contrast = 0.1d0
      CTF_B_factor = 50
      CTF_Magnification = 1.0d0 !not used yet ;
      CTF_Defocus_max = 26462.15039099d0 !maximum possible value of defocus ;
      CTF_Defocus_min = 8112.450195d0 !minimum possible value of defocus ;
      CTF_Astigmatism_max = 457.03125d0 !max astigmatism (CTF_Defocus_U-CTF_Defocus_V) ;
      CTF_Astigmatism_min = 10.509765d0 !min astigmatism ;

      do nctf=0,n_ctf-1
c$$$  generate a value between the min and max defocus ;
         CTF_Defocus_V = adi_rand_f(rseed)*(CTF_Defocus_max
     $        -CTF_Defocus_min)
         CTF_Defocus_V = CTF_Defocus_V + CTF_Defocus_min
c$$$  generate a value between max and min astigmatism ;
         CTF_Astigmatism_value = adi_rand_f(rseed)*(CTF_Astigmatism_max
     $        -CTF_Astigmatism_min)
         CTF_Astigmatism_value = CTF_Astigmatism_value +
     $        CTF_Astigmatism_min
c$$$  ensure that CTF_Defocus_U > CTF_Defocus_U ;
         CTF_Defocus_U = CTF_Defocus_V+CTF_Astigmatism_value
         CTF_Voltage_kV_(nctf) = CTF_Voltage_kV
         CTF_Defocus_U_(nctf) = CTF_Defocus_U
         CTF_Defocus_V_(nctf) = CTF_Defocus_V
         CTF_Defocus_Angle_(nctf) = CTF_Defocus_Angle
         CTF_Spherical_Aberration_(nctf) = CTF_Spherical_Aberration
         CTF_Detector_Pixel_Size_(nctf) = CTF_Detector_Pixel_Size
         CTF_Magnification_(nctf) = CTF_Magnification
         CTF_Amplitude_Contrast_(nctf) = CTF_Amplitude_Contrast
      enddo !do nctf=0,n_ctf-1

c-----------------------------------------------------------------------
c$$$ 5) generate ctf functions to associate with slices. ;
c-----------------------------------------------------------------------
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      do nctf=0,n_ctf-1
c$$$  spherical aberation of the lens in mm ;
         CTF_Spherical_Aberration = CTF_Spherical_Aberration_(nctf)
c$$$  convert into Angstroms ;
         CTF_Spherical_Aberration=CTF_Spherical_Aberration*(10.0d0
     $        **7.0d0)
c$$$  voltage in kVolts ;
         CTF_Voltage_kV = CTF_Voltage_kV_(nctf)
c$$$  convert into Volts ;
         CTF_Voltage_1V=CTF_Voltage_kV*1000.0 
c$$$  electron wavelength in Angstroms ;
         CTF_lambda = 12.2643247/dsqrt(CTF_Voltage_1V+CTF_Voltage_1V**2
     $        *0.978466d-6)
c$$$  defocus values (in Angstroms) ;
         CTF_Defocus_U = CTF_Defocus_U_(nctf)
         CTF_Defocus_V = CTF_Defocus_V_(nctf)
c$$$  angle of astigmatism ;
         CTF_Defocus_Angle = CTF_Defocus_Angle_(nctf)
c$$$  convert into radians ;
         CTF_Defocus_Angle = CTF_Defocus_Angle*pi/180.0d0
c$$$  CTF_Amplitude Contrast ;
         CTF_Amplitude_Contrast = CTF_Amplitude_Contrast_(nctf)
c$$$  weights for the amplitude and phase contrasts in CTF ;
         tmp_w1=dsqrt(1.0d0-CTF_Amplitude_Contrast**2)
         tmp_w2=CTF_Amplitude_Contrast
c$$$  pixel size of the scanner in physical space (not magnified) in Angstroms ;
c$$$   CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size/CTF_Magnification
      CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size_(nctf)
c$$$  n_x_c_max*CTF_Object_Pixel_Size is the box size in Angstroms ;
      CTF_lambda_per_box = CTF_lambda/(n_x_c_max*CTF_Object_Pixel_Size)
c$$$   call envelope_fxn(ngridr,xnodesr/pi,D,envelope)
      na=0
      do nk = 0,n_k_p_max-1
         do nw = 0,n_w_(nk)-1
            tmp_theta = (2.0d0*pi*nw)/n_w_(nk)
            tmp_k_c_1 = grid_k_p_(nk)*dcos(tmp_theta)
            tmp_k_c_2 = grid_k_p_(nk)*dsin(tmp_theta)
            call niko_ctf(CTF_Spherical_Aberration,CTF_lambda,tmp_w1
     $           ,tmp_w2,CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle
     $           ,CTF_lambda_per_box,tmp_k_c_1/pi,tmp_k_c_2/pi
     $           ,tmp_ctf_value)
            CTF_k_p__(na + nctf*ld_CTF) = -tmp_ctf_value !Note change in sign
            na = na+1
         enddo !do nw = 0,n_w_(nk)-1
      enddo !do nk = 0,n_k_p_max-1
      enddo !do nctf=0,n_ctf-1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
