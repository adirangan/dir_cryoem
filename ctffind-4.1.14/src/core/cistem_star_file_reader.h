class  cisTEMStarFileReader {

	int current_position_in_stack;
	int current_column;

	int 	position_in_stack_column;
	int 	image_is_active_column;
	int	 	psi_column;
	int		theta_column;
	int		phi_column;
	int		x_shift_column;
	int		y_shift_column;
	int		defocus_1_column;
	int		defocus_2_column;
	int		defocus_angle_column;
	int		phase_shift_column;
	int		occupancy_column;
	int		logp_column;
	int		sigma_column;
	int		score_column;
	int		score_change_column;
	int		pixel_size_column;
	int		microscope_voltage_kv_column;
	int		microscope_spherical_aberration_mm_column;
	int		beam_tilt_x_column;
	int		beam_tilt_y_column;

public:

	wxString    filename;
	wxTextFile *input_file;

	bool using_external_array;

	ArrayOfcisTEMParameterLines *cached_parameters;

	cisTEMStarFileReader();
	~cisTEMStarFileReader();

	cisTEMStarFileReader(wxString wanted_filename, ArrayOfcisTEMParameterLines *alternate_cached_parameters_pointer = NULL);

	void Open(wxString wanted_filename, ArrayOfcisTEMParameterLines *alternate_cached_parameters_pointer = NULL);
	void Close();
	bool ReadFile(wxString wanted_filename, wxString *error_string = NULL, ArrayOfcisTEMParameterLines *alternate_cached_parameters_pointer = NULL);

	bool ExtractParametersFromLine(wxString &wanted_line, wxString *error_string = NULL);

	inline int   ReturnPositionInStack(int line_number) { return cached_parameters->Item(line_number).position_in_stack;}
	inline int   ReturnImageIsActive(int line_number) { return cached_parameters->Item(line_number).image_is_active;}
	inline float ReturnPhi(int line_number) { return cached_parameters->Item(line_number).phi;}
	inline float ReturnTheta(int line_number) { return cached_parameters->Item(line_number).theta;}
	inline float ReturnPsi(int line_number) { return cached_parameters->Item(line_number).psi;}
	inline float ReturnXShift(int line_number) { return cached_parameters->Item(line_number).x_shift;}
	inline float ReturnYShift(int line_number) { return cached_parameters->Item(line_number).y_shift;}
	inline float ReturnDefocus1(int line_number) { return cached_parameters->Item(line_number).defocus_1;}
	inline float ReturnDefocus2(int line_number) { return cached_parameters->Item(line_number).defocus_2;}
	inline float ReturnDefocusAngle(int line_number) { return cached_parameters->Item(line_number).defocus_angle;}
	inline float ReturnPhaseShift(int line_number) { return cached_parameters->Item(line_number).phase_shift;}
	inline int   ReturnLogP(int line_number) { return cached_parameters->Item(line_number).logp;}
	inline float ReturnSigma(int line_number) { return cached_parameters->Item(line_number).sigma;}
	inline float ReturnScore(int line_number) { return cached_parameters->Item(line_number).score;}
	inline float ReturnScoreChange(int line_number) { return cached_parameters->Item(line_number).score_change;}
	inline float ReturnPixelSize(int line_number) { return cached_parameters->Item(line_number).pixel_size;}
	inline float ReturnMicroscopekV(int line_number) { return cached_parameters->Item(line_number).microscope_voltage_kv;}
	inline float ReturnMicroscopeCs(int line_number) { return cached_parameters->Item(line_number).microscope_spherical_aberration_mm;}
	inline float ReturnBeamTiltX(int line_number) { return cached_parameters->Item(line_number).beam_tilt_x;}
	inline float ReturnBeamTiltY(int line_number) { return cached_parameters->Item(line_number).beam_tilt_y;}


};
