
void		swapbytes(unsigned char* v, size_t n);
void		swapbytes(size_t size, unsigned char* v, size_t n);

bool GetMRCDetails(const char *filename, int &x_size, int &y_size, int &number_of_images);

inline void ZeroBoolArray(bool *array_to_zero, int size_of_array)
{
	for (int counter = 0; counter < size_of_array; counter++)
	{
		array_to_zero[counter] = false;
	}
}

int sizeCanBe4BitK2SuperRes(int nx, int ny); // provided by David Mastronarde

inline void ZeroIntArray(int *array_to_zero, int size_of_array)
{
	for (int counter = 0; counter < size_of_array; counter++)
	{
		array_to_zero[counter] = 0;
	}
}

inline void ZeroLongArray(long *array_to_zero, int size_of_array)
{
	for (int counter = 0; counter < size_of_array; counter++)
	{
		array_to_zero[counter] = 0;
	}
}

inline void ZeroFloatArray(float *array_to_zero, int size_of_array)
{
	for (int counter = 0; counter < size_of_array; counter++)
	{
		array_to_zero[counter] = 0.0;
	}
}

int ReturnSafeBinnedBoxSize(int original_box_size, float bin_factor);


float ReturnMagDistortionCorrectedPixelSize(float original_pixel_size, float major_axis_scale, float minor_axis_scale);

wxString ReturnSocketErrorText(wxSocketBase *socket_to_check);

inline void WriteToSocket	(	wxSocketBase *socket, const void * 	buffer, wxUint32 nbytes, bool die_on_error = false)
{

	bool should_abort = false;
	if (socket != NULL)
	{
		if (socket->IsOk() == true && socket->IsConnected() == true)
		{

#ifdef DEBUG
			//	socket->SetFlags(wxSOCKET_WAITALL | wxSOCKET_BLOCK);
			//if (socket->GetFlags() != (wxSOCKET_WAITALL) && socket->GetFlags() != (wxSOCKET_WAITALL | wxSOCKET_BLOCK)) 	{MyPrintWithDetails("Wait all flag not set!"); should_abort = true;}
			if (socket->GetFlags() != (wxSOCKET_WAITALL | wxSOCKET_BLOCK)) 	{MyPrintWithDetails("Wait all / block flag not set!"); should_abort = true;}
#endif

			//socket->SetTimeout(60);
			socket->WaitForWrite();
			socket->Write(buffer, nbytes);

			int number_of_retries = 0;
			while (socket->LastWriteCount() == 0 && number_of_retries < 10)
			{
				wxMilliSleep(100);
				socket->WaitForWrite();
				socket->Write(buffer, nbytes);
				number_of_retries++;
			}

			if (socket->LastWriteCount() != nbytes) {MyDebugPrintWithDetails("Socket didn't write all bytes! (%u / %u) with %i retries ", socket->LastWriteCount(), nbytes, number_of_retries); should_abort = true;}
			if (socket->Error() == true) {MyDebugPrintWithDetails("Socket has an error (%s) ", ReturnSocketErrorText(socket)); should_abort = true;}

			if (should_abort == true)
			{
#ifdef DEBUG
				wxIPV4address peer_address;
				socket->GetPeer(peer_address);

				wxPrintf("Failed socket is connected to : %s, (%s)\n", peer_address.Hostname(), peer_address.IPAddress());
				socket->Destroy();
				socket = NULL;
				wxFAIL;
				DEBUG_ABORT;
#else
				exit(-1);
#endif
			}
		}
		else
		{
			if (die_on_error == true) DEBUG_ABORT;
		}
	}
}

inline void ReadFromSocket	(	wxSocketBase *socket, void * 	buffer, wxUint32 nbytes, bool die_on_error = false)
{
	bool should_abort = false;
	if (socket != NULL)
	{
		if (socket->IsOk() == true && socket->IsConnected() == true)
		{
#ifdef DEBUG
			//	socket->SetFlags(wxSOCKET_WAITALL | wxSOCKET_BLOCK);
			//if (socket->GetFlags() != (wxSOCKET_WAITALL) && socket->GetFlags() != (wxSOCKET_WAITALL | wxSOCKET_BLOCK)) 	{MyPrintWithDetails("Wait all flag not set!"); should_abort = true;}
			if (socket->GetFlags() != (wxSOCKET_WAITALL | wxSOCKET_BLOCK)) 	{MyPrintWithDetails("Wait all / block flag not set!"); should_abort = true;}
#endif

			//socket->SetTimeout(60);
			socket->WaitForRead();
			socket->Read(buffer, nbytes);

			int number_of_retries = 0;
			while (socket->LastReadCount() == 0 && number_of_retries < 10)
			{
				wxMilliSleep(100);
				socket->WaitForRead();
				socket->Read(buffer, nbytes);
				number_of_retries++;
			}

			if (socket->LastReadCount() != nbytes) {MyDebugPrintWithDetails("Socket didn't read all bytes! (%u / %u) with %i retries ", socket->LastReadCount(), nbytes, number_of_retries); should_abort = true;}
			if (socket->Error() == true) {MyDebugPrintWithDetails("Socket has an error (%s) ", ReturnSocketErrorText(socket)); should_abort = true;}
			if (should_abort == true)
			{
#ifdef DEBUG
				wxIPV4address peer_address;
				socket->GetPeer(peer_address);

				wxPrintf("Failed socket is connected to : %s, (%s)\n", peer_address.Hostname(), peer_address.IPAddress());
				socket->Destroy();
				socket = NULL;
				wxFAIL;
				DEBUG_ABORT;
#else
				exit(-1);
#endif

			}
		}
		else
		{
			if (die_on_error == true) DEBUG_ABORT;
		}
	}
}



inline void ZeroDoubleArray(double *array_to_zero, int size_of_array)
{
	for (int counter = 0; counter < size_of_array; counter++)
	{
		array_to_zero[counter] = 0.0;
	}
}

inline bool IsEven(int number_to_check)
{
	  if ( number_to_check % 2== 0 ) return true;
	  else return false;
}

inline bool DoesFileExist(wxString filename)
{
    std::ifstream file_to_check (filename.c_str());

    if(file_to_check.is_open()) return true;
    return false;
}

inline bool DoesFileExistWithWait(wxString filename, int max_wait_time_in_seconds)
{

    if ( ! DoesFileExist(filename) )
	{
		for (int wait_counter = 0; wait_counter < max_wait_time_in_seconds; wait_counter ++)
		{
			wxSleep(1);
			if (DoesFileExist(filename)) break;
		}
	}

    return DoesFileExist(filename);

}

inline float rad_2_deg(float radians)
{
  return radians / (PI / 180.);
}

inline float deg_2_rad(float degrees)
{
  return degrees * PI / 180.;
}

inline float sinc(float radians)
{
	if (radians == 0.0) return 1.0;
	if (radians >= 0.01) return sinf(radians) / radians;
	float temp_float = radians * radians;
	return 1.0 - temp_float / 6.0 + temp_float * temp_float / 120.0;
}

inline double myround(double a)
{
	if (a > 0) return double(long(a + 0.5));
	else return double(long(a - 0.5));
}

inline float myround(float a)
{
	if (a > 0) return float(int(a + 0.5));	else return float(int(a - 0.5));
}

inline int myroundint(double a)
{
	if (a > 0) return int(a + 0.5); else return int(a - 0.5);
}

inline int myroundint(float a)
{
	if (a > 0) return int(a + 0.5);	else return int(a - 0.5);
}

inline bool IsOdd(int number)
{
	if ((number & 1) == 0) return false;
	else return true;
}

wxArrayString ReturnIPAddress();
wxString ReturnIPAddressFromSocket(wxSocketBase *socket);

void SendwxStringToSocket(wxString *string_to_send, wxSocketBase *socket);
wxString ReceivewxStringFromSocket(wxSocketBase *socket);

inline float ReturnPhaseFromShift(float real_space_shift, float distance_from_origin, float dimension_size)
{
	return real_space_shift * distance_from_origin * 2.0 * PI / dimension_size;
}

inline std::complex<float> Return3DPhaseFromIndividualDimensions( float phase_x, float phase_y, float phase_z)
{
	float temp_phase = -phase_x-phase_y-phase_z;

	return cosf(temp_phase) + sinf(temp_phase) * I;
}

inline bool DoublesAreAlmostTheSame(double a, double b)
{
	return (fabs(a-b) < 0.000001);
}

inline bool FloatsAreAlmostTheSame(float a, float b)
{
	return (fabs(a-b) < 0.0001);
}

inline bool InputIsATerminal()
{
   return isatty(fileno(stdin));
};

inline bool OutputIsAtTerminal()
{
    return isatty(fileno(stdout));
}

float CalculateAngularStep(float required_resolution, float radius_in_angstroms);

int ReturnClosestFactorizedUpper(int wanted_int, int largest_factor, bool enforce_even = false);
int ReturnClosestFactorizedLower(int wanted_int, int largest_factor, bool enforce_even = false);

/*
 *
 * String manipulations
 *
 */
std::string FilenameReplaceExtension(std::string filename, std::string new_extension);
std::string FilenameAddSuffix(std::string filename, std::string suffix_to_add);
void SplitFileIntoDirectoryAndFile(wxString &input_file, wxString &output_directory, wxString &output_file);

void Allocate2DFloatArray(float **&array, int dim1, int dim2);
void Deallocate2DFloatArray(float **&array, int dim1);

void CheckSocketForError(wxSocketBase *socket_to_check);

inline wxString BoolToYesNo(bool b)
{
  return b ? "Yes" : "No";
}

long ReturnFileSizeInBytes(wxString filename);

inline float kDa_to_Angstrom3(float kilo_daltons)
{
  return kilo_daltons * 1000.0 / 0.81;
}

inline bool IsAValidSymmetry(wxString *string_to_check)
{
	long junk;
	wxString buffer_string = *string_to_check;
	buffer_string.Trim();
	buffer_string.Trim(false);


	if (string_to_check->StartsWith("C") == false && string_to_check->StartsWith("c") == false && string_to_check->StartsWith("D") == false && string_to_check->StartsWith("d") == false && string_to_check->StartsWith("I") == false && string_to_check->StartsWith("i") == false && string_to_check->StartsWith("O") == false && string_to_check->StartsWith("o") == false && string_to_check->StartsWith("T") == false && string_to_check->StartsWith("t") == false) return false;
	if (buffer_string.Mid(1).IsEmpty() == true && (buffer_string.StartsWith("I") == true || buffer_string.StartsWith("i") == true || buffer_string.StartsWith("T") == true || buffer_string.StartsWith("t") == true || buffer_string.StartsWith("O") == true || buffer_string.StartsWith("o") == true )) return true;


	return string_to_check->Mid(1).ToLong(&junk);

}

float ReturnSumOfLogP(float logp1, float logp2, float log_range);

int ReturnNumberofAsymmetricUnits(wxString symmetry);

inline float ConvertProjectionXYToThetaInDegrees(float x, float y) // assumes that max x,y is 1
{
	return rad_2_deg(asin(sqrtf(powf(x,2) + powf(y,2))));
}

inline float ConvertXYToPhiInDegrees(float x, float y)
{
	if (x == 0 && y==0) return 0;
	else return rad_2_deg(atan2f(y, x));
}

std::vector<size_t> rankSort(const float* v_temp, const size_t size);
std::vector<size_t> rankSort(const std::vector<float>& v_temp);

wxString StringFromSocketCode(unsigned char *socket_input_buffer);

int CheckNumberOfThreads(int number_of_threads);

// From David Mastronarde
int ReturnAppropriateNumberOfThreads(int optimalThreads);
int ReturnThreadNumberOfCurrentThread();

double cisTEM_erfinv(double x);
double cisTEM_erfcinv(double x);
