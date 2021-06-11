#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

inline void fill_uchar_zero(unsigned char* uchar_array_, size_t size) { unsigned char*iptr = &(uchar_array_[size]); while (uchar_array_<iptr){ *uchar_array_++=0;}}
inline void fill_uchar_ones(unsigned char* uchar_array_, size_t size) { unsigned char*iptr = &(uchar_array_[size]); while (uchar_array_<iptr){ *uchar_array_++=255;}}
inline void fill_long_zero(long* long_array_, size_t size) { long* lptr = &(long_array_[size]); while (long_array_ < lptr) { *long_array_++ = 0;} }
