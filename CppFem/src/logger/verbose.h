#ifndef verbose_h
#define verbose_h

namespace fem
{
#define VERBOSE             // please activate or de-activate this line

#define VERBOSE_PRINTS(str, str1) FEM_LOG_INFO("%-30s %6s\n", str, str1);
#define VERBOSE_PRINT0(str, number) FEM_LOG_DEBUG("%-30s %6d\n", str, number);

#define TIME_REPORT        // please activate or de-activate this line

#ifndef DETAILED_REPORT
	//#define DETAILED_REPORT  // please activate or de-activate this line
	//#define VERBOSE          // please activate or de-activate this line
#define TIME_REPORT       // please activate or de-activate this line
#endif
} // end namespace fem
#endif // verbose_h
