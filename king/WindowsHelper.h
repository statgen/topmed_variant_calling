#ifndef __WINDOWSHELPER_H__
#define __WINDOWSHELPER_H__

#ifndef __WIN32__
inline void WildCardArguments(int argc, char ** argv) { }
#else
void WildCardArguments(int & argc, char ** & argv);
#endif

#endif


