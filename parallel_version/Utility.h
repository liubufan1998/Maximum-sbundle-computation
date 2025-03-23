#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <tuple>
#include <queue>
#include <set>
#include <omp.h>
#include <chrono>

#define NDEBUG 
#include <cassert>
using namespace std;



#define _HEURISTIC_FIRST_STAGE_	 
#define _HEURISTIC_SECOND_STAGE_ 
#define _HEURISTIC_THIRD_STAGE_	 
#define _LOCAL_VC_DELETE_		 
#define _VP_BOUND_				 
#define _BETTER_SEQEUNCE_		 

#define _SUPPORT_BOUND_		
#define _SBUNDLE_BRANCHING_ 

#define _IE_FRAMEWORK_ 

using ui = unsigned int;  
using ept = unsigned int; 

#define pb push_back
#define mp make_pair

#define mmax(a, b) ((a) > (b) ? (a) : (b))
#define mmin(a, b) ((a) < (b) ? (a) : (b))

class Utility
{
public:
	static FILE *open_file(const char *file_name, const char *mode)
	{
		FILE *f = fopen(file_name, mode);
		if (f == nullptr)
		{
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}

		return f;
	}

	static std::string integer_to_string(long long number)
	{
		std::vector<ui> sequence;
		if (number == 0)
			sequence.pb(0);
		while (number > 0)
		{
			sequence.pb(number % 1000);
			number /= 1000;
		}

		char buf[5];
		std::string res;
		for (ui i = sequence.size(); i > 0; i--)
		{
			if (i == sequence.size())
				sprintf(buf, "%u", sequence[i - 1]);
			else
				sprintf(buf, ",%03u", sequence[i - 1]);
			res += std::string(buf);
		}
		return res;
	}
};

#endif
