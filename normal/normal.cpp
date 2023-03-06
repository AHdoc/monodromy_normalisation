// g++ normal.cpp -O2 -o normal -I includes -I algos

#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<vector>
#include<set>
#include<map>
#include<ctime>
#include<algorithm>

using namespace std;

#include"myassert.h"
#include"group.h" 
#include"tuple.h"
#include"iteratedtuple.h"

#include"contraction_restoration.h"

#include"normal_form.cpp"
#include"randomm.cpp"
#include"short.cpp"
#include"almost_short.cpp"
#include"test.cpp"


int main(){
	//test_Ga3b2();
	//test_tuple();
	//test_CR();

	//multi_test_normalize_short(false);
	multi_test_normalize_almost_short(false);
}
