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

#include"myassert.h"
#include"group.h" 
#include"tuple.h"
#include"iteratedtuple.h"

#include"contraction_restoration.h"


using namespace std;

#include"normal_form.cpp"
#include"short.cpp"
#include"almost_short.cpp"
#include"test.cpp"


int main(){
	//test_Ga3b2();
	//test_tuple();
	//test_CR();

	//multi_test_normalize_short(false);
	test_normalize_almost_short(10,true);
}
