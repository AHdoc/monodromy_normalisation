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
#include"test.cpp"


int main(){
	//test_Ga3b2();
	//test_tuple();
	//test_CR();

	multi_test_normalize_short(false);

	//auto ret=search(Tuple<Ga3b2>({b,t0,t2,t0,a2,a2,t1,t2,a,a,a}),Tuple<Ga3b2>({a,a2}));
	//cout<<"ret.first="<<ret.first<<"\n";
}

// a^2b ba^2 a^2ba^2 a^2b ba^2 ba^2 a^2b ba^2
