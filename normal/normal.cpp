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


int main(int argc, char *argv[]){
	srand (time(NULL));

	if(argc==1 || string(argv[1])=="-h"){
		cout<<"Options:\n";
		cout<<"  -h             Display this information.\n";
		cout<<"\n";
		cout<<"  -test\n";
		cout<<"    -Ga3b2       Tests for Ga3b2.\n";
		cout<<"    -Tuple       Tests for Tuple.\n";
		cout<<"    -CR          Tests for CR.\n";
		cout<<"\n";
		cout<<"    -short       Tests for normalize_short.\n";
		cout<<"    -almostshort Tests for normalize_almost_short.\n";
	}
	if(argc>=2){
		if(string(argv[1])=="-test"){
			if(argc>=3){
				if(string(argv[2])=="-Ga3b2")
					test_Ga3b2();
				if(string(argv[2])=="-Tuple")
					test_tuple();
				if(string(argv[2])=="-CR")
					test_CR();
				
				bool details=(argc>=4 && string(argv[3])=="true");
				if(string(argv[2])=="-short"){
					multi_test_normalize_short(details);
				}
				if(string(argv[2])=="-almostshort"){
					multi_test_normalize_almost_short(details);
				}
			}
		}
	}
}