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

using namespace std;

//void normal_Ga3b2(Tuple g){
//}

void test_Ga3b2(){
	cout<<"This is a quick test for struct Ga3b2.\n";
	cout<<"Input two elements in Ga3b2, for instance, aba^2 a^2ba^2.\n";

	Ga3b2 g, h;
	cin>>g>>h;
	cout<<"g="<<g<<"\n";
	cout<<"h="<<h<<"\n";
	cout<<"g.inv()="<<g.inv()<<"\n";
	cout<<"h.inv()="<<h.inv()<<"\n";
	cout<<"g*h="<<g*h<<"\n";
	cout<<"pow(g,2)*pow(h,3)="<<pow(g,2)*pow(h,3)<<"\n";
}

void test_tuple(){
	cout<<"This is a quick test for struct Tuple and iteratedTuple.\n";
	int n;
	cin>>n;
	vector<Ga3b2> gvec;
	gvec.clear();
	for(int i=0;i<n;i++){
		Ga3b2 gveci;
		cin>>gveci;
		gvec.push_back(gveci);
	}
	Tuple<Ga3b2> g;
	g.init(gvec);
	cout<<g<<"\n";
}

int main(){
	//test_Ga3b2();
	test_tuple();
}