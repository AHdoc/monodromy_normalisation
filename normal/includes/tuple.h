#ifndef TUPLE_H
#define TUPLE_H

#include<iostream>
#include<cstdio>
#include<cstring>

#include"myassert.h"
#include"group.h"

using namespace std;

template<typename G>
struct Tuple{
	vector<G> e;
	int len() const{
		return e.size();
	}
	void init(vector<G> e0){
		e=e0;
	}
	string print() const{
		string s="(";
		for(int i=0;i<len();i++){
			if(i>0) s+=",";
			s+=e[i].print();
		}
		s+=")";
		return s;
	}
};

template<typename G>
std::ostream& operator<<(std::ostream& os, const Tuple<G>& g){
	return os<<g.print();
}

#endif
