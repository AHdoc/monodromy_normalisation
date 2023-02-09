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
	
	Tuple(vector<G> e0): e(e0){}
	Tuple(): e({}){}

	int len() const{
		return e.size();
	}
	
	void clear(){
		e.clear();
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

	void Elementary_transformation(int i,int epsilon){
		myassert(1<=i && i<=len()-1 && (epsilon==1 || epsilon==-1),"a proper Elementary transformation");
		if(epsilon==1){
			e[i-1].conjugation(e[i]);
			swap(e[i-1],e[i]);
		}else{
			e[i].conjugation(e[i-1].inv());
			swap(e[i-1],e[i]);
		}
	}
};

template<typename G>
std::ostream& operator<<(std::ostream& os, const Tuple<G>& g){
	return os<<g.print();
}

template<typename G>
bool operator==(const Tuple<G> g1, const Tuple<G> g2){
	return g1.e==g2.e;
}

template<typename G>
Tuple<G> bullet(Tuple<G> g1,Tuple<G> g2){
	Tuple<G> g;
	g.e=g1.e;
	g.e.insert(g.e.end(),g2.e.begin(),g2.e.end());
	return g;
}

#endif
