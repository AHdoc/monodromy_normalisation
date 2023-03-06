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
		if(s.size()<=130) return s;
		else              return string(s,0,130)+"... of length "+to_string(s.size());
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

/*************************************************************/

Tuple<Ga3b2> exceptional_tuples[18]={
	Tuple<Ga3b2>({a2,s0,s2}), Tuple<Ga3b2>({a,t2,t0}), Tuple<Ga3b2>({a,s0,s0,s2,s0}), Tuple<Ga3b2>({a2,t0,t2,t0,t0}), Tuple<Ga3b2>({b,s0,s2,s0}), Tuple<Ga3b2>({b,t0,t2,t0}),
	Tuple<Ga3b2>({a,b,s2}), Tuple<Ga3b2>({a2,b,t0}), Tuple<Ga3b2>({a,t2,t0,b,t0,t2,t0}), Tuple<Ga3b2>({a2,s0,s2,b,s0,s2,s0}),
	Tuple<Ga3b2>({a,a,s0,s2}), Tuple<Ga3b2>({a2,a2,t2,t0}), Tuple<Ga3b2>({a2,a2,s0,s0,s2,s0}), Tuple<Ga3b2>({a,a,t0,t2,t0,t0}),
	Tuple<Ga3b2>({a2,a2,b,s2}), Tuple<Ga3b2>({a,a,b,t0}), Tuple<Ga3b2>({a2,a2,t2,t0,b,t0,t2,t0}), Tuple<Ga3b2>({a,a,s0,s2,b,s0,s2,s0})
};

/*************************************************************/

bool each_component_is_short(Tuple<Ga3b2> g){
	for(Ga3b2 x:g.e)
		if(is_short(x)>=0);
		else return false;
	return true;
}

int S_complexity(Tuple<Ga3b2> g){
	int ret=0;
	for(int i=0;i<g.len();i++) ret+=S_complexity(g.e[i]);
	return ret;
}

/*************************************************************/

bool each_component_is_almost_short(Tuple<Ga3b2> g){
	for(Ga3b2 x:g.e)
		if(is_almost_short(x)>=0);
		else return false;
	return true;
}

int S2_complexity(Tuple<Ga3b2> g){
	int ret=0;
	for(int i=0;i<g.len();i++) ret+=S2_complexity(g.e[i]);
	return ret;
}
#endif
