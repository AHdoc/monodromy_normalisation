#ifndef ITERATEDTUPLE_H
#define ITERATEDTUPLE_H

#include<iostream>
#include<cstdio>
#include<cstring>

#include"myassert.h"
#include"group.h"

using namespace std;

template<typename G>
struct iteratedTuple{
	vector<struct iteratedTuple*> e;
	G val;

	int h() const{
		int maxh=0;
		for(auto ei:e) maxh=max(maxh,1+ei->h());
		return maxh;
	}
	
	G ev() const{
		if(e.empty()) return val;
		else{
			G ret;
			ret.be_identity();
			for(auto ei:e) ret=ret*ei->ev();
			return ret;
		}
	}
	
	int len() const{
		return e.size();
	}
	
	void init(G val0){
		e.clear();
		val=val0;
	}
	
	void init(Tuple<G> g){
		e.clear();
		for(int i=0;i<g.len();i++){
			struct iteratedTuple<G>* component = new struct iteratedTuple<G>();
			component->init(g.e[i]);
			e.push_back(component);
		}
	}

	void conjugation(G v){
		if(e.empty()) val.conjugation(v);
		else{
			for(auto ei:e) ei->conjugation(v);
		}
	}
	
	string print() const{
		if(h()==0) return val.print();
		else{
			string s="(";
			for(int i=0;i<len();i++){
				if(i>0) s+=",";
				s+=e[i]->print();
			}
			s+=")";
			return s;
		}
	}

	void Elementary_transformation(int i,int epsilon){
		myassert(1<=i && i<=len()-1 && (epsilon==1 || epsilon==-1),"a proper Elementary transformation");
		if(epsilon==1){
			e[i-1]->conjugation(e[i]->ev());
			swap(e[i-1],e[i]);
		}else{
			e[i]->conjugation((e[i-1]->ev()).inv());
			swap(e[i-1],e[i]);
		}
	}
};

template<typename G>
std::ostream& operator<<(std::ostream& os, const iteratedTuple<G>& g){
	return os<<g.print();
}

#endif
