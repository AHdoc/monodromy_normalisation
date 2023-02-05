#ifndef CONTRACTION_RESTORATION_H
#define CONTRACTION_RESTORATION_H

#include<iostream>
#include<cstdio>
#include<cstring>
#include<list>

#include"myassert.h"
#include"group.h"
#include"tuple.h"
#include"iteratedtuple.h"

using namespace std;

template<typename G>
struct CR{
	Tuple<G> h;
	iteratedTuple<G> H;

	list<vector<int>> F;

	void init(Tuple<G> g){
		h=g;
		H.init(g);
		F.clear();
	}

	void Elementary_transformation(int i,int epsilon){
		h.Elementary_transformation(i,epsilon);
		H.Elementary_transformation(i,epsilon);

		F.push_back({1,h.len(),i,epsilon});
	}

	void Contraction(int l,int r){
		myassert(h.len()==H.len(),"h.len()==H.len()");
		int m=h.len();
		myassert(1<=l && l<r && r<=m,"For contraction, 1<=l<r<=m");

		G prod;
		prod.be_identity();
		for(int i=l;i<=r;i++){
			prod=prod*h.e[i-1];
		}
		h.e.erase(h.e.begin()+l,h.e.begin()+r);
		h.e[l-1]=prod;

		struct iteratedTuple<G>* H2 = new struct iteratedTuple<G>();
		H2->e.clear();
		for(int i=l;i<=r;i++){
			H2->e.push_back(H.e[i-1]);
		}
		H.e.erase(H.e.begin()+l,H.e.begin()+r);
		H.e[l-1]=H2;

		F.push_back({2,l,r});
	}
};

#endif
