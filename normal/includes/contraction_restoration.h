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

	int Restoration(){
		list<vector<int>> F2;
		int l=-1,r=-1;
		for(auto it:F){
			if(it[0]==2){
				l=it[1];
				r=it[2];
				F2.clear();
			}else{
				F2.push_back(it);
			}
		}
		if(l==-1 && r==-1) return -1;
		else{
			while(F.back()!=vector<int>{2,l,r}) F.pop_back();
			F.pop_back();
		}

		int k=l,m=h.len();
		int m2=m+(r-l);
		for(auto it:F2){ // recall: F.push_back({1,h.len(),i,epsilon});
			myassert(it[0]==1,"All elements in F2 are elementary transformations");
			myassert(it[1]==m,"All elementary transformations in F2 are associated to the same length");
			int i=it[2],epsilon=it[3];
			if(1<=i && i<=k-2)
				F.push_back({1,m2,i,epsilon});
			else if(k+1<=i && i<=m)
				F.push_back({1,m2,i+(r-l),epsilon});
			else if(i==k-1){
				for(int j=k-1;j<=k-1+(r-l);j++)
					F.push_back({1,m2,j,epsilon});
				k=k-1;
			}else if(i==k){
				for(int j=k+(r-l);j>=k;j--)
					F.push_back({1,m2,j,epsilon});
				k=k+1;
			}
		}

		iteratedTuple<G> H2=*H.e[k-1];
		myassert(H2.h()>=1,"Restore an iterated tuple of height at least 2");

		h.e.erase(h.e.begin()+k-1);
		H.e.erase(H.e.begin()+k-1);
		for(int i=H2.len();i>=1;i--){
			auto it=H2.e[i-1];
			h.e.insert(h.e.begin()+k-1,it->ev());
			H.e.insert(H.e.begin()+k-1,it);
		}

		return k;
	}
};

#endif
