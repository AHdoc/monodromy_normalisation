#ifndef GROUP_H
#define GROUP_H

#include<iostream>
#include<cstdio>
#include<cstring>

#include"myassert.h"

using namespace std;

/*
typedef long long LL;

struct SL2Z{
	LL e[2][2];
	void print(){
		cout<<" / "<<e[0][0]<<" , "<<e[0][1]<<" \\ \n";
		cout<<" \\ "<<e[1][0]<<" , "<<e[1][1]<<" / \n";
	}
};

SL2Z operator *(const SL2Z m1,const SL2Z m2){
	SL2Z m3;
	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++){
			m3.e[i][j]=0;
			for(int k=0;k<2;k++)
				m3.e[i][j]+=m1.e[i][k]*m2.e[k][j];
		}
	return m3;
}

SL2Z get(int a,int b,int c,int d){
	SL2Z ret;
	ret.e[0][0]=a;
	ret.e[0][1]=b;
	ret.e[1][0]=c;
	ret.e[1][1]=d;
	return ret;
}
*/

/*************************************************************/

struct Ga3b2{
	vector<string> e;

	Ga3b2(vector<string> e0): e(e0){}
	Ga3b2(): e({}){}

	void init(string s){
		e.clear();
		if(s!="1"){
			for(int i=0;i<s.size();i++){
				if(s[i]=='b'){
					if(i+1<s.size()) myassert(s[i+1]=='a',"'b' must be followed by 'a'");
					e.push_back("b");
				}if(s[i]=='a'){
					if(i+1<s.size()){
						myassert(s[i+1]=='b' or s[i+1]=='^',"'a' must be followed by either 'b' or '^'");
						if(s[i+1]=='b') e.push_back("a");
						else if(s[i+1]=='^'){
							myassert(i+2<s.size() && s[i+2]=='2',"'a^' must be followed by '2'");
							if(i+3<s.size()) myassert(s[i+3]=='b',"'a^2' must be followed by 'b'");
							e.push_back("a^2");
							i+=2;
						}
					}else e.push_back("a");
				}
			}
		}
	}
	
	void be_identity(){
		init("1");
	}
	
	Ga3b2 inv(){
		Ga3b2 g;
		g.e.clear();
		for(int i=e.size()-1;i>=0;i--){
			if(e[i]=="b") g.e.push_back("b");
			if(e[i]=="a") g.e.push_back("a^2");
			if(e[i]=="a^2") g.e.push_back("a"); 
		}
		myassert(e.size()==g.len(),"g.len == g.inv().len()");
		return g;
	}
	
	friend bool operator==(const Ga3b2 g1, const Ga3b2 g2){
		return g1.e==g2.e;
	}

	friend bool operator!=(const Ga3b2 g1, const Ga3b2 g2){
		return g1.e!=g2.e;
	}
	
	friend Ga3b2 operator *(const Ga3b2 g1,const Ga3b2 g2){
		int i=g1.len()-1,j=0;
		string r="1";
		while(i>=0 && j<g2.len()){
			if(g1.e[i]=="b" || g2.e[j]=="b"){
				if(g1.e[i]=="b" && g2.e[j]=="b"){
					--i; ++j;
				}else break;
			}else{
				if(g1.e[i]!=g2.e[j]){
					--i; ++j;
				}else{
					if(g1.e[i]=="a") r="a^2";
					else r="a";
					--i; ++j;
					break;
				}
			}
		}
		Ga3b2 g3;
		g3.e.clear();
		for(int k=0;k<=i;k++) g3.e.push_back(g1.e[k]);
		if(r!="1") g3.e.push_back(r);
		for(int k=j;k<g2.len();k++) g3.e.push_back(g2.e[k]);
		return g3;
	}

	void conjugation(Ga3b2 v){ // take the conjugation of this with v
		*this = v.inv() * (*this) * v;
	}

	int len() const{
		return e.size();
	}
	
	bool is_power_of_a() const{
		return len()==1 && (e[0]=="a" || e[0]=="a^2");
	}

	string print() const{
		if(e.size()==0) return "1";
		else{
			string s;
			for(int i=0;i<len();i++) s+=e[i];
			return s;
		}
	}
};

istream& operator>>(istream& is, Ga3b2& g){
	string s;
	is >> s;
	g.init(s);
    return is;
}

ostream& operator<<(ostream& os, const Ga3b2& g){
	return os<<g.print();
}

Ga3b2 pow(const Ga3b2 g,int n){
	Ga3b2 h;
	h.e.clear();
	for(int i=0;i<n;i++) h=h*g;
	return h;
}

/*************************************************************/

Ga3b2 a({"a"}),a2({"a^2"}),b({"b"});
Ga3b2 s0({"a^2","b"}),s1({"a","b","a"}),s2({"b","a^2"}),t0({"b","a"}),t1({"a^2","b","a^2"}),t2({"a","b"});

Ga3b2 short_elems[9]={a, a2, b, s0, s1, s2, t0, t1, t2};
Ga3b2 almost_short_elems[19]={
	a, a2, b, s0, s1, s2, t0, t1, t2,
	Ga3b2({"b","a","b"}), Ga3b2({"b","a^2","b"}), Ga3b2({"a^2","b","a"}), Ga3b2({"a","b","a^2"}),
	Ga3b2({"a^2","b","a","b"}), Ga3b2({"a","b","a","b","a"}), Ga3b2({"b","a","b","a^2"}),
	Ga3b2({"b","a^2","b","a"}), Ga3b2({"a^2","b","a^2","b","a^2"}), Ga3b2({"a","b","a^2","b"}),
};

/*************************************************************/
/* short elements */

bool short_table[9][9]={
	{1, 1, 0, 1, 0, 0, 1, 1, 1},
	{1, 1, 0, 1, 1, 1, 0, 0, 1},
	{0, 0, 1, 0, 0, 1, 1, 0, 0},
	{0, 1, 1, 0, 0, 1, 1, 0, 0},
	{0, 1, 0, 1, 0, 0, 0, 1, 0},
	{1, 1, 0, 0, 1, 0, 0, 0, 1},
	{1, 1, 0, 1, 0, 0, 0, 1, 0},
	{1, 0, 0, 0, 1, 0, 0, 0, 1},
	{1, 0, 1, 0, 0, 1, 1, 0, 0}
};

int is_short(Ga3b2 g){ // return -1 if not
	for(int i=0;i<9;i++)
		if(short_elems[i]==g)
			return i;
	return -1;
}

pair<Ga3b2,Ga3b2> get_expression(Ga3b2 g){
	if(is_short(g)!=-1 || g.len()==0)
		return make_pair(g,Ga3b2());
	if(g.len()>=3){
		Ga3b2 tau=Ga3b2({g.e[g.len()/2]});
		Ga3b2 q=Ga3b2(vector<string>(g.e.begin(),g.e.begin()+g.len()/2)).inv();
		if(q.inv()*tau*q==g) return make_pair(tau, q);
	}
	if(g.len()>=5){
		Ga3b2 tau=Ga3b2({g.e[g.len()/2-1],g.e[g.len()/2],g.e[g.len()/2+1]});
		Ga3b2 q=Ga3b2(vector<string>(g.e.begin(),g.e.begin()+g.len()/2-1)).inv();
		if(q.inv()*tau*q==g) return make_pair(tau, q);
	}
}

Ga3b2 get_tau_in_expression(Ga3b2 g){return get_expression(g).first;}
Ga3b2 get_Q_in_expression(Ga3b2 g){return get_expression(g).second;}

int S_complexity(Ga3b2 g){
	Ga3b2 tau=get_tau_in_expression(g);
	Ga3b2 Q=get_Q_in_expression(g);
	return Q.len();
}

/*************************************************************/
/* almost short elements */

bool almost_short_table[19][19]={
{1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0,0,0,0},
{1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,1},
{1,1,1,1,0,1,1,0,1,1,1,0,0,0,0,1,1,0,0},
{1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0},
{1,1,0,1,0,0,0,1,1,0,0,1,0,1,0,0,0,1,0},
{1,1,1,1,1,0,0,1,1,0,0,0,1,0,1,0,0,0,1},
{1,1,1,1,1,0,0,1,1,0,0,1,0,1,0,0,0,1,0},
{1,1,0,1,1,0,0,0,1,0,0,0,1,0,1,0,0,0,1},
{1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0},
{0,0,1,0,0,1,1,0,0,1,1,0,0,0,0,0,1,0,0},
{0,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,0},
{1,1,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0},
{1,1,0,0,1,0,0,0,1,0,0,1,1,0,1,0,0,0,1},
{0,0,1,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0},
{0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0},
{1,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,1},
{0,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0},
{0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,1},
{0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,1,1,0,0}
};

int is_almost_short(Ga3b2 g){ // return -1 if not
	for(int i=0;i<19;i++)
		if(almost_short_elems[i]==g)
			return i;
	return -1;
}

pair<Ga3b2,Ga3b2> get_expression_2(Ga3b2 g){
	if(is_almost_short(g)!=-1 || g.len()==0)
		return make_pair(g,Ga3b2());
	if(g.len()>=5){
		Ga3b2 tau=Ga3b2({g.e[g.len()/2-1],g.e[g.len()/2],g.e[g.len()/2+1]});
		Ga3b2 q=Ga3b2(vector<string>(g.e.begin(),g.e.begin()+g.len()/2-1)).inv();
		if(q.inv()*tau*q==g) return make_pair(tau, q);
	}
	if(g.len()>=7){
		Ga3b2 tau=Ga3b2({g.e[g.len()/2-2],g.e[g.len()/2-1],g.e[g.len()/2],g.e[g.len()/2+1],g.e[g.len()/2+2]});
		Ga3b2 q=Ga3b2(vector<string>(g.e.begin(),g.e.begin()+g.len()/2-2)).inv();
		if(q.inv()*tau*q==g) return make_pair(tau, q);
	}
	myassert(false,"get_expression_2 can only handle conjugates of almost short elements");
}

Ga3b2 get_tau_in_expression_2(Ga3b2 g){return get_expression_2(g).first;}
Ga3b2 get_Q_in_expression_2(Ga3b2 g){return get_expression_2(g).second;}

int S2_complexity(Ga3b2 g){
	Ga3b2 tau=get_tau_in_expression_2(g);
	Ga3b2 Q=get_Q_in_expression_2(g);
	if(tau.len()==5 && Q.len()==0) return 1;
	else return Q.len()*2;
}

/*************************************************************/

#endif
