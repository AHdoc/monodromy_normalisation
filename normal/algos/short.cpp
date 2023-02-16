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

bool each_component_is_short(Tuple<Ga3b2> g){
	for(Ga3b2 x:g.e)
		if(is_short(x)>=0);
		else return false;
	return true;
}

int S_complexity(Ga3b2 g){
	if(is_short(g)!=-1 || g.len()==0) return 0;
	myassert(g.len()%2==1,"len() of a long element in Ga3b2 is odd");
	if(g.e[g.len()/2]=="a" || g.e[g.len()/2]=="a^2") return g.len()/2;
	else{
		if(g.len()>=3 && g.e[g.len()/2-1]==g.e[g.len()/2+1]) return g.len()/2-1;
		else                                                 return g.len()/2;
	}
}

int S_complexity(Tuple<Ga3b2> g){
	int ret=0;
	for(int i=0;i<g.len();i++) ret+=S_complexity(g.e[i]);
	return ret;
}

pair<Ga3b2,Ga3b2> get_expression(Ga3b2 g){
	if(is_short(g)!=-1 || g.len()==0)
		return make_pair(g,Ga3b2());
	if(g.e[g.len()/2]=="a" || g.e[g.len()/2]=="a^2")
		return make_pair(
			Ga3b2({g.e[g.len()/2]}),
			Ga3b2(vector<string>(g.e.begin(),g.e.begin()+g.len()/2)).inv()
		);
	else{
		if(g.len()>=3 && g.e[g.len()/2-1]==g.e[g.len()/2+1])
			return make_pair(
				Ga3b2({g.e[g.len()/2-1],g.e[g.len()/2],g.e[g.len()/2+1]}),
				Ga3b2(vector<string>(g.e.begin(),g.e.begin()+g.len()/2-1)).inv()
			);
		else
			return make_pair(
				Ga3b2({g.e[g.len()/2]}),
				Ga3b2(vector<string>(g.e.begin(),g.e.begin()+g.len()/2)).inv()
			);
	}
}

Ga3b2 get_tau_in_expression(Ga3b2 g){return get_expression(g).first;}
Ga3b2 get_Q_in_expression(Ga3b2 g){return get_expression(g).second;}

/*****/

/*--Operation 1--------------*/
bool Operation1(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=1;i+1<=n;i++){
		Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i];
		Ga3b2 tau_1=get_tau_in_expression(g1),Q_1=get_Q_in_expression(g1);
		Ga3b2 tau_2=get_tau_in_expression(g2),Q_2=get_Q_in_expression(g2);
		int row=is_short(tau_1), col=is_short(tau_2); // row==-1 iff tau_1==identity
		if(Q_1==Q_2 && row!=-1 && col!=-1 && short_table[row][col]){
			if(Q_1.len()==0 || (tau_1*tau_2).len()==0 || (tau_1.is_power_of_a() && tau_2.is_power_of_a())){
				M->Contraction(i,i+1);
				return true;
			}
		}
	}
	return false;
}

/*--Operation 2--------------*/
bool Operation2(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=1;i<n;i++){
		if(M->h.e[i-1].len()==0 && M->h.e[i].len()>=1){
			for(int j=i;j<n;j++){
				M->Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/*--Elementary transformation avoiding that a short becomes long--------------*/
bool Operation3(CR<Ga3b2>* M,int il,int ir){ // restrict that 1<=il<=i<=ir<=n-1
	for(int i=il;i<=ir;i++) for(int epsilon=-1;epsilon<=1;epsilon+=2){
		Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i,epsilon);
		if(S_complexity(h1)<S_complexity(M->h)){
			M->Elementary_transformation(i,epsilon);
			return true;
		}
	}
	/*To handle cases including (s1, a2ba) and (t1, aba2),
	  we have to consider a sequence of elementary transformtaions of length 2.*/
	for(int i1=il;i1<=ir;i1++) for(int epsilon1=-1;epsilon1<=1;epsilon1+=2)
		for(int i2=il;i2<=ir;i2++) for(int epsilon2=-1;epsilon2<=1;epsilon2+=2){
			Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i1,epsilon1);
			Tuple<Ga3b2> h2=h1  ; h2.Elementary_transformation(i2,epsilon2);
			if(S_complexity(M->h)==S_complexity(h1) && S_complexity(h1)>S_complexity(h2)){
				M->Elementary_transformation(i1,epsilon1);
				M->Elementary_transformation(i2,epsilon2);
				return true;
			}
		}
	return false;
}

/*--Handle the case h_i=Q^{-1}aQ and l(h_{i-1})>l(h_i)--------------*/
bool Operation4(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=2;i<=n;i++){
		Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i-1,1);
		if(S_complexity(h1)==S_complexity(M->h) && h1.e[i-2].len()<M->h.e[i-2].len()){
			M->Elementary_transformation(i-1,1);
			return true;
		}
	}
	return false;
}

/*--find a triple (l,l,l) s.t. l^3=1--------------*/
bool find_triple(CR<Ga3b2>* M,Tuple<Ga3b2>* g_non_inverse_free,list<vector<int>>* F){
	int n=M->h.len();
	for(int i1=1;i1<=n;i1++) for(int i2=i1+1;i2<=n;i2++) for(int i3=i2+1;i3<=n;i3++){
		Ga3b2 g1=M->h.e[i1-1];
		Ga3b2 g2=M->h.e[i2-1];
		Ga3b2 g3=M->h.e[i3-1];
		if(g1==g2 && g1==g3 && pow(g1,3).len()==0){
			for(int j=i3;j<=n-1;j++) M->Elementary_transformation(j,-1);
			for(int j=i2;j<=n-2;j++) M->Elementary_transformation(j,-1);
			for(int j=i1;j<=n-3;j++) M->Elementary_transformation(j,-1);

			g_non_inverse_free->e.insert(g_non_inverse_free->e.begin(),g3);
			g_non_inverse_free->e.insert(g_non_inverse_free->e.begin(),g2);
			g_non_inverse_free->e.insert(g_non_inverse_free->e.begin(),g1);
			M->h.e.pop_back(); M->h.e.pop_back(); M->h.e.pop_back();
			M->H.e.pop_back(); M->H.e.pop_back(); M->H.e.pop_back();

			F->insert(F->end(),M->F.begin(),M->F.end());
			M->F.clear();
			return true;
		}
	}
	return false;
}

/*--find a pair (g1,g2) of inverse elements--------------*/
bool find_pair(CR<Ga3b2>* M,Tuple<Ga3b2>* g_non_inverse_free,list<vector<int>>* F){
	int n=M->h.len();
	for(int i1=1;i1<=n;i1++) for(int i2=i1+1;i2<=n;i2++){
		Ga3b2 g1=M->h.e[i1-1],g2=M->h.e[i2-1];
		if((g1*g2).len()==0){
			for(int j=i2;j<=n-1;j++) M->Elementary_transformation(j,-1);
			for(int j=i1;j<=n-2;j++) M->Elementary_transformation(j,-1);

			g_non_inverse_free->e.insert(g_non_inverse_free->e.begin(),g2);
			g_non_inverse_free->e.insert(g_non_inverse_free->e.begin(),g1);
			M->h.e.pop_back(); M->h.e.pop_back();
			M->H.e.pop_back(); M->H.e.pop_back();

			F->insert(F->end(),M->F.begin(),M->F.end());
			M->F.clear();
			return true;
		}
	}
	return false;
}

/* Restore a short element to a pair of short elements */
void careful_restorations(CR<Ga3b2>* M){
	for(;;){
		int k=M->Restoration();
		if(k==-1) break;
		while(Operation3(M,k,k)) continue;
	}
}

/*--an induction that--------------*/
/*----transform (g1,...,gn) in Ga3b2 whose components are conjuagates of short elements--*/
/*----into-(h1,...,hm) bullet g_non_inverse_free--------------------------------------------*/
/*--s.t. each component of (h1,...,hm) is short--*/
/*--return mp(h,F)--*/
pair<Tuple<Ga3b2>,list<vector<int>>> shorten_induction(Tuple<Ga3b2> g_input,bool details,string namestr){
	cout<<"| +====start the shorten_induction on "+namestr+" (use h to denote the tuple in this process)===\n";
	cout<<"| | h="<<g_input<<"\n";

	CR<Ga3b2> M;
	Tuple<Ga3b2> g_non_inverse_free;
	list<vector<int>> F;

	M.init(g_input);
	g_non_inverse_free.clear();
	F.clear();

	for(int t=1;;t++){
		if(details) cout<<"| | ===the "<<t<<"-th induction on "<<namestr<<" starts===\n";
		if(details) cout<<"| | h="<<M.h<<"\n";
		if(details) cout<<"| | H="<<M.H<<"\n";
		if(details) cout<<"| | g_non_inverse_free="<<g_non_inverse_free<<"\n";

		myassert(M.F.empty(),"F is empty");

		/***********/

		for(;;){
			if(Operation1(&M))              continue;
			if(Operation2(&M))              continue;
			if(Operation3(&M,1,M.h.len()-1)) continue;
			if(Operation4(&M))              continue;
			break;
		}
		if(details) cout<<"| | ---the "<<t<<"-th induction is done---\n";
		if(details) cout<<"| | h="<<M.h<<"\n";
		if(details) cout<<"| | H="<<M.H<<"\n";

		/***********/
		
		careful_restorations(&M);
		if(details) cout<<"| | ---restoration is done---\n";
		if(details) cout<<"| | h="<<M.h<<"\n";
		if(details) cout<<"| | H="<<M.H<<"\n";

		/***********/

		if(find_pair(&M,&g_non_inverse_free,&F))   continue;
		if(find_triple(&M,&g_non_inverse_free,&F)) continue;

		/***********/
		// Fron now on, each component of M.h is short;
		myassert(each_component_is_short(M.h),"fron now on, each component of M.h is short");
		// there are at most 2 a's, at most 2 a's, at most 1 b;
		// a and a^2 cannot occur together; si and ti cannot occur together.

		break;
	}
	F.insert(F.end(),M.F.begin(),M.F.end()); M.F.clear();

	cout<<"| | ===conclusion after shorten_induction===\n";
	cout<<"| | h="<<M.h<<"\n";
	cout<<"| | g_non_inverse_free="<<g_non_inverse_free<<"\n";

	/***********/

	if(details) cout<<"| | ---check elementary transformations in F---\n";
	if(details) cout<<"| | F.size()="<<F.size()<<"\n";
	Tuple<Ga3b2> g=g_input;
	for(auto it:F){
		int i=it[2],epsilon=it[3];
		g.Elementary_transformation(i,epsilon);
	}
	myassert(g==bullet(M.h,g_non_inverse_free),"g==h dot g_non_inverse_free");
	cout<<"| +==="+namestr+" is transformed into  \"h bullet g_non_inverse_free\"  by "<<F.size()<<" elementary transformations\n";
	// return something
	return make_pair(M.h,F);
}

/***************************/
void get_cnts(Tuple<Ga3b2> g,int* cnt_a,int* cnt_a2,int* cnt_b,int* cnt_s,int* cnt_t){
	*cnt_a=0,*cnt_a2=0,*cnt_b=0,*cnt_s=0,*cnt_t=0;
	for(Ga3b2 x:g.e){
		if(x==a)  ++*cnt_a;
		if(x==a2) ++*cnt_a2;
		if(x==b)  ++*cnt_b;
		if(x==s0||x==s1||x==s2) ++*cnt_s;
		if(x==t0||x==t1||x==t2) ++*cnt_t;
	}
}

vector<int> get_positions(Tuple<Ga3b2> g,Ga3b2 x){
	vector<int> poss; poss.clear();
	int n=g.len();
	for(int i=1;i<=n;i++)
		if(g.e[i-1]==x)
			poss.push_back(i);
	return poss;
}

void cyclic_permutation(CR<Ga3b2>* M){ // transform (h1,...,hn) into (h2,...,hn,h1); transform (H1,...,Hn) into (H2,...,Hn,H1).
	int n=M->h.len();
	for(int i=1;i<=n-1;i++) M->Elementary_transformation(i,1);
}

/* Suppose that h1***h_{n2}=1.
   We transform (h1,...,h_{n2},...) into (h_{n2},h1,...,h_{n2-1},...);
      transform (H1,...,H_{n2},...) into (H_{n2},H1,...,H_{n2-1},...).
}
*/
void cyclic_permutation_inv(CR<Ga3b2>* M,int n2){
	Ga3b2 prod; prod.be_identity();
	for(int i=1;i<=n2;i++) prod=prod*M->h.e[i-1];
	myassert(prod.len()==0,"prod=1");
	
	for(int i=n2-1;i>=1;i--) M->Elementary_transformation(i,-1);
}

void cyclic_permutation_inv(CR<Ga3b2>* M){ // transform (h1,...,hn) into (hn,h1,...,h_{n-1}); transform (H1,...,Hn) into (Hn,H1,...,H_{n-1}).
	int n=M->h.len();
	for(int i=n-1;i>=1;i--) M->Elementary_transformation(i,-1);
}

/******************/

namespace tuple_of_short_elements{
	/**investigate a tuple of short elements such that**/
	/**cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0)**/

	bool eliminate_a_and_a(CR<Ga3b2>* M,bool details){ // (a,a)->a^2 or (a^2,a^2)->a
		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(M->h,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
		myassert(cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0),"cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0)");
		
		if(cnt_a==2 || cnt_a2==2){
			Ga3b2 x=(cnt_a==2?a:a2);
			auto poss=get_positions(M->h,x);
			for(int j=poss[0]-1;j>=1;j--) M->Elementary_transformation(j,1);
			for(int j=poss[1]-1;j>=2;j--) M->Elementary_transformation(j,1);
			
			cout<<"| ---combine "<<M->h.e[0]<<" and "<<M->h.e[1]<<"---\n";
			M->Contraction(1,2);
			cout<<"| h="<<M->h<<"\n";
			if(details) cout<<"H="<<M->H<<"\n";
			return true;
		}
		return false;
	}

	bool eliminate_a_and_b(CR<Ga3b2>* M,bool details){ // I_a=1=I_b
		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(M->h,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
		myassert(cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0),"cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0)");
		
		if(cnt_a+cnt_a2==1 && cnt_b==1){ 
			Ga3b2 x=(cnt_a==1?a:a2);
			for(;M->h.e[0]!=x;) cyclic_permutation(M); // h1=x
			cout<<"| ---cyclic_permutation--\n";
			cout<<"| h="<<M->h<<"\n";
			for(int i=2;i<=M->h.len();i++) // attemmpt to make (a,a^2b)->b, (a,ba)->aba, (a,a^2ba^2)->ba^2 or (a,ab)->a^2b
				if(is_short(M->h.e[i-1])>=3 && short_table[is_short(M->h.e[0])][is_short(M->h.e[i-1])]){
					for(int j=i-1;j>=2;j--) M->Elementary_transformation(j,1); // (a,ba,...)

					cout<<"| ---combine "<<M->h.e[0]<<" and "<<M->h.e[1]<<"---\n";
					M->Contraction(1,2);
					cout<<"| h="<<M->h<<"\n";
					if(details) cout<<"| H="<<M->H<<"\n";
					return true;
				}
			cyclic_permutation(M); // hn=x
			cout<<"| ---cyclic_permutation--\n";
			cout<<"| h="<<M->h<<"\n";
			for(int i=1;i<=M->h.len()-1;i++) // attemmpt to make (ba,a)->ba^2
				if(is_short(M->h.e[i-1])>=3 && short_table[is_short(M->h.e[i-1])][is_short(M->h.e[M->h.len()-1])]){
					for(int j=i;j<=M->h.len()-1;j++) M->Elementary_transformation(j,-1); // (...,ba,a)
					
					cout<<"| ---combine "<<M->h.e[M->h.len()-2]<<" and "<<M->h.e[M->h.len()-1]<<"---\n";
					M->Contraction(M->h.len()-1,M->h.len());
					cout<<"| h="<<M->h<<"\n";
					if(details) cout<<"| H="<<M->H<<"\n";
					return true;
				}
			// Otherwise, A={a^e,b,a^eba^e}.
			myassert(M->h.e[M->h.len()-1]==x,"hn=x");
			for(int i=M->h.len()-1;i>=1;i--) // attemmpt to make (aba,a,aba)->(aba,a^2b,a)->(a,a)->a^2
				if((x==a && M->h.e[i-1]==s1)||(x==a2 && M->h.e[i-1]==t1)){
					for(int j=i;j<=M->h.len()-2;j++) M->Elementary_transformation(j,-1);
					break;
				}
			for(auto poss=get_positions(M->h,x);poss[0]!=2;poss=get_positions(M->h,x)) cyclic_permutation(M); // h2=x
			for(int i=3;i<=M->h.len();i++)
				if((x==a && M->h.e[i-1]==s1)||(x==a2 && M->h.e[i-1]==t1)){
					for(int j=i-1;j>=3;j--) M->Elementary_transformation(j,1);
					break;
				}
			if(x==a)  myassert(M->h.e[0]==s1 && M->h.e[1]==a  && M->h.e[2]==s1,"(aba,a,aba)");
			if(x==a2) myassert(M->h.e[0]==t1 && M->h.e[1]==a2 && M->h.e[2]==s2,"(aabaa,aa,aabaa)");
			M->Elementary_transformation(2,1); M->Contraction(1,2); M->Contraction(1,2);
			
			cout<<"| ---transform and contract (a^eba^e,a^e,a^eba^e) into a^2---\n";
			cout<<"| h="<<M->h<<"\n";
			if(details) cout<<"| H="<<M->H<<"\n";
			return true;
		}
		return false;
	}

	/*--Suppose that h1,...,hn are short s.t. at most 1 of them is equal to one of a,a^2,b--*/
	/*--Rearrange them and transform h into a tuple of the form either--*/
	/*----(a/a^2/b,s---s) * (s_0,t_0)^k or (a/a^2/b,t---t) * (s_0,t_0)^k--*/
	int Rearrangement(CR<Ga3b2>* M,bool details){
		int n=M->h.len();
		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(M->h,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
		myassert(cnt_a+cnt_a2+cnt_b<=1,"at most 1 of h1,...,hn is equal to one of a,a^2,b");

		if(cnt_a+cnt_a2+cnt_b==1)
			for(;!(M->h.e[0]==a || M->h.e[0]==a2 || M->h.e[0]==b);) cyclic_permutation(M); // h1=a/a^2/b
		for(;;){
			bool ok=true;
			for(int i=1;ok && i<=n-1;i++){ // avoid the appearance of (ti,sj)
				if((M->h.e[i-1]==t0 || M->h.e[i-1]==t1 || M->h.e[i-1]==t2)&&(M->h.e[i]==s0 || M->h.e[i]==s1 || M->h.e[i]==s2)){
					if(is_short(M->h.e[i])+3==is_short(M->h.e[i-1])){ // (t0,s0), (t1,s1) or (s2,s2)
						for(int j=i+1;j<=n-1;j++) M->Elementary_transformation(j,-1);
						for(int j=i  ;j<=n-2;j++) M->Elementary_transformation(j,-1);
						return n-2;
					}
					if((is_short(M->h.e[i])-3+2)%3==(is_short(M->h.e[i-1])-6)%3){ // (t_{k+1},s_{k-1}) -> (s_k, t_{k+1})
						M->Elementary_transformation(i,-1);
						ok=false; break;
					}
					if((is_short(M->h.e[i])-3+1)%3==(is_short(M->h.e[i-1])-6)%3){ // (t_{k+1},s_{k}) -> (s_k, t_{k-1})
						M->Elementary_transformation(i,1);
						ok=false; break;
					}
				}
			}
			if(ok){
				cout<<"| ---Rearrange (h1,...,hn) of short elements at most 1 of which is equal to a/a^2/b---\n";
				cout<<"| h="<<M->h<<"\n";
				// s_i,...,t_i
				for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++){
					if((M->h.e[i-1]==s0 && M->h.e[j-1]==t0)||
					   (M->h.e[i-1]==s1 && M->h.e[j-1]==t1)||
					   (M->h.e[i-1]==s2 && M->h.e[j-1]==t2)){
						for(int k=j;k<=n-1;k++) M->Elementary_transformation(k,-1);
						for(int k=i;k<=n-2;k++) M->Elementary_transformation(k,-1);
						return n-2;
					}
				}
				// s_i,...,t_{i+1},...,t_{i+2}
				for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++) for(int k=j+1;k<=n;k++){
					if((M->h.e[i-1]==s0 && M->h.e[j-1]==t1 && M->h.e[k-1]==t2)||
					   (M->h.e[i-1]==s1 && M->h.e[j-1]==t2 && M->h.e[k-1]==t0)||
					   (M->h.e[i-1]==s2 && M->h.e[j-1]==t0 && M->h.e[k-1]==t1)){
						for(int l=k;l<=n-1;l++) M->Elementary_transformation(l,-1);
						for(int l=j;l<=n-2;l++) M->Elementary_transformation(l,-1);
						for(int l=i;l<=n-3;l++) M->Elementary_transformation(l,-1);
						M->Elementary_transformation(n-1,1); // (s_i,t_{i+2},t_i)
						M->Elementary_transformation(n-2,-1); // (*,s_i,t_i)
						return n-2;
					}
				}
				// s_{i+2},...,s_{i+1},...,t_{i}
				for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++) for(int k=j+1;k<=n;k++){
					if((M->h.e[i-1]==s2 && M->h.e[j-1]==s1 && M->h.e[k-1]==t0)||
					   (M->h.e[i-1]==s0 && M->h.e[j-1]==s2 && M->h.e[k-1]==t1)||
					   (M->h.e[i-1]==s1 && M->h.e[j-1]==s0 && M->h.e[k-1]==t2)){
						for(int l=k;l<=n-1;l++) M->Elementary_transformation(l,-1);
						for(int l=j;l<=n-2;l++) M->Elementary_transformation(l,-1);
						for(int l=i;l<=n-3;l++) M->Elementary_transformation(l,-1);
						M->Elementary_transformation(n-2,1); // (s_{i+1},s_i,t_i)
						return n-2;
					}
				}
				// a,...,s_{i+1},...,t_i
				for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++){
					if(M->h.e[0]!=a) continue;
					if((M->h.e[i-1]==s1 && M->h.e[j-1]==t0)||
					   (M->h.e[i-1]==s2 && M->h.e[j-1]==t1)||
					   (M->h.e[i-1]==s0 && M->h.e[j-1]==t2)){
						for(int l=j;l<=n-1;l++) M->Elementary_transformation(l,-1);
						for(int l=i;l<=n-2;l++) M->Elementary_transformation(l,-1);
						for(int l=1;l<=n-3;l++) M->Elementary_transformation(l,-1);
						M->Elementary_transformation(n-2,-1); // (s_i,a,t_i)
						M->Elementary_transformation(n-2,-1); // (*,s_i,t_i)
						return n-2;
					}
				}
				// a^2,...,s_i,...,t_{i+1}
				for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++){
					if(M->h.e[0]!=a2) continue;
					if((M->h.e[i-1]==s0 && M->h.e[j-1]==t1)||
					   (M->h.e[i-1]==s1 && M->h.e[j-1]==t2)||
					   (M->h.e[i-1]==s2 && M->h.e[j-1]==t0)){
						for(int l=j;l<=n-1;l++) M->Elementary_transformation(l,-1);
						for(int l=i;l<=n-2;l++) M->Elementary_transformation(l,-1);
						for(int l=1;l<=n-3;l++) M->Elementary_transformation(l,-1);
						M->Elementary_transformation(n-2,-1); // (s_{i+1},a^2,t_{i+1})
						M->Elementary_transformation(n-2,-1); // (*,s_{i+1},t_{i+1})
						return n-2;
					}
				}
				// b,...,s0,...,t2 or b,...,s2,...,t0
				for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++){
					if(M->h.e[0]!=b) continue;
					if((M->h.e[i-1]==s0 && M->h.e[j-1]==t2)||(M->h.e[i-1]==s2 && M->h.e[j-1]==t0)){
						for(int l=j;l<=n-1;l++) M->Elementary_transformation(l,-1);
						for(int l=i;l<=n-2;l++) M->Elementary_transformation(l,-1);
						for(int l=1;l<=n-3;l++) M->Elementary_transformation(l,-1);
						M->Elementary_transformation(n-2,-1); // (s2,b,t2) or (s0,b,t0)
						M->Elementary_transformation(n-2,-1); // (*,s_i,t_i)
						return n-2;
					}
				}

				set<int> As,At; As.clear(); At.clear();
				for(int i=1;i<=n;i++){
					if     (M->h.e[i-1]==s0 || M->h.e[i-1]==s1 || M->h.e[i-1]==s2) As.insert(is_short(M->h.e[i-1]));
					else if(M->h.e[i-1]==t0 || M->h.e[i-1]==t1 || M->h.e[i-1]==t2) At.insert(is_short(M->h.e[i-1]));
				}
				if(As.size()==0 || At.size()==0) return n;

				myassert((As.size()==1 && At.size()==2)||(As.size()==2 && At.size()==1),"either |As|==1, |At|==2 or |As|==2, |At|==1");
				myassert(M->h.e[0]==b,"the tuple must start with b");
				if(As.size()==1 && At.size()==2){
					// b,s1,...,s1,t0,...,t0,t2,...,t2
					int u=0; while(1+    (u+1)<=n && M->h.e[1+    (u+1)-1]==s1) ++u;
					int v=0; while(1+u+  (v+1)<=n && M->h.e[1+u+  (v+1)-1]==t0) ++v;
					int w=0; while(1+u+v+(w+1)<=n && M->h.e[1+u+v+(w+1)-1]==t2) ++w;
					myassert(1+u+v+w==n && u>0 && v>0 && w>0,"the tuple must be of the form (b,s1,...,s1,t0,...,t0,t2,...,t2");

					for(int i=1;i<=v+w;i++) cyclic_permutation_inv(M); // t0,...,t0,t2,...,t2,b,s1,...,s1
					for(int i=v+w;i>=1;i--) M->Elementary_transformation(i,1); // b,t2,...,t2,t0,...t0,s1,...s1
					// t_{i+1},...,t_{i+2},...,s_i
					for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++) for(int k=j+1;k<=n;k++){
						if((M->h.e[i-1]==t2 && M->h.e[j-1]==t1 && M->h.e[k-1]==s0)||
						   (M->h.e[i-1]==t0 && M->h.e[j-1]==t2 && M->h.e[k-1]==s1)||
						   (M->h.e[i-1]==t1 && M->h.e[j-1]==t0 && M->h.e[k-1]==s2)){
							for(int l=k;l<=n-1;l++) M->Elementary_transformation(l,-1);
							for(int l=j;l<=n-2;l++) M->Elementary_transformation(l,-1);
							for(int l=i;l<=n-3;l++) M->Elementary_transformation(l,-1);
							M->Elementary_transformation(n-2,1); // (t_{i+2},t_i,s_i)
							return n-2;
						}
					}
					myassert(true,"ERROR");
				}else if(As.size()==2 && At.size()==1){
					// b,s2,...,s2,s0,...,s0,t1,...,t1
					int u=0; while(1+    (u+1)<=n && M->h.e[1+    (u+1)-1]==s2) ++u;
					int v=0; while(1+u+  (v+1)<=n && M->h.e[1+u+  (v+1)-1]==s0) ++v;
					int w=0; while(1+u+v+(w+1)<=n && M->h.e[1+u+v+(w+1)-1]==t1) ++w;
					myassert(1+u+v+w==n && u>0 && v>0 && w>0,"the tuple must be of the form (b,s2,...,s2,s0,...,s0,t1,...,t1");

					for(int i=1;i<=1+u+v;i++) cyclic_permutation(M); // t1,...,t1,b,s2,...,s2,s0,...,s0
					for(int i=w+1;i<n;i++) M->Elementary_transformation(i,-1); // t1,...,t1,s0,...,s0,s2,...,s2,b
					// t_{i},...,s_{i+2},...,s_{i+1}
					for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++) for(int k=j+1;k<=n;k++){
						if((M->h.e[i-1]==t0 && M->h.e[j-1]==s2 && M->h.e[k-1]==s1)||
						   (M->h.e[i-1]==t1 && M->h.e[j-1]==s0 && M->h.e[k-1]==s2)||
						   (M->h.e[i-1]==t2 && M->h.e[j-1]==s1 && M->h.e[k-1]==s0)){
							for(int l=k;l<=n-1;l++) M->Elementary_transformation(l,-1);
							for(int l=j;l<=n-2;l++) M->Elementary_transformation(l,-1);
							for(int l=i;l<=n-3;l++) M->Elementary_transformation(l,-1);
							M->Elementary_transformation(n-1,1); // (t_i,s_{i+1},s_i)
							M->Elementary_transformation(n-2,-1); // (*,t_i,s_i)
							return n-2;
						}
					}
					myassert(true,"ERROR");
				}
			}
		}
	}
}

/*--transform/contract--*/
/*----a tuple in Ga3b2 whose components are conjuagates of short elements--------------*/
/*--into an inverse-free tuple--*/
pair<Tuple<Ga3b2>,list<vector<int>>> transform_into_inverse_free(Tuple<Ga3b2> g_input,bool details){
	cout<<"+====transform a tuple into an inverse-free tuple===\n";
	CR<Ga3b2> MM; MM.init(g_input);
	Tuple<Ga3b2> g; g=g_input; // g must be a prefix of M.h

	cout<<"| g="<<g<<"\n";
	for(;;){
		auto ret=shorten_induction(g,details,"h");
		MM.Apply(ret.second); g=ret.first;
		myassert(each_component_is_short(g),"after the shorten_induction, each component must be short");
		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(g,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
		myassert(cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0),"after the shorten_induction, cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0)");

		cout<<"| ===get a tuple of short elements==\n";
		cout<<"| h="<<g<<"\n";

		CR<Ga3b2> M; M.init(g);
		if(tuple_of_short_elements::eliminate_a_and_a(&M,details)){ // (a,a)->a^2 or (a^2,a^2)->a
			MM.Apply(M.F); g=M.h; continue;
		}
		if(tuple_of_short_elements::eliminate_a_and_b(&M,details)){ // I_a=1=I_b
			MM.Apply(M.F); g=M.h; continue;
		}
		int newn=tuple_of_short_elements::Rearrangement(&M,details);
		{
			MM.Apply(M.F); g.init(vector<Ga3b2>(M.h.e.begin(),M.h.e.begin()+newn));
			if(newn<M.h.len()) continue;
			else break;
		}
	}

	cout<<"| ===conclusion===\n";
	cout<<"| g="<<MM.h<<"\n";
	cout<<"| G="<<MM.H<<"\n";
	cout<<"| F.size()="<<MM.F.size()<<"\n";
	cout<<"+====\n";
	return make_pair(g,MM.F);
}

/******************/

namespace Moishezon{

	/* normalize_s is an extension of the procedure introduced by Moishezon

	   Given a tuple (g1,...,gn) of {a,b,s0,s1,s2} that has at most 1 component being a or b,
	   it is transformed into (h1,...,hn) s.t. either
	   (1) (h1,...,hn) = (s0,s1,s0,s1,s0,s1)^m, or
	   (2) (h1,h2) = (a,s0) or (hn,h1) = (s2,a) or (hn,h1,h2) = (s1,a,s1), or
	   (3) h starts with (b,s2) or (b,s0,...,s0,s2,s0), or
	       (the cyclic_permutation of) h ends with (s0,b) or (s2,s0,s2,...,s2,b).
	*/
	void normalize_s(CR<Ga3b2>* M,bool details){ // by Moishezon
		cout<<"+====normalize \"a tuple of {s0,s1,s2}\" (in fact, a tuple of {a,b,s0,s1,s2}) by Moishezon===\n";
		cout<<"| g="<<M->h<<"\n";

		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(M->h,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
		while(cnt_a+cnt_b==1 && M->h.e[0]!=a && M->h.e[0]!=b) cyclic_permutation(M);

		Tuple<Ga3b2> g=M->h; // g must be a prefix of M.h, g starts with a or b when cnt_a+cnt_b==1
		cout<<"| ->"<<g<<"\n";

		for(;;){
			if(details) cout<<"| g="<<M->h<<"\n";
			if(details) cout<<"| h="<<g<<"\n";
			bool ok=true;
			// (s1,s0)=(aba,a^2b)->(a^2b,ba aba a^2b)=(a^2b,ba^2)=(s0,s2), which decreases # of s1
			for(int i=1;ok && i+1<=g.len();i++)
				if(g.e[i-1]==s1 && g.e[i]==s0){
					M->Elementary_transformation(i,1);
					g.Elementary_transformation(i,1);
					ok=false;
				}
			// (s2,s1)=(ba^2,aba)->(ba^2 aba ab,ba^2)=(a^2b,ba^2)=(s0,s2), which decreases # of s1
			for(int i=1;ok && i+1<=g.len();i++)
				if(g.e[i-1]==s2 && g.e[i]==s1){
					M->Elementary_transformation(i,-1);
					g.Elementary_transformation(i,-1);
					ok=false;
				}
			// (s2,s0,s2)->(s2,s1,s0)->(s0,s2,s0), which decreases the lexicographical order
			for(int i=1;ok && i+2<=g.len();i++)
				if(g.e[i-1]==s2 && g.e[i]==s0 && g.e[i+1]==s2){
					M->Elementary_transformation(i+1,-1); M->Elementary_transformation(i,-1);
					g.Elementary_transformation(i+1,-1); g.Elementary_transformation(i,-1);
					ok=false;
				}
			// (s2,s0,s0,s2,s0,s0)->(s2,s0,s2,s0,s2,s0)->(s0,s2,s0,s2,s0,s2)
			for(int i=1;ok && i+5<=g.len();i++){
				if(g.e[i-1]==s2 && g.e[i]==s0 && g.e[i+1]==s0 && g.e[i+2]==s2 && g.e[i+3]==s0 && g.e[i+4]==s0){
					M->Elementary_transformation(i+3,-1); M->Elementary_transformation(i+2,-1);
					g.Elementary_transformation(i+3,-1); g.Elementary_transformation(i+2,-1);
					for(int j=i;j<=i+4;j++){
						M->Elementary_transformation(j,1);
						g.Elementary_transformation(j,1);
					}
					for(int j=i+5;j>=i;j--)
						for(int k=j;k+(i+6-j)<=g.len();k++){
							M->Elementary_transformation(k,-1);
							g.Elementary_transformation(k,-1);
						}
					for(int j=1;j<=6;j++) g.e.pop_back();
					ok=false;
				}
			}
			// (s0,s0,s2,s0,s0,s2)->(s0,s2,s1,s0,s0,s2)->(s0,s2,s0,s2,s0,s2)
			for(int i=1;ok && i+5<=g.len();i++){
				if(g.e[i-1]==s0 && g.e[i]==s0 && g.e[i+1]==s2 && g.e[i+2]==s0 && g.e[i+3]==s0 && g.e[i+4]==s2){
					M->Elementary_transformation(i+1,1); M->Elementary_transformation(i+2,1);
					g.Elementary_transformation(i+1,1); g.Elementary_transformation(i+2,1);
					for(int j=i+5;j>=i;j--)
						for(int k=j;k+(i+6-j)<=g.len();k++){
							M->Elementary_transformation(k,-1);
							g.Elementary_transformation(k,-1);
						}
					for(int j=1;j<=6;j++) g.e.pop_back();
					ok=false;
				}
			}
			// (s0,s2,s0,...(k>=0)...,s0,s2,s0)->(s0,s2,s0,s2,s0,s2,...(k>=0)...)
			for(int i=1;ok && i+5<=g.len();i++){
				if(g.e[i-1]==s0 && g.e[i]==s2 && g.e[i+1]==s0 && g.e[i+2]==s0){
					int k=0;
					while(i+2+(k+1)+1<=g.len() && g.e[i+2+(k+1)]==s0) ++k;
					if(i+2+(k+2)+1<=g.len() && g.e[i+2+(k+1)]==s2 && g.e[i+2+(k+2)]==s0){
						for(int j=i+2+(k+2)+1;j>=i+5;j--){
							M->Elementary_transformation(j-2,1); M->Elementary_transformation(j-1,1);
							g.Elementary_transformation(j-2,1); g.Elementary_transformation(j-1,1);
						}
						for(int jj=i+5;jj>=i;jj--)
							for(int kk=jj;kk+(i+6-jj)<=g.len();kk++){
								M->Elementary_transformation(kk,-1);
								g.Elementary_transformation(kk,-1);
							}
						for(int jj=1;jj<=6;jj++) g.e.pop_back();
						ok=false;
					}
				}
			}
			if(ok) break;
		}
		cout<<"|----conclusion---\n";
		cout<<"| g="<<M->h<<"\n";
		cout<<"| h="<<g<<"\n";
		
		// --- (h1,h2) = (a,s0) or (hn,h1) = (s2,a) or (hn,h1,h2) = (s1,a,s1), or
		// --- h starts with (b,s2) or (b,s0,...,s0,s2,s0), or
		// --- (the cyclic_permutation of) h ends with (s0,b) or (s2,s0,s2,...,s2,b).

		// Notice: g.len <= M.h.len()
		if(g.len()>0){
			if(g.e[0]==a){
				if(g.e[1]==s0) // (h1,h2) = (a,s0) ---> (b)
					M->Contraction(1,2);
				else if(g.e[g.len()-1]==s2){ // (hn,h1) = (s2,a) ---> (b)
					cyclic_permutation_inv(M, g.len()); // (h1,h2) = (s2,a)
					M->Contraction(1,2); // ---> (b)
				}else if(g.e[g.len()-1]==s1 && g.e[1]==s1){
					cyclic_permutation_inv(M, g.len()); // (h1,h2,h3) = (s1,a,s1)
					M->Elementary_transformation(2,-1); // (h1,h2,h3) = (s1,s0,a)
					M->Contraction(1,2); // ---> (a,a)
					M->Contraction(1,2); // ---> (a)
				}else
					myassert(true,"when h1=a, there is no more possibility");
			}else if(g.e[0]==b){
				if(g.e[1]==s2) // (h1,h2) = (b,s2) ---> (a^2)
					M->Contraction(1,2);
				else if(g.e[g.len()-1]==s0){ // (hn,h1) = (s0,b)
					cyclic_permutation_inv(M, g.len()); // (h1,h2) = (s0,b)
					M->Contraction(1,2); // ---> (a^2)
				}else{
					int n2=g.len(); bool chk=false;
					if(g.e[1]==s0){
						int k=1;
						while(1+(k+1)<=n2 && g.e[k+1]==s0) ++k;
						if(1+(k+2)<=n2 && g.e[k+1]==s2 && g.e[k+2]==s0){ // (h1,h2,...,h_{k+3}) = (b,s0,...,s0,s2,s0)
							for(int i=k+1;i>=2;i--){ // (h_i,h_{i+1},h_{i+2}) = (s0,s2,s0) -> (s2,s1,s0) -> (s2,s0,s2)
								M->Elementary_transformation(i,1); M->Elementary_transformation(i+1,1);
							}
							M->Contraction(1,2); // (h1,h2) = (b,s2) ---> (a^2)
							chk=true;
						}
					}
					if(!chk && g.e[n2-1]==s2){
						int k=1;
						while(1+(k+1)<=n2 && g.e[n2-(k+1)]==s2) ++k;
						if(1+(k+2)<=n2 && g.e[n2-(k+1)]==s0 && g.e[n2-(k+2)]==s2){ // (h1,...,h_{n2},...) = (b,...,s2,s0,s2,...,...)
							for(int i=k+1;i>=2;i--){ // (h_{n2-i},h_{n2-i+1},h_{n2-i+2}) = (s2,s0,s2)
								M->Elementary_transformation(n2-i,1); M->Elementary_transformation(n2-i+1,1);
							}
							cyclic_permutation_inv(M, n2); // (h1,h2) = (s0,b)
							M->Contraction(1,2); // ---> (a^2)
							chk=true;
						}
					}
					myassert(chk,"when h1=b, there is no more possibility");
				}
			}

			cout<<"|----contraction---\n";
			cout<<"| g="<<M->h<<"\n";
			cout<<"| G="<<M->H<<"\n";
		}
		cout<<"+=======\n";
	}

	void normalize_t(CR<Ga3b2>* M,bool details){
		vector<Ga3b2> g; g.clear();
		for(Ga3b2 x:M->h.e) g.push_back(x.inv());
		reverse(g.begin(),g.end());
		CR<Ga3b2> M2; M2.init(Tuple<Ga3b2>(g));
		normalize_s(&M2,details);

		// R_1:     (g,h)        ---> (h,h^{-1}gh)
		// R_1^{-1}:(h^{-1},g^{-1})  ---> (h^{-1}g^{-1}h,h^{-1})
		// R_i ~~~> R_{n-i}^{-1} and R_i^{-1} ~~~> R_{n-i}
		
		cout<<"+====normalize \"a tuple of {t0,t1,t2}\" (in fact, a tuple of {a,b,s0,s1,s2}) by Moishezon===\n";
		cout<<"| g="<<M->h<<"\n";
		for(auto it:M2.F){
			if(it[0]==1){
				int n=M->h.len();
				M->Elementary_transformation(n-it[2],-it[3]);
			}else{
				int n=M->h.len();
				M->Contraction(n+1-it[2],n+1-it[1]);
			}

		}
		cyclic_permutation(M); // (ab,ba,ab,ba,...) -> (ba,ab,ba,ab,...)
		cout<<"| ->"<<M->h<<"\n";
		cout<<"+=======\n";
	}
}

/*****************/

void combine_a_t(CR<Ga3b2>* M,bool details){
	while(M->h.e[0]!=a) cyclic_permutation(M);
	cout<<"---combine "<<M->h.e[0]<<" and some t_i="<<M->h.e[1]<<"---\n";
	M->Contraction(1,2);
}

void combine_a2_s(CR<Ga3b2>* M,bool details){
	while(M->h.e[0]!=a2) cyclic_permutation(M);
	cout<<"---combine "<<M->h.e[0]<<" and some t_i="<<M->h.e[1]<<"---\n";
	M->Contraction(1,2);
}

list<vector<int>> normalize_inverse_free_tuple(Tuple<Ga3b2> h_input,bool details){
	cout<<"===normalize an inverse-free tuple===\n";
	CR<Ga3b2> MM; MM.init(h_input);

	for(Tuple<Ga3b2> g=h_input;;){ // g must be a prefix of MM.h
		cout<<"g="<<g<<"\n";

		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(g,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
		myassert(each_component_is_short(g),"the tuple is of short elements");
		myassert(cnt_a+cnt_a2+cnt_b<=1 && (cnt_s==0 || cnt_t==0),"the tuple is inverse-free and contains at most 1 a/a^2/b");

		CR<Ga3b2> M; M.init(g);
		bool cont=false;

		if(!cont && cnt_a==1 && cnt_s==0){ // (a,t,...) ---> (t,...,t) * (a,ti,tj)
			combine_a_t(&M,details);
			auto ret=transform_into_inverse_free(M.h,details);
			M.Apply(ret.second); g=ret.first; cont=true;
		}

		if(!cont && cnt_a2==1 && cnt_t==0){ // (a^2,s,...) ---> (s,...,s) * (a^2,si,sj)
			combine_a2_s(&M,details);
			auto ret=transform_into_inverse_free(M.h,details);
			M.Apply(ret.second); g=ret.first; cont=true;
		}

		if(!cont && (cnt_a==1 || cnt_b==1) && cnt_t==0){ // (a,s,...) or (b,s,...)
			Moishezon::normalize_s(&M,details);
			g=M.h; cont=true;
		}

		if(!cont && (cnt_a2==1 || cnt_b==1) && cnt_s==0){ // (a^2,t,...) or (b,t,...)
			Moishezon::normalize_t(&M,details);
			g=M.h; cont=true;
		}

		if(!cont && cnt_a+cnt_a2+cnt_b==0){
			if(cnt_t==0){ // a tuple of s0,s1,s2
				Moishezon::normalize_s(&M,details);
			}else{ // a tuple of t0,t1,t2
				Moishezon::normalize_t(&M,details);
			}
			g=M.h; // no need to continue anymore
		}

		MM.Apply(M.F);
		if(!cont) break;
	}

	return MM.F;
}


/***gotcha begin***/

/**find (Q^{-1}s_1Q,Q^{-1}t_1Q) and (Q^{-1}t_1Q,Q^{-1}s_1Q)**/
bool gotcha_s_and_t(CR<Ga3b2>* M,int nh){
	for(int i=1,j;i<=nh;i=j+1){
		j=i; Ga3b2 prod=M->h.e[i-1];
		while(prod.len()>0) prod=prod*M->h.e[j++];
		if(i+1==j){
			Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i];
			if((g1*g2).len()!=0) continue;
			if(3<=is_short(get_tau_in_expression(g2)) && is_short(get_tau_in_expression(g2))<6){
				M->Elementary_transformation(i,1);
				swap(g1,g2);
			}
			if(3<=is_short(get_tau_in_expression(g1)) && is_short(get_tau_in_expression(g1))<6){
				for(int j=i+1;j<=nh-1;j++) M->Elementary_transformation(j,-1);
				for(int j=i  ;j<=nh-2;j++) M->Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/**find (Q^{-1}aQ,Q^{-1}a^2Q) and (Q^{-1}a^2Q,Q^{-1}aQ)**/
bool gotcha_a_and_a(CR<Ga3b2>* M,int nh){
	for(int i=1,j;i<=nh;i=j+1){
		j=i; Ga3b2 prod=M->h.e[i-1];
		while(prod.len()>0) prod=prod*M->h.e[j++];
		if(i+1==j){
			Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i];
			if(is_short(get_tau_in_expression(g2))==0){
				M->Elementary_transformation(i,1);
				swap(g1,g2);
			}
			if(is_short(get_tau_in_expression(g1))==0){
				for(int j=i+1;j<=nh-1;j++) M->Elementary_transformation(j,-1);
				for(int j=i  ;j<=nh-2;j++) M->Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/**find (Q^{-1}bQ,Q^{-1}bQ)**/
bool gotcha_b_and_b(CR<Ga3b2>* M,int nh){
	for(int i=1,j;i<=nh;i=j+1){
		j=i; Ga3b2 prod=M->h.e[i-1];
		while(prod.len()>0) prod=prod*M->h.e[j++];
		if(i+1==j){
			Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i];
			if(is_short(get_tau_in_expression(g1))==2){
				for(int j=i+1;j<=nh-1;j++) M->Elementary_transformation(j,-1);
				for(int j=i  ;j<=nh-2;j++) M->Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/**find (Q^{-1}aQ,Q^{-1}aQ,Q^{-1}aQ)**/
bool gotcha_a_and_a_and_a(CR<Ga3b2>* M,int nh){
	for(int i=1,j;i<=nh;i=j+1){
		j=i; Ga3b2 prod=M->h.e[i-1];
		while(prod.len()>0) prod=prod*M->h.e[j++];
		if(i+2==j){
			Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i],g3=M->h.e[i+1];
			if(g1==g2 && g2==g3 && is_short(get_tau_in_expression(g1))==0){
				for(int j=i+2;j<=nh-1;j++) M->Elementary_transformation(j,-1);
				for(int j=i+1;j<=nh-2;j++) M->Elementary_transformation(j,-1);
				for(int j=i  ;j<=nh-3;j++) M->Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/**find (Q^{-1}a^2Q,Q^{-1}a^2Q,Q^{-1}a^2Q)**/
bool gotcha_a2_and_a2_and_a2(CR<Ga3b2>* M,int nh){
	for(int i=1,j;i<=nh;i=j+1){
		j=i; Ga3b2 prod=M->h.e[i-1];
		while(prod.len()>0) prod=prod*M->h.e[j++];
		if(i+2==j){
			Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i],g3=M->h.e[i+1];
			if(g1==g2 && g2==g3 && is_short(get_tau_in_expression(g1))==1){
				for(int j=i+2;j<=nh-1;j++) M->Elementary_transformation(j,-1);
				for(int j=i+1;j<=nh-2;j++) M->Elementary_transformation(j,-1);
				for(int j=i  ;j<=nh-3;j++) M->Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/**find (s0,s2,s0,s2,s0,s2)**/
bool gotcha_s020202(CR<Ga3b2>* M,int nh){
	for(int i=1,j;i<=nh;i=j+1){
		j=i; Ga3b2 prod=M->h.e[i-1];
		while(prod.len()>0) prod=prod*M->h.e[j++];
		if(i+5==j){
			Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i],g3=M->h.e[i+1],g4=M->h.e[i+2],g5=M->h.e[i+3],g6=M->h.e[i+4];
			if(g1==s0 && g2==s2 && g3==s0 && g4==s2 && g5==s0 && g6==s2){
				for(int j=i+5;j<=nh-1;j++) M->Elementary_transformation(j,-1);
				for(int j=i+4;j<=nh-2;j++) M->Elementary_transformation(j,-1);
				for(int j=i+3;j<=nh-3;j++) M->Elementary_transformation(j,-1);
				for(int j=i+2;j<=nh-4;j++) M->Elementary_transformation(j,-1);
				for(int j=i+1;j<=nh-5;j++) M->Elementary_transformation(j,-1);
				for(int j=i  ;j<=nh-6;j++) M->Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/**find (t0,t2,t0,t2,t0,t2)**/
bool gotcha_t020202(CR<Ga3b2>* M,int nh){
	for(int i=1,j;i<=nh;i=j+1){
		j=i; Ga3b2 prod=M->h.e[i-1];
		while(prod.len()>0) prod=prod*M->h.e[j++];
		if(i+5==j){
			Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i],g3=M->h.e[i+1],g4=M->h.e[i+2],g5=M->h.e[i+3],g6=M->h.e[i+4];
			if(g1==t0 && g2==t2 && g3==t0 && g4==t2 && g5==t0 && g6==t2){
				for(int j=i+5;j<=nh-1;j++) M->Elementary_transformation(j,-1);
				for(int j=i+4;j<=nh-2;j++) M->Elementary_transformation(j,-1);
				for(int j=i+3;j<=nh-3;j++) M->Elementary_transformation(j,-1);
				for(int j=i+2;j<=nh-4;j++) M->Elementary_transformation(j,-1);
				for(int j=i+1;j<=nh-5;j++) M->Elementary_transformation(j,-1);
				for(int j=i  ;j<=nh-6;j++) M->Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/***gotcha end***/



/***search begin***/

map<vector<int>,pair<vector<int>,list<vector<int>>>> all_tuples;
bool cont_search;

void search_tuples_of_short_elements(int n,Tuple<Ga3b2> g,int m,vector<int> idx_h,vector<int> idx_lastg,list<vector<int>> F,int height,bool details){
	if(!cont_search) return;
	vector<int> idx_g; idx_g.clear(); for(int i=0;i<n;i++) idx_g.push_back(is_short(g.e[i]));

	if(all_tuples.find(idx_g)!=all_tuples.end()) return;
	else all_tuples[idx_g]=make_pair(idx_lastg,F);

	if(vector<int>(idx_g.begin()+n-m,idx_g.begin()+n)==idx_h){
		cont_search=false;
		return;
	}
	if(height>=1000) return; // This is heuristic

	if(details) cout<<"search: "<<g<<"   height="<<height<<"\n";

	Tuple<Ga3b2> g2;
	list<vector<int>> F2;

	for(int i=1;i<=n-1;i++) for(int epsilon=-1;epsilon<=1;epsilon+=2){
		F2.clear(); g2=g; g2.Elementary_transformation(i,epsilon); F2.push_back({1,n,i,epsilon});
		if(each_component_is_short(g2))
			search_tuples_of_short_elements(n,g2,m,idx_h,idx_g,F2,height+1,details);
	}

	for(int ii=1,jj;ii<=n;ii=jj+1){
		jj=ii; Ga3b2 prod=g.e[jj-1];
		while(prod.len()>0 && jj+1<=n) prod=prod*g.e[jj++];
		if(prod.len()==0){
			if(ii==1 && jj==n);
			else if(jj<n){
				F2.clear(); g2=g;
				for(int i=jj;i>=ii;i--){ // move from i to i+1
					g2.Elementary_transformation(i,-1);
					F2.push_back({1,n,i,-1});
				}
				search_tuples_of_short_elements(n,g2,m,idx_h,idx_g,F2,height+1,details);
			}else{
				F2.clear(); g2=g;
				for(int i=ii;i<=jj;i++) // move from i to i-ii+1
					for(int j=i-1;j>=i-ii+1;j--){ // move from j+1 to j
						g2.Elementary_transformation(j,1);
						F2.push_back({1,n,j,1});
					}
				search_tuples_of_short_elements(n,g2,m,idx_h,idx_g,F2,height+1,details);
			}

			F2.clear(); g2=g;
			for(int i=ii;i<jj;i++){ // move from i+1 to i
				g2.Elementary_transformation(i,1);
				F2.push_back({1,n,i,1});
			}
			//cout<<"cyclic permutation from "<<g<<" to "<<g2<<"\n";
			search_tuples_of_short_elements(n,g2,m,idx_h,idx_g,F2,height+1,details);
		}
	}
}

pair<bool,list<vector<int>>> search(Tuple<Ga3b2> g,Tuple<Ga3b2> h){
	int n=g.len(),m=h.len();
	vector<int> idx_h; idx_h.clear(); for(int i=0;i<m;i++) idx_h.push_back(is_short(h.e[i]));

	all_tuples.clear(); cont_search=true;
	search_tuples_of_short_elements(n,g,m,idx_h,vector<int>(),list<vector<int>>(),0,false);

	for(auto it:all_tuples){
		vector<int> idx_g=it.first;
		if(vector<int>(idx_g.begin()+n-m,idx_g.begin()+n)==idx_h){
			//cout<<"---find a way from "<<g<<" to (...) * "<<h<<"\n";
			list<vector<int>> F; F.clear();
			for(vector<int> idx_now_g=idx_g;idx_now_g.size()>0;idx_now_g=all_tuples[idx_now_g].first)
				F.insert(F.begin(),all_tuples[idx_now_g].second.begin(),all_tuples[idx_now_g].second.end());
			return make_pair(true,F);
		}
	}
	return make_pair(false,list<vector<int>>());
}

list<vector<int>> get_minimal_lexicographical_order(Tuple<Ga3b2> g){
	int n=g.len();
	all_tuples.clear(); cont_search=true;
	search_tuples_of_short_elements(n,g,1,{-1},vector<int>(),list<vector<int>>(),0,false);

	for(auto it:all_tuples){
		vector<int> idx_g=it.first;
		for(int i=0;i<18;i++){
			Tuple<Ga3b2> h=exceptional_tuples[i];
			int m=h.len();
			vector<int> idx_h; idx_h.clear(); for(int i=0;i<m;i++) idx_h.push_back(is_short(h.e[i]));

			if(idx_g==idx_h){
				list<vector<int>> F; F.clear();
				for(vector<int> idx_now_g=idx_g;idx_now_g.size()>0;idx_now_g=all_tuples[idx_now_g].first)
					F.insert(F.begin(),all_tuples[idx_now_g].second.begin(),all_tuples[idx_now_g].second.end());
				return F;
			}
		}
	}

	cout<<"---ERROR---\n";
	cout<<"g="<<g<<"\n";
	myassert(true,"it is impossible that the resulting tuple is not a known exceptional tuple");
}

/***search end***/

pair<Tuple<Ga3b2>,int> sort_concatenation(Tuple<Ga3b2> g_input,bool details){
	CR<Ga3b2> MM; MM.init(g_input);
	int n=g_input.len();

	int m_s=0,m_t=0,m_st=0,m_a=0,m_b=0,n_0=0,n_1=0,m=n;
	while(gotcha_a2_and_a2_and_a2(&MM,m)) ++n_1, m-=3;
	while(gotcha_a_and_a_and_a(&MM,m)) ++n_0, m-=3;
	while(gotcha_b_and_b(&MM,m)) ++m_b, m-=2;
	while(gotcha_a_and_a(&MM,m)) ++m_a, m-=2;
	while(gotcha_s_and_t(&MM,m)) ++m_st, m-=2;
	while(gotcha_t020202(&MM,m)) ++m_t, m-=6;
	while(gotcha_s020202(&MM,m)) ++m_s, m-=6;

	MM.Apply(shorten_induction(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m)),details,"h").second);
	cout<<"g="<<MM.h<<"\n";
	int length_exceptional_part=m;

	{
		Tuple<Ga3b2> h(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m));
		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(h,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
		if(m>0){
			myassert(cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0),"h.len()>0");
			myassert(cnt_a+cnt_a2+cnt_b>0 || abs(cnt_s-cnt_t)>0,"h.len()>0");
		}
	}

	/* apply diagonal conjugacies on subtuples with prod=1 */
	
	Tuple<Ga3b2> h(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m));
	Tuple<Ga3b2> g; g=h;
	for(int i=1;i<=m_s ;i++) g=bullet(g, Tuple<Ga3b2>({s0,s2,s0,s2,s0,s2}));
	for(int i=1;i<=m_t ;i++) g=bullet(g, Tuple<Ga3b2>({t0,t2,t0,t2,t0,t2}));
	for(int i=1;i<=m_st;i++) g=bullet(g, Tuple<Ga3b2>({s0,t0}));
	for(int i=1;i<=m_a ;i++) g=bullet(g, Tuple<Ga3b2>({a,a2}));
	for(int i=1;i<=m_b ;i++) g=bullet(g, Tuple<Ga3b2>({b,b}));
	for(int i=1;i<=n_0 ;i++) g=bullet(g, Tuple<Ga3b2>({a,a,a}));
	for(int i=1;i<=n_1 ;i++) g=bullet(g, Tuple<Ga3b2>({a2,a2,a2}));
	MM.init(g);

	cout<<"===after diagonal conjugacies on subtuples===\n";
	cout<<"g="<<MM.h<<"\n";
	
	/* -------- */
	
	while(gotcha_a_and_a_and_a(&MM,n) && gotcha_a2_and_a2_and_a2(&MM,n)){
		// (g_{n-5}, g_{n-4}, g_{n-3}, g_{n-2}, g_{n-1}, g_n) = (a, a, a, a^2, a^2, a^2)
		MM.Elementary_transformation(n-3,-1); // -> (a, a, a^2, a, a^2, a^2)
		MM.Elementary_transformation(n-2,-1); // -> (a, a, a^2, a^2, a, a^2)
		MM.Elementary_transformation(n-4,-1); // -> (a, a^2, a, a^2, a, a^2)
		--n_0; --n_1; m_a+=3;
	}
	cout<<"---avoid (a,a,a) and (a^2,a^2,a^2) occurring together---\n";
	cout<<"g="<<MM.h<<"\n";

	/* -------- */

	myassert(m_s==0 || m_t==0,"(s0,s2,s0,s2,s0,s2) and (t0,t2,t0,t2,t0,t2) cannot occur together");
	myassert(n_0==0 || n_1==0,"(a,a,a) and (a^2,a^2,a^2) cannot occur together");

	cout<<"===current concatenation form of the resulting tuple===\n";

	for(;;){
		if(m==0) break;

		h.init(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m));
		pair<bool,list<vector<int>>> ret;

		cout<<"---> ";
		cout<<h;
		if(m_s >0) cout<<" * (s0,s2,s0,s2,s0,s2)^"<<m_s;
		if(m_t >0) cout<<" * (t0,t2,t0,t2,t0,t2)^"<<m_t;
		if(m_st>0) cout<<" * (s0,t0)^"<<m_st;
		if(m_a >0) cout<<" * (a,a^2)^"<<m_a;
		if(m_b >0) cout<<" * (b,b)^"<<m_b;
		if(n_0 >0) cout<<" * (a,a,a)^"<<n_0;
		if(n_1 >0) cout<<" * (a^2,a^2,a^2)^"<<n_1;
		cout<<"\n";

		/***increase m_a and m_b***/

		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(h,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);

		// if h can end with (a,a^2)
		if(cnt_a>=1 && cnt_a2>=1){ MM.Apply(search(h, Tuple<Ga3b2>({a,a2})).second); ++m_a; m-=2; continue; }

		// if h can end with (b,b)
		if(cnt_b>=2){ MM.Apply(search(h, Tuple<Ga3b2>({b,b})).second); ++m_b; m-=2; continue; }

		/***increase m_st***/

		if(cnt_s>=1 && cnt_t>=1){
			// if h can end with (s0,t0)
			ret=search(h, Tuple<Ga3b2>({s0,t0}));
			if(ret.first){ MM.Apply(ret.second); ++m_st; m-=2; continue; }
			// if h can end with (s1,t1)
			ret=search(h, Tuple<Ga3b2>({s1,t1}));
			if(ret.first){ MM.Apply(ret.second); ++m_st; m-=2; continue; }
			// if h can end with (s2,t2)
			ret=search(h, Tuple<Ga3b2>({s2,t2}));
			if(ret.first){ MM.Apply(ret.second); ++m_st; m-=2; continue; }
		}

		if(cnt_s>=1 && m_t>=1){
			// if h can end with s0 and m_t>0
			ret=search(h, Tuple<Ga3b2>({s0}));
			if(ret.first){
				MM.Apply(ret.second); gotcha_t020202(&MM,MM.h.len());
				for(int i=n-5;i<=n;i++) // move from i to i-(n-6)+m
					for(int j=i-1;j>=i-(n-6)+m;j--) // move from j+1 to j
						MM.Elementary_transformation(j,1);
				// Now, M.h starts with (x,x,x,x,s0) * (t0,t2,t0,t2,t0,t2)
				ret=search(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m+6)), Tuple<Ga3b2>({s0,t0}));
				MM.Apply(ret.second); --m_t; ++m_st; m+=4; continue;
			}
			// if h can end with s1 and m_t>0
			ret=search(h, Tuple<Ga3b2>({s1}));
			if(ret.first){
				MM.Apply(ret.second); gotcha_t020202(&MM,MM.h.len());
				for(int i=n-5;i<=n;i++) // move from i to i-(n-6)+m
					for(int j=i-1;j>=i-(n-6)+m;j--) // move from j+1 to j
						MM.Elementary_transformation(j,1);
				// Now, M.h starts with (x,x,x,x,s1) * (t0,t2,t0,t2,t0,t2)
				ret=search(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m+6)), Tuple<Ga3b2>({s1,t1}));
				MM.Apply(ret.second); --m_t; ++m_st; m+=4; continue;
			}
			// if h can end with s2 and m_t>0
			ret=search(h, Tuple<Ga3b2>({s2}));
			if(ret.first){
				MM.Apply(ret.second); gotcha_t020202(&MM,MM.h.len());
				for(int i=n-5;i<=n;i++) // move from i to i-(n-6)+m
					for(int j=i-1;j>=i-(n-6)+m;j--) // move from j+1 to j
						MM.Elementary_transformation(j,1);
				// Now, M.h starts with (x,x,x,x,s2) * (t0,t2,t0,t2,t0,t2)
				ret=search(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m+6)), Tuple<Ga3b2>({s2,t2}));
				MM.Apply(ret.second); --m_t; ++m_st; m+=4; continue;
			}
		}

		if(cnt_t>=1 && m_s>=1){
			// if h can end with t0 and m_s>0
			ret=search(h, Tuple<Ga3b2>({t0}));
			if(ret.first){
				MM.Apply(ret.second); gotcha_s020202(&MM,MM.h.len());
				for(int i=n-5;i<=n;i++) // move from i to i-(n-6)+m
					for(int j=i-1;j>=i-(n-6)+m;j--) // move from j+1 to j
						MM.Elementary_transformation(j,1);
				// Now, M.h starts with (x,x,x,x,t0) * (s0,s2,s0,s2,s0,s2)
				ret=search(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m+6)), Tuple<Ga3b2>({s0,t0}));
				MM.Apply(ret.second); --m_s; ++m_st; m+=4; continue;
			}
			// if h can end with s1 and m_t>0
			ret=search(h, Tuple<Ga3b2>({t1}));
			if(ret.first){
				MM.Apply(ret.second); gotcha_s020202(&MM,MM.h.len());
				for(int i=n-5;i<=n;i++) // move from i to i-(n-6)+m
					for(int j=i-1;j>=i-(n-6)+m;j--) // move from j+1 to j
						MM.Elementary_transformation(j,1);
				// Now, M.h starts with (x,x,x,x,t1) * (s0,s2,s0,s2,s0,s2)
				ret=search(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m+6)), Tuple<Ga3b2>({s1,t1}));
				MM.Apply(ret.second); --m_s; ++m_st; m+=4; continue;
			}
			// if h can end with s2 and m_t>0
			ret=search(h, Tuple<Ga3b2>({t2}));
			if(ret.first){
				MM.Apply(ret.second); gotcha_s020202(&MM,MM.h.len());
				for(int i=n-5;i<=n;i++) // move from i to i-(n-6)+m
					for(int j=i-1;j>=i-(n-6)+m;j--) // move from j+1 to j
						MM.Elementary_transformation(j,1);
				// Now, M.h starts with (x,x,x,x,t2) * (s0,s2,s0,s2,s0,s2)
				ret=search(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m+6)), Tuple<Ga3b2>({s2,t2}));
				MM.Apply(ret.second); --m_s; ++m_st; m+=4; continue;
			}
		}

		/***increase m_a***/
		if(cnt_a>=1 && n_1>=1){
			MM.Apply(search(h, Tuple<Ga3b2>({a})).second); gotcha_a2_and_a2_and_a2(&MM,MM.h.len());
			for(int i=n-2;i<=n;i++) // move from i to i-(n-3)+m
				for(int j=i-1;j>=i-(n-3)+m;j--) // move from j+1 to j
					MM.Elementary_transformation(j,1);
			// Now, M.h starts with (x,x,x,x,a) * (a^2,a^2,a^2)
			ret=search(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m+3)), Tuple<Ga3b2>({a,a2}));
			MM.Apply(ret.second); --n_1; ++m_a; ++m; continue;
		}
		if(cnt_a2>=1 && n_0>=1){
			MM.Apply(search(h, Tuple<Ga3b2>({a2})).second); gotcha_a_and_a_and_a(&MM,MM.h.len());
			for(int i=n-2;i<=n;i++) // move from i to i-(n-3)+m
				for(int j=i-1;j>=i-(n-3)+m;j--) // move from j+1 to j
					MM.Elementary_transformation(j,1);
			// Now, M.h starts with (x,x,x,x,a^2) * (a,a,a)
			ret=search(Tuple<Ga3b2>(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+m+3)), Tuple<Ga3b2>({a,a2}));
			MM.Apply(ret.second); --n_0; ++m_a; ++m; continue;
		}
		
		/***increase m_s and m_t***/
		if(cnt_s>=6){
			ret=search(h, Tuple<Ga3b2>({s0,s2,s0,s2,s0,s2}));
			if(ret.first){ MM.Apply(ret.second); ++m_s; m-=6; continue; }
		}
		if(cnt_t>=6){
			ret=search(h, Tuple<Ga3b2>({t0,t2,t0,t2,t0,t2}));
			if(ret.first){ MM.Apply(ret.second); ++m_t; m-=6; continue; }
		}

		/***the exceptional tuple***/
		auto retF=get_minimal_lexicographical_order(h);
		MM.Apply(retF); break;
	}

	int Yn=MM.h.len();
	int Ym_s=0,Ym_t=0,Ym_st=0,Ym_a=0,Ym_b=0,Yn_0=0,Yn_1=0,Ym=Yn;
	
	while(gotcha_a2_and_a2_and_a2(&MM,Ym)) ++Yn_1, Ym-=3;
	while(gotcha_a_and_a_and_a(&MM,Ym)) ++Yn_0, Ym-=3;
	while(gotcha_b_and_b(&MM,Ym)) ++Ym_b, Ym-=2;
	while(gotcha_a_and_a(&MM,Ym)) ++Ym_a, Ym-=2;
	while(gotcha_s_and_t(&MM,Ym)) ++Ym_st, Ym-=2;
	while(gotcha_t020202(&MM,Ym)) ++Ym_t, Ym-=6;
	while(gotcha_s020202(&MM,Ym)) ++Ym_s, Ym-=6;
	//cout<< n<<" "<< m<<" "<< m_s<<" "<< m_t<<" "<< m_st<<" "<< m_a<<" "<< m_b<<" "<< n_0<<" "<< n_1<<"\n";
	//cout<<Yn<<" "<<Ym<<" "<<Ym_s<<" "<<Ym_t<<" "<<Ym_st<<" "<<Ym_a<<" "<<Ym_b<<" "<<Yn_0<<" "<<Yn_1<<"\n";

	myassert(n==Yn && m==Ym && m_s==Ym_s && m_t==Ym_t && m_st==Ym_st && m_a==Ym_a && m_b==Ym_b && n_0==Yn_0 & n_1==Yn_1,"same concatenation");

	cout<<"===final normal form===\n";
	h.init(vector<Ga3b2>(MM.h.e.begin(),MM.h.e.begin()+Ym));

	cout<<"---> ";
	cout<<h;
	if(Ym_s >0) cout<<" * (s0,s2,s0,s2,s0,s2)^"<<Ym_s;
	if(Ym_t >0) cout<<" * (t0,t2,t0,t2,t0,t2)^"<<Ym_t;
	if(Ym_st>0) cout<<" * (s0,t0)^"<<Ym_st;
	if(Ym_a >0) cout<<" * (a,a^2)^"<<Ym_a;
	if(Ym_b >0) cout<<" * (b,b)^"<<Ym_b;
	if(Yn_0 >0) cout<<" * (a,a,a)^"<<Yn_0;
	if(Yn_1 >0) cout<<" * (a^2,a^2,a^2)^"<<Yn_1;
	cout<<"\n";
	cout<<" = "<<MM.h<<"\n";


	/***check***/
	{
		myassert(each_component_is_short(h),"each component of h is short");
		int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t; get_cnts(h,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
		myassert(cnt_a<=2 && cnt_a2<=2 && cnt_b<=1,"cnt_a<=2 && cnt_a2<=2 && cnt_b<=1");
		myassert(cnt_a==0 || cnt_a2==0,"cnt_a==0 || cnt_a2==0");
		myassert(cnt_s==0 || cnt_t==0,"cnt_s==0 || cnt_t==0");
	}

	return make_pair(MM.h,length_exceptional_part);
}
