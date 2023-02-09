Ga3b2 short_elems[9]={
	Ga3b2({"a"}), Ga3b2({"a^2"}), Ga3b2({"b"}),
	Ga3b2({"a^2","b"}), Ga3b2({"a","b","a"}), Ga3b2({"b","a^2"}),
	Ga3b2({"b","a"}), Ga3b2({"a^2","b","a^2"}), Ga3b2({"a","b"})
};

Ga3b2 a({"a"}),a2({"a^2"}),b({"b"});
Ga3b2 s0({"a^2","b"}),s1({"a","b","a"}),s2({"b","a^2"}),t0({"b","a"}),t1({"a^2","b","a^2"}),t2({"a","b"});

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
		int x=i-(epsilon+1)/2,y=i+(epsilon-1)/2;
		Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i,epsilon);
		if(S_complexity(h1)<S_complexity(M->h) && !(is_short(M->h.e[x])>=0 && is_short(h1.e[y])==-1)){
			M->Elementary_transformation(i,epsilon);
			return true;
		}
	}
	/*To handle cases including (s1, a2ba) and (t1, aba2),
	  we have to consider a sequence of elementary transformtaions of length 2.*/
	for(int i1=il;i1<=ir;i1++) for(int epsilon1=-1;epsilon1<=1;epsilon1+=2)
		for(int i2=il;i2<=ir;i2++) for(int epsilon2=-1;epsilon2<=1;epsilon2+=2){
			int x1=i1-(epsilon1+1)/2,y1=i1+(epsilon1-1)/2;
			int x2=i2-(epsilon2+1)/2,y2=i2+(epsilon2-1)/2;
			Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i1,epsilon1);
			Tuple<Ga3b2> h2=h1 ; h2.Elementary_transformation(i2,epsilon2);
			if(S_complexity(h2)<S_complexity(M->h) &&
				!(is_short(M->h.e[x1])>=0 && is_short(h1.e[y1])==-1) &&
				!(is_short(h1.e[x2])>=0 && is_short(h2.e[y2])==-1)
			){
				M->Elementary_transformation(i1,epsilon1);
				M->Elementary_transformation(i2,epsilon2);
				return true;
			}
		}
	return false;
}

/*--Handle the case h_i=Q^{-1}aQ and l(h_{-1})>l(h_i)--------------*/
bool Operation4(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=2;i<=n;i++){
		Ga3b2 gi=M->h.e[i-1];
		if(gi.len()==0) continue;
		else if(is_short(gi)==-1 && (gi.e[gi.len()/2]=="a" || gi.e[gi.len()/2]=="a^2")); // gi=Q^{-1}aQ or Q^{-1}a^2Q
		else if(is_short(gi)==0 || is_short(gi)==1); // gi=a or gi=a^2
		else continue;

		int x=i-1-1,y=i-1;
		Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i-1,1);
		if(S_complexity(h1)<=S_complexity(M->h) &&
			!(is_short(M->h.e[x])>=0 && is_short(h1.e[y])==-1) &&
			h1.e[i-2].len()<M->h.e[i-2].len()
		){
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

/*--an induction that--------------*/
/*----transform (g1,...,gn) in Ga3b2 whose components are conjuagates of short elements--*/
/*----into-(h1,...,hm) bullet g_non_inverse_free--------------------------------------------*/
/*--s.t. each component of (h1,...,hm) is short--*/
/*--return mp(h,F)--*/
pair<Tuple<Ga3b2>,list<vector<int>>> shorten_induction(Tuple<Ga3b2> g_input,bool details,string namestr){
	cout<<"+====start the shorten_induction on "+namestr+" (use h to denote the tuple in this process)===\n";
	cout<<"| h="<<g_input<<"\n";

	CR<Ga3b2> M;
	Tuple<Ga3b2> g_non_inverse_free;
	list<vector<int>> F;

	M.init(g_input);
	g_non_inverse_free.clear();
	F.clear();

	for(int t=1;;t++){
		if(details) cout<<"| ===the "<<t<<"-th induction on "<<namestr<<" starts===\n";
		if(details) cout<<"| h="<<M.h<<"\n";
		if(details) cout<<"| H="<<M.H<<"\n";
		if(details) cout<<"| g_non_inverse_free="<<g_non_inverse_free<<"\n";

		myassert(M.F.empty(),"F is empty");

		/***********/

		for(;;){
			if(Operation1(&M))              continue;
			if(Operation2(&M))              continue;
			if(Operation3(&M,1,M.h.len()-1)) continue;
			if(Operation4(&M))              continue;
			break;
		}
		if(details) cout<<"| ---the "<<t<<"-th induction is done---\n";
		if(details) cout<<"| h="<<M.h<<"\n";
		if(details) cout<<"| H="<<M.H<<"\n";

		/***********/
		
		for(;;){
			int k=M.Restoration();
			if(k==-1) break;
			while(Operation3(&M,k,k)) continue;
		}
		if(details) cout<<"| ---restoration is done---\n";
		if(details) cout<<"| h="<<M.h<<"\n";
		if(details) cout<<"| H="<<M.H<<"\n";

		/***********/

		if(find_triple(&M,&g_non_inverse_free,&F)) continue;
		if(find_pair(&M,&g_non_inverse_free,&F))   continue;

		/***********/
		// Fron now on, each component of M.h is short;
		myassert(each_component_is_short(M.h),"fron now on, each component of M.h is short");
		// there are at most 2 a's, at most 2 a's, at most 1 b;
		// a and a^2 cannot occur together; si and ti cannot occur together.

		break;
	}
	F.insert(F.end(),M.F.begin(),M.F.end()); M.F.clear();

	cout<<"| ===conclusion after shorten_induction===\n";
	cout<<"| h="<<M.h<<"\n";
	cout<<"| g_non_inverse_free="<<g_non_inverse_free<<"\n";

	/***********/

	if(details) cout<<"| ---check elementary transformations in F---\n";
	if(details) cout<<"| F.size()="<<F.size()<<"\n";
	Tuple<Ga3b2> g=g_input;
	for(auto it:F){
		int i=it[2],epsilon=it[3];
		g.Elementary_transformation(i,epsilon);
	}
	myassert(g==bullet(M.h,g_non_inverse_free),"g==h dot g_non_inverse_free");
	cout<<"+==="+namestr+" is transformed into  \"h bullet g_non_inverse_free\"  by "<<F.size()<<" elementary transformations\n";
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

/*--Suppose that h1,...,hn are short s.t. at most 1 of them is equal to one of a,a^2,b--*/
/*--Rearrange them and transform h into a tuple of the form (a/a^2/b,s---s,t---t)--*/
void Rearrangement(CR<Ga3b2>* M){
	int n=M->h.len(),c_unexcep=0;
	for(int i=1;i<=n;i++)
		if(M->h.e[i-1]==a || M->h.e[i-1]==a2 || M->h.e[i-1]==b)
			++c_unexcep;
	myassert(c_unexcep==0 || c_unexcep==1,"at most 1 of h1,...,hn is equal to one of a,a^2,b");
	if(c_unexcep==1)
		for(;!(M->h.e[0]==a || M->h.e[0]==a2 || M->h.e[0]==b);) cyclic_permutation(M); // h1=a/a^2/b
	for(;;){
		int ok=1;
		for(int i=1;i<=n-1;i++){
			// avoid the appearance of (ti,sj)
			if((M->h.e[i-1]==t0 || M->h.e[i-1]==t1 || M->h.e[i-1]==t2)&&(M->h.e[i]==s0 || M->h.e[i]==s1 || M->h.e[i]==s2)){
				if(is_short(M->h.e[i])+3==is_short(M->h.e[i-1])){ // (t0,s0), (t1,s1) or (s2,s2)
					for(int j=i-1;j>=1;j--) M->Elementary_transformation(j,1);
					for(int j=i  ;j>=2;j--) M->Elementary_transformation(j,1);
					ok=3; break;
				}
				if((is_short(M->h.e[i])-3+2)%3==(is_short(M->h.e[i-1])-6)%3){ // (t_{k+1},s_{k-1}) -> (s_k, t_{k+1})
					M->Elementary_transformation(i,-1);
					ok=2; break;
				}
				if((is_short(M->h.e[i])-3+1)%3==(is_short(M->h.e[i-1])-6)%3){ // (t_{k+1},s_{k}) -> (s_k, t_{k-1})
					M->Elementary_transformation(i,1);
					ok=2; break;
				}
			}
		}
		if(ok==1 || ok==3) break;
	}
	cout<<"--Rearrange (h1,...,hn) of short elements s.t. at most 1 of them is equal to a/a^2/b--\n";
	cout<<"h="<<M->h<<"\n";
}

/*--make inverse-free--------------*/
/*--a tuple in Ga3b2 whose components are conjuagates of short elements--------------*/
void inverse_free(Tuple<Ga3b2> g_input,bool details){
	CR<Ga3b2> MM; MM.init(g_input);
	Tuple<Ga3b2> g; g=g_input; // g must be a prefix of M.h
	for(;;){
		if(each_component_is_short(g)){
			cout<<"===get a tuple of short elements==\n";
			cout<<"h="<<g<<"\n";
			
			int cnt_a,cnt_a2,cnt_b,cnt_s,cnt_t;
			get_cnts(g,&cnt_a,&cnt_a2,&cnt_b,&cnt_s,&cnt_t);
			
			cout<<"cnt_a, cnt_a2, cnt_b, cnt_s, cnt_t = "<<cnt_a<<", "<<cnt_a2<<", "<<cnt_b<<", "<<cnt_s<<", "<<cnt_t<<"\n";
			myassert(cnt_a<=2 && cnt_a2<=2 && cnt_b<=1 && (cnt_a==0 || cnt_a2==0),"no pair of inverse elements nor (l,l,l) with l^3=1");

			CR<Ga3b2> M; M.init(g);

			if(cnt_a==2 || cnt_a2==2){ // I_a=2
				Ga3b2 x=(cnt_a==2?a:a2);
				auto poss=get_positions(M.h,x);
				for(int j=poss[0]-1;j>=1;j--) M.Elementary_transformation(j,1);
				for(int j=poss[1]-1;j>=2;j--) M.Elementary_transformation(j,1);
				M.Contraction(1,2);
				
				cout<<"---combine two powers of a---\n";
				cout<<"h="<<M.h<<"\n";
				if(details) cout<<"H="<<M.H<<"\n";

				for(auto it:M.F){
					if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
					if(it[0]==2) MM.Contraction(it[1],it[2]);
				}

				auto ret=shorten_induction(M.h,details,"h");
				for(auto it:ret.second){
					if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
					if(it[0]==2) MM.Contraction(it[1],it[2]);
				}
				g=ret.first;

				continue;
			}

			if(cnt_a+cnt_a2==1 && cnt_b==1){ // I_a=1=I_b
				Ga3b2 x=(cnt_a==1?a:a2);
				for(;!(M.h.e[0]==x);) cyclic_permutation(&M); // h1=x
				cout<<"---cyclic_permutation--\n";
				cout<<"h="<<M.h<<"\n";

				bool ok=false;
				if(!ok){
					for(int i=2;!ok && i<=M.h.len();i++)
						if(is_short(M.h.e[i-1])>=3 && short_table[is_short(M.h.e[0])][is_short(M.h.e[i-1])]){
							
							for(int j=i-1;j>=2;j--) M.Elementary_transformation(j,1); // (a,ba,...)
							M.Contraction(1,2);

							cout<<"---combine a^e and one of si, tj---\n";
							cout<<"h="<<M.h<<"\n";
							if(details) cout<<"H="<<M.H<<"\n";

							for(auto it:M.F){
								if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
								if(it[0]==2) MM.Contraction(it[1],it[2]);
							}
							
							auto ret=shorten_induction(M.h,details,"h");
							for(auto it:ret.second){
								if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
								if(it[0]==2) MM.Contraction(it[1],it[2]);
							}
							g=ret.first;

							ok=true;
						}
				}

				if(!ok){
					cyclic_permutation(&M); // hn=x
					cout<<"---cyclic_permutation--\n";
					cout<<"h="<<M.h<<"\n";
					for(int i=1;!ok && i<=M.h.len()-1;i++)
						if(is_short(M.h.e[i-1])>=3 && short_table[is_short(M.h.e[i-1])][is_short(M.h.e[M.h.len()-1])]){

							for(int j=i;j<=M.h.len()-1;j++) M.Elementary_transformation(j,-1); // (...,ba,a)
							M.Contraction(M.h.len()-1,M.h.len());

							cout<<"---combine a^e and one of si, tj---\n";
							cout<<"h="<<M.h<<"\n";
							if(details) cout<<"H="<<M.H<<"\n";

							for(auto it:M.F){
								if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
								if(it[0]==2) MM.Contraction(it[1],it[2]);
							}

							auto ret=shorten_induction(M.h,details,"h");
							for(auto it:ret.second){
								if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
								if(it[0]==2) MM.Contraction(it[1],it[2]);
							}
							g=ret.first;

							ok=true;
						}
				}

				// Otherwise, A={a^e,b,a^eba^e}.
				if(!ok){
					myassert(M.h.e[M.h.len()-1]==x,"hn=x");
					for(int i=M.h.len()-1;i>=1;i--)
						if((x==a && M.h.e[i-1]==s1)||(x==a2 && M.h.e[i-1]==t1)){
							for(int j=i;j<=M.h.len()-2;j++) M.Elementary_transformation(j,-1);
							break;
						}
					for(auto poss=get_positions(M.h,x);poss[0]!=2;poss=get_positions(M.h,x)) cyclic_permutation(&M); // h2=x
					for(int i=3;i<=M.h.len();i++)
						if((x==a && M.h.e[i-1]==s1)||(x==a2 && M.h.e[i-1]==t1)){
							for(int j=i-1;j>=3;j--) M.Elementary_transformation(j,1);
							break;
						}
					if(x==a)  myassert(M.h.e[0]==s1 && M.h.e[1]==a  && M.h.e[2]==s1,"(aba,a,aba)");
					if(x==a2) myassert(M.h.e[0]==t1 && M.h.e[1]==a2 && M.h.e[2]==s2,"(aabaa,aa,aabaa)");

					M.Elementary_transformation(2,1);
					M.Contraction(1,2);
					M.Contraction(1,2);
					
					cout<<"---transform and contract (aba,a,aba) into a^2---\n";
					cout<<"h="<<M.h<<"\n";
					if(details) cout<<"H="<<M.H<<"\n";

					for(auto it:M.F){
						if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
						if(it[0]==2) MM.Contraction(it[1],it[2]);
					}

					auto ret=shorten_induction(M.h,details,"h");
					for(auto it:ret.second){
						if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
						if(it[0]==2) MM.Contraction(it[1],it[2]);
					}
					g=ret.first;

					ok=true;
				}
				continue;
			}

			if(cnt_a+cnt_a2+cnt_b<=1 && cnt_s>0 && cnt_t>0){
				Rearrangement(&M);
				for(auto it:M.F){
					if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
					if(it[0]==2) MM.Contraction(it[1],it[2]);
				}
				
				auto ret=shorten_induction(M.h,details,"h");
				for(auto it:ret.second){
					if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
					if(it[0]==2) MM.Contraction(it[1],it[2]);
				}
				g=ret.first;

				continue;
			}

			/*--Fron now on, g is inverse-free!--*/

			break;
		}else{
			auto ret=shorten_induction(MM.h,details,"g");
			for(auto it:ret.second){
				if(it[0]==1) MM.Elementary_transformation(it[2],it[3]);
				if(it[0]==2) MM.Contraction(it[1],it[2]);
			}
			g=ret.first;
		}
	}

	/*******************/

	for(;;){
		int k=MM.Restoration();
		if(k==-1) break;
		while(Operation3(&MM,k,k)) continue;
	}
	cout<<"===the final restoration is done===\n";
	cout<<"h="<<MM.h<<"\n";

	return;
}
