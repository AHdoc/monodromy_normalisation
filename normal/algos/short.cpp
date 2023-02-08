Ga3b2 short_elems[9]={
	Ga3b2({"a"}), Ga3b2({"a^2"}), Ga3b2({"b"}),
	Ga3b2({"a^2","b"}), Ga3b2({"a","b","a"}), Ga3b2({"b","a^2"}),
	Ga3b2({"b","a"}), Ga3b2({"a^2","b","a^2"}), Ga3b2({"a","b"})
};

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

CR<Ga3b2> M;
Tuple<Ga3b2> g_non_inverse_free;
list<vector<int>> F;

/*--Operation 1--------------*/
bool Operation1(){
	int n=M.h.len();
	for(int i=1;i+1<=n;i++){
		Ga3b2 g1=M.h.e[i-1],g2=M.h.e[i];
		Ga3b2 tau_1=get_tau_in_expression(g1),Q_1=get_Q_in_expression(g1);
		Ga3b2 tau_2=get_tau_in_expression(g2),Q_2=get_Q_in_expression(g2);
		int row=is_short(tau_1), col=is_short(tau_2); // row==-1 iff tau_1==identity
		if(Q_1==Q_2 && row!=-1 && col!=-1 && short_table[row][col]){
			if(is_short(g1) || (tau_1*tau_2).len()==0 || (tau_1.is_power_of_a() && tau_2.is_power_of_a())){
				M.Contraction(i,i+1);
				return true;
			}
		}
	}
	return false;
}

/*--Operation 2--------------*/
bool Operation2(){
	int n=M.h.len();
	for(int i=1;i<n;i++){
		if(M.h.e[i-1].len()==0 && M.h.e[i].len()>=1){
			for(int j=i;j<n;j++){
				M.Elementary_transformation(j,-1);
				return true;
			}
		}
	}
	return false;
}

/*--Elementary transformation avoiding that a short becomes long--------------*/
bool Operation3(int il,int ir){ // restrict that 1<=il<=i<=ir<=n-1
	for(int i=il;i<=ir;i++) for(int epsilon=-1;epsilon<=1;epsilon+=2){
		int x=i-(epsilon+1)/2,y=i+(epsilon-1)/2;
		Tuple<Ga3b2> h1=M.h; h1.Elementary_transformation(i,epsilon);
		if(S_complexity(h1)<S_complexity(M.h) && !(is_short(M.h.e[x])>=0 && is_short(h1.e[y])==-1)){
			M.Elementary_transformation(i,epsilon);
			return true;
		}
	}
	/*To handle cases including (s1, a2ba) and (t1, aba2),
	  we have to consider a sequence of elementary transformtaions of length 2.*/
	for(int i1=il;i1<=ir;i1++) for(int epsilon1=-1;epsilon1<=1;epsilon1+=2)
		for(int i2=il;i2<=ir;i2++) for(int epsilon2=-1;epsilon2<=1;epsilon2+=2){
			int x1=i1-(epsilon1+1)/2,y1=i1+(epsilon1-1)/2;
			int x2=i2-(epsilon2+1)/2,y2=i2+(epsilon2-1)/2;
			Tuple<Ga3b2> h1=M.h; h1.Elementary_transformation(i1,epsilon1);
			Tuple<Ga3b2> h2=h1 ; h2.Elementary_transformation(i2,epsilon2);
			if(S_complexity(h2)<S_complexity(M.h) &&
				!(is_short(M.h.e[x1])>=0 && is_short(h1.e[y1])==-1) &&
				!(is_short(h1.e[x2])>=0 && is_short(h2.e[y2])==-1)
			){
				M.Elementary_transformation(i1,epsilon1);
				M.Elementary_transformation(i2,epsilon2);
				return true;
			}
		}
	return false;
}

/*--Handle the case h_i=Q^{-1}aQ and l(h_{-1})>l(h_i)--------------*/
bool Operation4(){
	int n=M.h.len();
	for(int i=2;i<=n;i++){
		Ga3b2 gi=M.h.e[i-1];
		if(gi.len()==0) continue;
		else if(is_short(gi)==-1 && (gi.e[gi.len()/2]=="a" || gi.e[gi.len()/2]=="a^2")); // gi=Q^{-1}aQ or Q^{-1}a^2Q
		else if(is_short(gi)==0 || is_short(gi)==1); // gi=a or gi=a^2
		else continue;

		int x=i-1-1,y=i-1;
		Tuple<Ga3b2> h1=M.h; h1.Elementary_transformation(i-1,1);
		if(S_complexity(h1)<=S_complexity(M.h) &&
			!(is_short(M.h.e[x])>=0 && is_short(h1.e[y])==-1) &&
			h1.e[i-2].len()<M.h.e[i-2].len()
		){
			M.Elementary_transformation(i-1,1);
			return true;
		}
	}
	return false;
}

/*--find a triple (l,l,l) s.t. l^3=1--------------*/
bool find_triple(){
	int n=M.h.len();
	for(int i1=1;i1<=n;i1++) for(int i2=i1+1;i2<=n;i2++) for(int i3=i2+1;i3<=n;i3++){
		Ga3b2 g1=M.h.e[i1-1];
		Ga3b2 g2=M.h.e[i2-1];
		Ga3b2 g3=M.h.e[i3-1];
		if(g1==g2 && g1==g3 && pow(g1,3).len()==0){
			for(int j=i3;j<=n-1;j++) M.Elementary_transformation(j,-1);
			for(int j=i2;j<=n-2;j++) M.Elementary_transformation(j,-1);
			for(int j=i1;j<=n-3;j++) M.Elementary_transformation(j,-1);

			g_non_inverse_free.e.insert(g_non_inverse_free.e.begin(),g3);
			g_non_inverse_free.e.insert(g_non_inverse_free.e.begin(),g2);
			g_non_inverse_free.e.insert(g_non_inverse_free.e.begin(),g1);
			M.h.e.pop_back(); M.h.e.pop_back(); M.h.e.pop_back();
			M.H.e.pop_back(); M.H.e.pop_back(); M.H.e.pop_back();

			F.insert(F.end(),M.F.begin(),M.F.end());
			M.F.clear();
			return true;
		}
	}
	return false;
}

/*--find a pair (g1,g2) of inverse elements--------------*/
bool find_pair(){
	int n=M.h.len();
	for(int i1=1;i1<=n;i1++) for(int i2=i1+1;i2<=n;i2++){
		Ga3b2 g1=M.h.e[i1-1];
		Ga3b2 g2=M.h.e[i2-1];
		if((g1*g2).len()==0){
			for(int j=i2;j<=n-1;j++) M.Elementary_transformation(j,-1);
			for(int j=i1;j<=n-2;j++) M.Elementary_transformation(j,-1);

			g_non_inverse_free.e.insert(g_non_inverse_free.e.begin(),g2);
			g_non_inverse_free.e.insert(g_non_inverse_free.e.begin(),g1);
			M.h.e.pop_back(); M.h.e.pop_back();
			M.H.e.pop_back(); M.H.e.pop_back();

			F.insert(F.end(),M.F.begin(),M.F.end());
			M.F.clear();
			return true;
		}
	}
	return false;
}



/*--an induction to shorten--------------*/
/*--a tuple in Ga3b2 whose components are conjuagates of short elements--------------*/
void shorten_induction(Tuple<Ga3b2> g_input){
	M.init(g_input);
	g_non_inverse_free.clear();
	F.clear();

	for(int t=1;;t++){
		cout<<"===the "<<t<<"-th induction starts===\n";
		cout<<"h="<<M.h<<"\n";
		cout<<"H="<<M.H<<"\n";
		cout<<"g_non_inverse_free="<<g_non_inverse_free<<"\n";

		myassert(M.F.empty(),"F is empty");

		/***********/

		for(;;){
			//cout<<"h="<<M.h<<"\n";
			if(Operation1())              continue;
			if(Operation2())              continue;
			if(Operation3(1,M.h.len()-1)) continue;
			if(Operation4())              continue;
			break;
		}
		cout<<"---the "<<t<<"-th induction is done---\n";
		cout<<"h="<<M.h<<"\n";
		cout<<"H="<<M.H<<"\n";

		/***********/
		
		for(;;){
			int k=M.Restoration();
			if(k==-1) break;
			while(Operation3(k,k)) continue;
		}
		cout<<"---restoration is done---\n";
		cout<<"h="<<M.h<<"\n";
		cout<<"H="<<M.H<<"\n";

		/***********/

		if(find_triple()) continue;
		if(find_pair())   continue;

		break;
	}
	cout<<"===conclusion===\n";
	cout<<"h="<<M.h<<"\n";
	cout<<"g_non_inverse_free="<<g_non_inverse_free<<"\n";

	/***********/

	cout<<"---check elementary transformations in F---\n";
	cout<<"F.size()="<<F.size()<<"\n";
	Tuple<Ga3b2> g=g_input;
	for(auto it:F){
		int i=it[2],epsilon=it[3];
		g.Elementary_transformation(i,epsilon);
	}
	myassert(g==bullet(M.h,g_non_inverse_free),"g==h dot g_non_inverse_free");
	cout<<"resulting g=h bullet g_non_inverse_free\n";

	// return something
}
