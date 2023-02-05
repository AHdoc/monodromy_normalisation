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

void solve_short(Tuple<Ga3b2> g){ // g is a tuple in Ga3b2 whose components are conjuagates of short elements
	CR<Ga3b2> M; M.init(g);

	for(bool inducible;;){
		inducible=false;
		int n=M.h.len();
		/*----------------*/
		for(int i=1;i+1<=n;i++){
			int row=is_short(M.h.e[i-1]), col=is_short(M.h.e[i]);
			if(row!=-1 && col!=-1 && short_table[row][col]){
				M.Contraction(i,i+1);
				inducible=true;
				break;
			}
		}
		if(inducible) continue;
		/*----------------*/
		for(int i=1;i<n;i++){
			if(M.h.e[i-1].len()==0 && M.h.e[i].len()>=1){
				for(int j=i;j<n;j++)
					M.Elementary_transformation(j,-1);
				inducible=true;
				break;
			}
		}
		if(inducible) continue;
		/*----------------*/
		for(int i=1;i<n;i++){
			Tuple<Ga3b2> h2=M.h; h2.Elementary_transformation(i,-1);
			if(S_complexity(h2)<S_complexity(M.h)){
				M.Elementary_transformation(i,-1);
				inducible=true;
				break;
			}
			h2=M.h; h2.Elementary_transformation(i,1);
			if(S_complexity(h2)<S_complexity(M.h)){
				M.Elementary_transformation(i,1);
				inducible=true;
				break;
			}
		}
		if(inducible) continue;
		/*----------------*/


		if(!inducible) break;
	}
	cout<<"h="<<M.h<<"\n";
	cout<<"H="<<M.H<<"\n";
}