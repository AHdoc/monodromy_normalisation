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

bool each_component_is_almost_short(Tuple<Ga3b2> g){
	for(Ga3b2 x:g.e)
		if(is_almost_short(x)>=0);
		else return false;
	return true;
}

int S2_complexity(Ga3b2 g){
	if(g.len()==0) return 0;
	int q=0;
	for(;;){
		if(is_almost_short(g)){
			if(g==Ga3b2({"a","b","a","b","a"}) || g==Ga3b2({"a^2","b","a^2","b","a^2"})){
				if(q>0) return 2*q;
				else return 1;
			}else return 2*q;
		}
		pop_front_back(g.e); ++q;
	}
}

int S2_complexity(Tuple<Ga3b2> g){
	int ret=0;
	for(int i=0;i<g.len();i++) ret+=S2_complexity(g.e[i]);
	return ret;
}

/*--Operation 1--------------*/
bool OperationA(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=1;i+1<=n;i++){
		Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i];
		for(;;){
			int row=is_almost_short(g1), col=is_almost_short(g2);
			if(row!=-1 && col!=-1 && almost_short_table[row][col]){
				M->Contraction(i,i+1);
				return true;
			}
			int len1=g1.len(),len2=g2.len();
			if(len1>=2 && len2>=2){
				if(g1.e[0]=="b" && g1.e[len1-1]=="b" && g2.e[0]=="b" && g2.e[len2-1]=="b"){
					pop_front_back(g1.e); pop_front_back(g2.e);
				}else if(g1.e[0]=="a^2" && g1.e[len1-1]=="a" && g2.e[0]=="a^2" && g2.e[len2-1]=="a"){
					pop_front_back(g1.e); pop_front_back(g2.e);
				}else if(g1.e[0]=="a" && g1.e[len1-1]=="a^2" && g2.e[0]=="a" && g2.e[len2-1]=="a^2"){
					pop_front_back(g1.e); pop_front_back(g2.e);
				}else break;
			}
		}
	}
	return false;
}

/*--Operation 1'--------------*/
bool OperationB(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=1;i+2<=n;i++){
		Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i],g3=M->h.e[i+1];
		for(;;){
			if((g1==Ga3b2({"a^2","b","a","b"}) && g2==Ga3b2({"b","a","b"}) && g3==Ga3b2({"b","a","b","a^2"})) ||
			   (g1==Ga3b2({"a","b","a^2","b"}) && g2==Ga3b2({"b","a^2","b"}) && g3==Ga3b2({"b","a^2","b","a"}))){
				M->Contraction(i,i+2);
			}
			int len1=g1.len(),len2=g2.len(),len3=g3.len();
			if(len1>=2 && len2>=2 && len3>=2){
				if(g1.e[0]=="b" && g1.e[len1-1]=="b" && g2.e[0]=="b" && g2.e[len2-1]=="b" && g3.e[0]=="b" && g3.e[len3-1]=="b"){
					pop_front_back(g1.e); pop_front_back(g2.e); pop_front_back(g3.e);
				}else if(g1.e[0]=="a^2" && g1.e[len1-1]=="a" && g2.e[0]=="a^2" && g2.e[len2-1]=="a" && g3.e[0]=="a^2" && g3.e[len3-1]=="a"){
					pop_front_back(g1.e); pop_front_back(g2.e); pop_front_back(g3.e);
				}else if(g1.e[0]=="a" && g1.e[len1-1]=="a^2" && g2.e[0]=="a" && g2.e[len2-1]=="a^2" && g3.e[0]=="a" && g3.e[len3-1]=="a^2"){
					pop_front_back(g1.e); pop_front_back(g2.e); pop_front_back(g3.e);
				}else break;
			}
		}
	}
	return false;
}

/*--Operation 2--------------*/
bool OperationC(CR<Ga3b2>* M){
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

/*--Elementary transformation strictly-decrease S2_complexity--------------*/
bool OperationD(CR<Ga3b2>* M,int il,int ir){ // restrict that 1<=il<=i<=ir<=n-1
	for(int i=il;i<=ir;i++) for(int epsilon=-1;epsilon<=1;epsilon+=2){
		Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i,epsilon);
		if(S2_complexity(h1)<S2_complexity(M->h)){
			M->Elementary_transformation(i,epsilon);
			return true;
		}
	}
	for(int i1=il;i1<=ir;i1++) for(int epsilon1=-1;epsilon1<=1;epsilon1+=2)
		for(int i2=il;i2<=ir;i2++) for(int epsilon2=-1;epsilon2<=1;epsilon2+=2){
			Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i1,epsilon1);
			Tuple<Ga3b2> h2=h1  ; h2.Elementary_transformation(i2,epsilon2);
			if(S2_complexity(M->h)==S2_complexity(h1) && S2_complexity(h1)>S2_complexity(h2)){
				M->Elementary_transformation(i1,epsilon1);
				M->Elementary_transformation(i2,epsilon2);
				return true;
			}
		}
	return false;
}

/*--Handle the case h_i=Q^{-1}aQ or ... and l(h_{i-1})>l(h_i)--------------*/
bool OperationE(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=2;i<=n;i++){
		Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i-1,1);
		if(S2_complexity(h1)==S2_complexity(M->h) && h1.e[i-2].len()<M->h.e[i-2].len()){
			M->Elementary_transformation(i-1,1);
			return true;
		}
	}
	return false;
}

/******************/

/*--an induction that--------------*/
/*----transform (g1,...,gn) in Ga3b2 whose components are conjuagates of almost short elements--*/
/*--into--*/
/*----(h1,...,hm) * g_non_inverse_free--*/
/*--s.t.--*/
/*----each component of (h1,...,hm) is almost short--*/
/*--return mp(h,F)--*/
pair<Tuple<Ga3b2>,list<vector<int>>> almost_shorten_induction(Tuple<Ga3b2> g_input,bool details,string namestr){
	cout<<"| +====start the almost_shorten_induction on "+namestr+"===\n";
	cout<<"| | "<<namestr<<"="<<g_input<<"\n";

	CR<Ga3b2> M;
	Tuple<Ga3b2> g_non_inverse_free;
	list<vector<int>> F;

	M.init(g_input);
	/***********/
	for(;;){
		if(details){
			cout<<"| +-------\n";
			cout<<"| | M.h="<<M.h<<"\n";
			cout<<"| | M.H="<<M.H<<"\n";
		}
		if(OperationA(&M))               continue;
		if(OperationB(&M))               continue;
		if(OperationC(&M))               continue;
		if(OperationD(&M,1,M.h.len()-1)) continue;
		if(OperationE(&M))               continue;
		break;
	}
	/***********/
	cout<<"| +----after the induction---\n";
	cout<<"| | M.h="<<M.h<<"\n";
	cout<<"| | M.H="<<M.H<<"\n";

	//careful_restorations(&M,false);

	/*
	g_non_inverse_free.clear();
	F.clear();
	while(find_pair(&M,&g_non_inverse_free,&F))   continue;
	while(find_triple(&M,&g_non_inverse_free,&F)) continue;
	*/

	// Fron now on, each component of M.h is short;
	//myassert(each_component_is_short(M.h),"fron now on, each component of M.h is short");

	//F.insert(F.end(),M.F.begin(),M.F.end()); M.F.clear();

	//cout<<"| | ===conclusion after shorten_induction===\n";
	//cout<<"| | "<<namestr<<"_h="<<M.h<<"\n";
	//cout<<"| | "<<namestr<<"_non_inverse_free="<<g_non_inverse_free<<"\n";

	/***********/

	// if(details) cout<<"| | ---check elementary transformations in F---\n";
	// if(details) cout<<"| | F.size()="<<F.size()<<"\n";
	// Tuple<Ga3b2> g=g_input;
	// for(auto it:F){
	// 	int i=it[2],epsilon=it[3];
	// 	g.Elementary_transformation(i,epsilon);
	// }
	// myassert(g==bullet(M.h,g_non_inverse_free),"g==h dot g_non_inverse_free");
	//cout<<"| +==="+namestr+" is transformed into  \""<<namestr<<"_h * "<<namestr<<"_non_inverse_free\"  by "<<F.size()<<" elementary transformations\n";

	//myassert(each_component_is_short(M.h),"after the shorten_induction, each component must be short");
	return make_pair(M.h,F);
}

/***************************/
