/*--Operation 1--------------*/
bool OperationA(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=1;i+1<=n;i++){
		Ga3b2 g1=M->h.e[i-1],g2=M->h.e[i];

		if(g1.len()==0 || g2.len()==0) continue;
		
		int row=is_almost_short(g1), col=is_almost_short(g2);
		if(row!=-1 && col!=-1 && almost_short_table[row][col]){
			M->Contraction(i,i+1);
			return true;
		}
		if(row==-1 && col==-1){
			Ga3b2 Q1=get_Q_in_expression_2(g1),Q2=get_Q_in_expression_2(g2);
			if(Q1==Q2 && almost_short_table[is_almost_short(Q1*g1*Q1.inv())][is_almost_short(Q2*g2*Q2.inv())]){
				M->Contraction(i,i+1);
				return true;
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
		/*
			Q' ba^2ba Q, Q' a Q, Q' aba^2b Q
			Q' baba^2 Q, Q' a^2 Q, Q' a^2bab Q
		*/

		if(g2.len()%2==1 && (get_tau_in_expression_2(g2)==a || get_tau_in_expression_2(g2)==a2 ||
			                 get_tau_in_expression_2(g2)==b*a*b || get_tau_in_expression_2(g2)==b*a2*b)){
			Ga3b2 Q=Ga3b2(vector<string>(g2.e.begin(),g2.e.begin()+g2.len()/2)).inv();
			g1=Q.inv()*g1*Q; g2=Ga3b2({g2.e[g2.len()/2]}); g3=Q.inv()*g3*Q;

			if((g1==Ga3b2({"b","a^2","b","a"}) && g2==a && g3==Ga3b2({"a","b","a^2","b"})) ||
			   (g1==Ga3b2({"b","a","b","a^2"}) && g2==a2 && g3==Ga3b2({"a^2","b","a","b"}))){
				M->Contraction(i,i+2);
			   	return true;
			}
		}
	}
	return false;
}

/*--Operation 2--------------*/
bool OperationC(CR<Ga3b2>* M){
	int n=M->h.len();
	for(int i=1;i<n;i++){
		if(M->h.e[i-1].len()>=1 && M->h.e[i].len()==0){
			for(int j=i;j>=1;j--)
				M->Elementary_transformation(j,1);
			return true;
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
	for(int i=1;i<=M.h.len();i++) myassert(M.h.e[i-1].len()==0,"each component of the tuple after induction is 1");
	/***********/
	cout<<"| +----after the induction---\n";
	cout<<"| | M.h="<<M.h<<"\n";
	cout<<"| | M.H="<<M.H<<"\n";
	cout<<"| +=======\n";

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


/*

a^2b  a^2b  a^2ba^2baba  ba^2ba^2ba^2b  baba^2  a^2bab  a^2ba^2babababa

a^2b  a^2b  a^2ba^2baba  ba^2ba^2ba^2b  baba^2  a^2bab  a^2ba^2babababa

ba^2ba^2ba^2b  baba^2

baba^2     aba^2b ba^2ba^2ba^2b baba^2=ababa

*/
