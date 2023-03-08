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
			myassert(Ga3b2({g2.e[g2.len()/2]})==Q*g2*Q.inv(),"g2=Q^{-1}aQ or Q^{-1}a^2Q");
			g1=Q*g1*Q.inv(); g2=Q*g2*Q.inv(); g3=Q*g3*Q.inv();

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

/*--Elementary transformation decreases S2_complexity--------------*/
bool OperationD(CR<Ga3b2>* M,int il,int ir){ // restrict that 1<=il<=i<=ir<=n-1
	for(int smallest_i=il;smallest_i<=ir;smallest_i++){
		for(int i=smallest_i;i<=smallest_i+1 && i<=ir;i++) for(int epsilon=-1;epsilon<=1;epsilon+=2){
			Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i,epsilon);
			if(S2_complexity(M->h)>S2_complexity(h1)){
				M->Elementary_transformation(i,epsilon);
				return true;
			}
		}
	}
	for(int smallest_i=il;smallest_i<=ir;smallest_i++){
		for(int i1=smallest_i;i1<=smallest_i+1 && i1<=ir;i1++) for(int epsilon1=-1;epsilon1<=1;epsilon1+=2)
		for(int i2=smallest_i;i2<=smallest_i+1 && i2<=ir;i2++) for(int epsilon2=-1;epsilon2<=1;epsilon2+=2){
			Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i1,epsilon1);
			Tuple<Ga3b2> h2=h1  ; h2.Elementary_transformation(i2,epsilon2);
			//if(S2_complexity(M->h)==S2_complexity(h1) && S2_complexity(h1)>S2_complexity(h2)){
			if(S2_complexity(M->h)>S2_complexity(h2)){
				M->Elementary_transformation(i1,epsilon1);
				M->Elementary_transformation(i2,epsilon2);
				return true;
			}
		}
	}
	for(int smallest_i=il;smallest_i<=ir;smallest_i++){
		for(int i1=smallest_i;i1<=smallest_i+1 && i1<=ir;i1++) for(int epsilon1=-1;epsilon1<=1;epsilon1+=2)
		for(int i2=smallest_i;i2<=smallest_i+1 && i2<=ir;i2++) for(int epsilon2=-1;epsilon2<=1;epsilon2+=2)
		for(int i3=smallest_i;i3<=smallest_i+1 && i3<=ir;i3++) for(int epsilon3=-1;epsilon3<=1;epsilon3+=2){
			Tuple<Ga3b2> h1=M->h; h1.Elementary_transformation(i1,epsilon1);
			Tuple<Ga3b2> h2=h1  ; h2.Elementary_transformation(i2,epsilon2);
			Tuple<Ga3b2> h3=h2  ; h3.Elementary_transformation(i3,epsilon3);
			//if(S2_complexity(M->h)==S2_complexity(h1) && S2_complexity(h1)==S2_complexity(h2) && S2_complexity(h2)>S2_complexity(h3)){
			if(S2_complexity(M->h)>S2_complexity(h3)){
				M->Elementary_transformation(i1,epsilon1);
				M->Elementary_transformation(i2,epsilon2);
				M->Elementary_transformation(i3,epsilon3);
				return true;
			}
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

void careful_restorations_for_almost_short_conjugates(CR<Ga3b2>* M,bool details){
	if(details) cout<<"~~~careful_restorations_for_almost_short_conjugates M~~~\n";
	if(details) cout<<"h = "<<M->h<<"\n";
	if(details) cout<<"H = "<<M->H<<"\n";
	for(;;){
		pair<int,int> sg=M->Restoration(); // segment [l,r]
		int l=sg.first,r=sg.second;

		if(l==-1) break;

		if(l+1==r && M->h.e[l-1]==M->h.e[r-1].inv()){
			//cout<<"("<<M->h.e[l-1]<<","<<M->h.e[r-1]<<") ---> ";
			Ga3b2 Q=get_Q_in_expression_2(M->h.e[l-1]);
			for(int _=l;_<=r;_++){
				M->h.e[_-1].conjugation(Q.inv());
				M->H.e[_-1]->conjugation(Q.inv());
			}
			//cout<<"("<<M->h.e[l-1]<<","<<M->h.e[r-1]<<") \n";
		}else{
			//for(int _=l;_<=r;_++){
				//if(_==l) cout<<"(";
				//cout<<M->h.e[_-1];
				//if(_<r) cout<<","; else cout<<") ---> ";
			//}
			while(OperationD(M,l,r-1)) continue;
			//for(int _=l;_<=r;_++){
			//	if(_==l) cout<<"(";
			//	cout<<M->h.e[_-1];
			//	if(_<r) cout<<","; else cout<<")\n";
			//}
		}

		if(details) cout<<"--> "<<M->h<<"\n";
		if(details) cout<<"    "<<M->H<<"\n";
	}
	if(details) cout<<"~~~~~~\n";
}

/******************/

/*--an induction that--------------*/
/*----transform (g1,...,gn) in Ga3b2 whose components are conjuagates of almost short elements--*/
/*--into--*/
/*----(h1,...,hm) * g_non_inverse_free--*/
/*--s.t.--*/
/*----each component of (h1,...,hm) is almost short--*/
/*--return mp(h,F)--*/
void almost_shorten_induction(Tuple<Ga3b2> g_input,bool details,string namestr){
	cout<<"| +====start the almost_shorten_induction on "+namestr+"===\n";
	cout<<"| | "<<namestr<<"="<<g_input<<"\n";

	CR<Ga3b2> M;
	Tuple<Ga3b2> g_non_inverse_free;

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
	for(int i=1;i<=M.h.len();i++) myassert(M.h.e[i-1].len()==0,"each component of the tuple after induction is 1");
	/***********/
	
	careful_restorations_for_almost_short_conjugates(&M,false);
	
	// Fron now on, each component of M.h is almost short;
	myassert(each_component_is_almost_short(M.h),"fron now on, each component of M.h is almost short");
	cout<<"| +----a tuple of almost short elements---\n";
	cout<<"| | M.h="<<M.h<<"\n";
	cout<<"| +=======\n";
}

/***************************/
