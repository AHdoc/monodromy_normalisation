void test_Ga3b2(){
	cout<<"This is a quick test for struct Ga3b2.\n";
	cout<<"Input two elements in Ga3b2, for instance, aba^2 a^2ba^2.\n";

	Ga3b2 g, h;
	cin>>g>>h;
	cout<<"g="<<g<<"\n";
	cout<<"h="<<h<<"\n";
	cout<<"g.inv()="<<g.inv()<<"\n";
	cout<<"h.inv()="<<h.inv()<<"\n";
	cout<<"g*h="<<g*h<<"\n";
	cout<<"pow(g,2)*pow(h,3)="<<pow(g,2)*pow(h,3)<<"\n";
	g.conjugation(h);
	cout<<"after g.conjugation(h), g="<<g<<"\n";
}

void test_tuple(){
	cout<<"This is a quick test for structs Tuple and iteratedTuple.\n";
	
	int n; cin>>n;
	vector<Ga3b2> gvec; gvec.clear();
	for(int i=0;i<n;i++){
		Ga3b2 gveci; cin>>gveci;
		gvec.push_back(gveci);
	}
	Tuple<Ga3b2> g; g.init(gvec);
	cout<<g<<"\n";

	iteratedTuple<Ga3b2> it_g;
	it_g.init(g);
	cout<<"it_g="<<it_g<<"\n";
	cerr<<"it_g.h()="<<it_g.h()<<",it_g.ev()="<<it_g.ev()<<"\n";
}

void test_CR(){
	cout<<"This is a quick test for elementary transformations and contractions.\n";
	cout<<"Sample case:\n";
	cout<<"    20\n";
	cout<<"    a b a^2 b aba aba^2 a^2b a a a^2\n";
	cout<<"    a b a^2 b aba aba^2 a^2b a a a^2\n";
	cout<<"    8\n";
	cout<<"    1 2 1\n";
	cout<<"    1 3 -1\n";
	cout<<"    2 2 4\n";
	cout<<"    1 2 1\n";
	cout<<"    2 7 11\n";
	cout<<"    2 2 5\n";
	cout<<"    1 2 1\n";
	cout<<"    1 3 1\n";

	int n; cin>>n;
	vector<Ga3b2> gvec; gvec.clear();
	for(int i=0;i<n;i++){
		Ga3b2 gveci; cin>>gveci;
		gvec.push_back(gveci);
	}
	Tuple<Ga3b2> g; g.init(gvec);

	CR<Ga3b2> M;
	M.init(g);
	cout<<"h="<<M.h<<"\n";
	cout<<"H="<<M.H<<"\n";

	int m; cin>>m;
	for(int i=0;i<m;i++){
		int o; cin>>o; // 1 - elementary trasnformation; 2 - contraction.
		if(o==1){
			int j,epsilon; cin>>j>>epsilon;
			M.Elementary_transformation(j,epsilon);
		}else if(o==2){
			int l,r; cin>>l>>r;
			M.Contraction(l,r);
		}
		cout<<"h="<<M.h<<"\n";
		cout<<"H="<<M.H<<"\n";
	}
}

Tuple<Ga3b2> test_normalize_short(int n,bool details){
	if(details){
		cout<<"    - Step 0: generate a random tuple g=(g1,...,gn) of length n=pa+qa+nb+p+q="<<n<<" s.t.\n";
		cout<<"        - pa of them are conjugates of a;\n";
		cout<<"        - qa of them are conjugates of a^2;\n";
		cout<<"        - nb of them are conjugates of b;\n";
		cout<<"        - p of them are conjugates of aba;\n";
		cout<<"        - q of them are conjugates of a^2ba^2.\n";
		cout<<"    - Step 1: transform (g1,...,gn) into (h1,...,hm) * several pairs (x,x^{-1}) and triples (l,l,l) with l^3=1\n";
		cout<<"        - where (h1,...,hm) is a tuple of short elements--\n";
		cout<<"    - Step 2: transform/contract (h1,...,hm) into (k1,...,kl)\n";
		cout<<"        - where (k1,...,kl) is an inverse-free tuple of short elements containing <=1 component equal to a/a^2/b\n";
		cout<<"    - Step 3: contract/normalize (k1,...,kl)\n";
		cout<<"        - into (s0,s2,s0,s2,s0,s2)^? * (t0,t2,t0,t2,t0,t2)^? * (s_i,t_i)^?\n";
		cout<<"    - Step 4: write the resulting tuple as a desired concatenation\n";
		cout<<"        - Step 4.1: restore the resulting tuple\n";
		cout<<"        - Step 4.2: handle the exceptional part, which is a tuple of finite length\n";
		cout<<"\n";
	}

	cout<<"+--------+\n";
	cout<<"| Step 0 |\n";
	cout<<"+--------+\n";
	
	auto ret0=random_tuple_of_short_elements(n,true);
	Tuple<Ga3b2> g0=ret0.first,g=ret0.second,h,k,g_final;
	list<vector<int>> F; F.clear();
	cout<<"g="<<g<<"\n";

	cout<<"+--------+\n";
	cout<<"| Step 1 |\n";
	cout<<"+--------+\n";
	auto ret1=shorten_induction(g,details,"g");
	h=ret1.first; F.insert(F.end(),ret1.second.begin(),ret1.second.end());

	cout<<"+--------+\n";
	cout<<"| Step 2 |\n";
	cout<<"+--------+\n";
	auto ret2=transform_into_inverse_free(h,details,"h");
	k=ret2.first; F.insert(F.end(),ret2.second.begin(),ret2.second.end());
	
	cout<<"+--------+\n";
	cout<<"| Step 3 |\n";
	cout<<"+--------+\n";
	auto ret3=normalize_inverse_free_tuple(k,details,"k");
	F.insert(F.end(),ret3.begin(),ret3.end());

	cout<<"+--------+\n";
	cout<<"| Step 4 |\n";
	cout<<"+--------+\n";
	CR<Ga3b2> M; M.init(g); M.Apply(F); careful_restorations(&M,details);
	auto ret4=sort_concatenation(M.h,details);
	g_final=ret4.first;

	myassert(g0==g_final,"g0 == g_final");

	cout<<"\n";
	return ret4.second;
}

void multi_test_normalize_short(bool details){
	cout<<"This is a quick multi-test for short.cpp\n";
	cout<<"Input   nl   nr   s.t.   (g1,...,gn) with nl<=n<=nr   is an n-tuple in PSL(2,Z)=Ga3b2 satisfying g1***gn=1.\n";
	int n,nl,nr; cin>>nl>>nr;

	cout<<"Input t for the number of test cases.\n";
	int test; cin>>test;
	Tuple<Ga3b2> longest_exceptional_part(vector<Ga3b2>{});

	for(int t=1;t<=test;t++){
		string str1="+---------------------+";
		string str2="|  Test Case "+to_string(t);
		while(str2.size()+1<str1.size()) str2=str2+" "; str2=str2+"|";
		cout<<str1<<"\n"<<str2<<"\n"<<str1<<"\n";

		n=rand()%(nr-nl+1)+nl;
		Tuple<Ga3b2> g_excep=test_normalize_short(n,details);
		if(g_excep.len()>longest_exceptional_part.len()) longest_exceptional_part=g_excep;

		cout<<"longest_exceptional_part = "<<longest_exceptional_part<<" of length "<<longest_exceptional_part.len()<<"\n";
	}

	cout<<"longest_exceptional_part = "<<longest_exceptional_part<<" of length "<<longest_exceptional_part.len()<<"\n";
}

/******************************/

void test_normalize_almost_short(int n,bool details){
	Tuple<Ga3b2> g=random_tuple_of_almost_short_elements(n);
	almost_shorten_induction(g,details,"g");

	// almost_shorten_induction(Tuple<Ga3b2>(
	// 	{a*b, a2*b, a2*b*a, b*a*b*a*b*a2*b*a*b*a2*b*a2*b, b*a2*b*a2*b*a2*b*a2*b*a*b,
	// 	 a*b*a*b*a2, a*b*a*b*a, a2*b*a*b*a, a*b*a2*b*a2*b*a2*b*a2, a*b*a*b*a2*b*a*b*a*b*a2*b*a*b*a2*b*a2,
	// 	 a2*b*a2*b*a2, b*a2*b*a, b*a*b*a2}),true,"g");
}

void multi_test_normalize_almost_short(bool details){
	cout<<"This is a quick multi-test for almost_short.cpp\n";
	cout<<"Input   nl   nr   s.t.   (g1,...,gn) with nl<=n<=nr   is an n-tuple in PSL(2,Z)=Ga3b2 of almost short elements satisfying g1***gn=1.\n";
	int n,nl,nr; cin>>nl>>nr;

	cout<<"Input t for the number of test cases.\n";
	int test; cin>>test;

	for(int t=1;t<=test;t++){
		string str1="+---------------------+";
		string str2="|  Test Case "+to_string(t);
		while(str2.size()+1<str1.size()) str2=str2+" "; str2=str2+"|";
		cout<<str1<<"\n"<<str2<<"\n"<<str1<<"\n";

		n=rand()%(nr-nl+1)+nl;
		test_normalize_almost_short(n,details);
	}
}
