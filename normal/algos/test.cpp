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

int test_normalize_short(int testnum,int n,int pa,int qa,int nb,int p,int q,bool details){
	string str1="+---------------------+";
	string str2="|  Test Case "+to_string(testnum);
	while(str2.size()+1<str1.size()) str2=str2+" "; str2=str2+"|";
	cout<<str1<<"\n"<<str2<<"\n"<<str1<<"\n";

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

	cout<<"+--------+\n";
	cout<<"| Step 0 |\n";
	cout<<"+--------+\n";
	
	Tuple<Ga3b2> g0(get_normal_form(pa,qa,nb,p,q).second.e),g,h,k,g_final;
	g=g0;
	list<vector<int>> F; F.clear();
	for(int j=1;j<=100;j++){
		int i=rand()%(n-1)+1;
		int epsilon=(rand()%2)*2-1;
		g.Elementary_transformation(i,epsilon);
	}
	cout<<"n, pa, qa, nb, p, q = "<<n<<", "<<pa<<", "<<qa<<", "<<nb<<", "<<p<<", "<<q<<"\n";
	cout<<"g="<<g<<"\n";

	cout<<"+--------+\n";
	cout<<"| Step 1 |\n";
	cout<<"+--------+\n";
	auto ret2=shorten_induction(g,details,"g");
	h=ret2.first; F.insert(F.end(),ret2.second.begin(),ret2.second.end());

	cout<<"+--------+\n";
	cout<<"| Step 2 |\n";
	cout<<"+--------+\n";
	auto ret3=transform_into_inverse_free(h,details,"h");
	k=ret3.first; F.insert(F.end(),ret3.second.begin(),ret3.second.end());
	
	cout<<"+--------+\n";
	cout<<"| Step 3 |\n";
	cout<<"+--------+\n";
	auto ret4=normalize_inverse_free_tuple(k,details,"k");
	F.insert(F.end(),ret4.begin(),ret4.end());

	cout<<"+--------+\n";
	cout<<"| Step 4 |\n";
	cout<<"+--------+\n";
	CR<Ga3b2> M; M.init(g); M.Apply(F); careful_restorations(&M,true);
	auto ret5=sort_concatenation(M.h,details);
	g_final=ret5.first;

	myassert(g0==g_final,"g0 == g_final");

	cout<<"\n";
	return ret5.second;
}

void multi_test_normalize_short(bool details){
	cout<<"This is a quick multi-test for short.cpp\n";
	cout<<"Input   nl   nr   s.t.   (g1,...,gn) with nl<=n<=nr   is an n-tuple in PSL(2,Z)=Ga3b2 satisfying g1***gn=1.\n";
	int n,nl,nr; cin>>nl>>nr;

	cout<<"Input t for the number of test cases.\n";
	int test; cin>>test;
	int max_length_exceptional_part=0;

	for(int t=1;t<=test;t++){
		n=rand()%(nr-nl+1)+nl;
		int pa,qa,nb,p,q;
		for(;;){
			pa=0,qa=0,nb=0,p=0,q=0;
			for(int i=1;i<=n;i++){
				int o=rand()%10;
				if(o==0) ++pa;
				if(o==1) ++qa;
				if(o==2) ++nb;
				if(o==3) ++p;
				if(o>=4) ++q; // I prefer more components conjugate to ba=t0
			}
			if(get_normal_form(pa,qa,nb,p,q).first) break;
		}
		int length_exceptional_part=test_normalize_short(t,n,pa,qa,nb,p,q,details);
		if(length_exceptional_part>max_length_exceptional_part) max_length_exceptional_part=length_exceptional_part;
	}

	cout<<"max_length_exceptional_part = "<<max_length_exceptional_part<<"\n";
}
