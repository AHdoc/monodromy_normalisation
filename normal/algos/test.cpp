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

void test_solve_short(){
	cout<<"This is a quick test for solve_short.\n";
	cout<<"Sample case:\n";
	cout<<"    28\n";
	cout<<"    a^2b ba^2 a^2b ba^2 a^2b ba^2 a^2b ba^2 a^2b ba^2 a^2b ba^2 a^2 a^2b ba^2 a^2b ba a^2b ba a a^2 b b b b a^2 a^2 a^2\n";

	int n; cin>>n;
	vector<Ga3b2> gvec; gvec.clear();
	Ga3b2 prod; prod.be_identity();
	for(int i=0;i<n;i++){
		Ga3b2 gveci; cin>>gveci;
		prod=prod*gveci;
		gvec.push_back(gveci);
	}
	Tuple<Ga3b2> g; g.init(gvec);
	cout<<"prod="<<prod<<"\n";
	
	for(int t=1;t<=100;t++){
		int i=rand()%(n-1)+1;
		int epsilon=(rand()%2)*2-1;
		g.Elementary_transformation(i,epsilon);
	}
	cout<<"g="<<g<<"\n";

	solve_short(g);
}