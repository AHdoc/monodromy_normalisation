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

void test_solve_short(bool details){
	cout<<"This is a quick test for solve_short.\n";

	cout<<"Input   pa   qa   nb   p   q   such that\n";
	cout<<"    - (g1,...,gn) is an n-tuple in PSL(2,Z)=Ga3b2 satisfying g1***gn=1;\n";
	cout<<"    - n = pa+qa+nb+p+q;\n";
	cout<<"    - g1, ..., gn are conjugates of a, a^2, b, aba, a^2ba^2;\n";
	cout<<"    - pa of them are conjugates of a;\n";
	cout<<"    - qa of them are conjugates of a^2;\n";
	cout<<"    - nb of them are conjugates of b;\n";
	cout<<"    - p of them are conjugates of aba;\n";
	cout<<"    - q of them are conjugates of a^2ba^2.\n";
	cout<<"Sample case:\n";
	cout<<"    3 5 7 14 7\n"; // -> 0 2 1 7 0
	int pa,qa,nb,p,q;
	Tuple<Ga3b2> g0;
	for(;;){
		cin>>pa>>qa>>nb>>p>>q;
		auto ret=get_normal_form(pa,qa,nb,p,q);
		if(ret.first){
			g0=ret.second;
			break;
		}else{
			cout<<"Impossible pattern! Please input again.\n";
		}
	}
	/******************************************************************/

	cout<<"Input t for the number of test cases.\n";
	int test;
	cin>>test;

	for(int t=1;t<=test;t++){
		Tuple<Ga3b2> g(g0.e),h;
		list<vector<int>> F1,F2;
		int n=g.len();

		string str1="+---------------------+";
		string str2="|  Test Case "+to_string(t);
		while(str2.size()+1<str1.size()) str2=str2+" "; str2=str2+"|";
		cout<<str1<<"\n"<<str2<<"\n"<<str1<<"\n";

		cout<<"    - Step 1: generate a random tuple g=(g1,...,gn) of length n=pa+pq+nb+p+q="<<n<<"\n";
		cout<<"    - Step 2: shorten g=(g1,...,gn) and make it inverse-free#\n";
		cout<<"        transform/contract (g1,...,gn) into (h1,...,hm) * (y,y^{-1}) * (x,x^{-1}) * ... * (l,l,l) * ..., where:\n";
		cout<<"        - (h1,...,hm) and (y,y^{-1}) are iterated tuples of short elements\n";
		cout<<"        - one of h and k is of height 3, the other is of height 1\n";
		cout<<"        - (h1,...,hm) is inverse-free and contains at most 1 component being a, a^2 or b\n";
		cout<<"        - all pairs of the form (x,x^{-1}) and triples of the form (l,l,l) with l^3=1 are of height 1\n";
		cout<<"    - Step 3: normalize (h1,...hm)\n";

		cout<<"\n";

		cout<<"+--------+\n";
		cout<<"| Step 1 |\n";
		cout<<"+--------+\n";
		for(int j=1;j<=100;j++){
			int i=rand()%(n-1)+1;
			int epsilon=(rand()%2)*2-1;
			g.Elementary_transformation(i,epsilon);
		}
		cout<<"g="<<g<<"\n";

		cout<<"+--------+\n";
		cout<<"| Step 2 |\n";
		cout<<"+--------+\n";
		auto ret=transform_into_inverse_free(g,details);
		h=ret.first; F1=ret.second;
		
		cout<<"+--------+\n";
		cout<<"| Step 3 |\n";
		cout<<"+--------+\n";
		F2=normalize_inverse_free_tuple(h,details);

		cout<<"+--------+\n";
		cout<<"| Output |\n";
		cout<<"+--------+\n";
		CR<Ga3b2> M; M.init(g);
		for(auto it:F1) if(it[0]==1) M.Elementary_transformation(it[2],it[3]); else M.Contraction(it[1],it[2]);
		for(auto it:F2) if(it[0]==1) M.Elementary_transformation(it[2],it[3]); else M.Contraction(it[1],it[2]);
		for(;;){
			int k=M.Restoration();
			if(k==-1) break;
			while(Operation3(&M,k,k)) continue;
		}
		cout<<"g="<<M.h<<"\n";

		cout<<"\n";
	}
}