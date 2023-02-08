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
	cout<<"    1 5 5 14 3\n"; // -> 1 2 1 11 0 -> 0 1 1 5 0
	int pa,qa,nb,p,q;
	cin>>pa>>qa>>nb>>p>>q;
	vector<Ga3b2> gvec; gvec.clear();
	Ga3b2 a({"a"}),a2({"a^2"}),b({"b"}),s0({"a^2","b"}),s2({"b","a^2"}),t0({"b","a"}),t2({"a","b"});
	while(pa>=3){
		gvec.push_back(a); gvec.push_back(a); gvec.push_back(a);
		pa-=3;
	}
	while(qa>=3){
		gvec.push_back(a2); gvec.push_back(a2); gvec.push_back(a2);
		qa-=3;
	}
	while(pa>=1 && qa>=1){
		gvec.push_back(a); gvec.push_back(a2);
		pa--; qa--;
	}
	while(nb>=2){
		gvec.push_back(b); gvec.push_back(b);
		nb-=2;
	}
	while(p>=1 && q>=1){
		gvec.push_back(s0); gvec.push_back(t0);
		p--; q--;
	}
	while(p>=6){
		gvec.push_back(s0); gvec.push_back(s2); gvec.push_back(s0);
		gvec.push_back(s2); gvec.push_back(s0); gvec.push_back(s2);
		p-=6;
	}
	while(q>=6){ // ba ab ba ab ba ab = 1
		gvec.push_back(t0); gvec.push_back(t2); gvec.push_back(t0);
		gvec.push_back(t2); gvec.push_back(t0); gvec.push_back(t2);
		q-=6;
	}
	if(qa==1 && p==2){
		gvec.push_back(a2); gvec.push_back(s0); gvec.push_back(s2);
		qa--; p-=2;
	}else if(pa==1 && q==2){
		gvec.push_back(a); gvec.push_back(t2); gvec.push_back(t0);
		pa--; q-=2;
	}else if(pa==1 && p==4){
		gvec.push_back(a); gvec.push_back(s0); gvec.push_back(s0); gvec.push_back(s2); gvec.push_back(s0);
		pa--; p-=4;
	}else if(qa==1 && q==4){
		gvec.push_back(a2); gvec.push_back(t0); gvec.push_back(t2); gvec.push_back(t0); gvec.push_back(t0);
		qa--; q-=4;
	}else if(nb==1 && p==3){ // b aab baa aab
		gvec.push_back(b); gvec.push_back(s0); gvec.push_back(s2); gvec.push_back(s0);
		nb--; p-=3;
	}else if(nb==1 && q==3){ // b ba ab ba
		gvec.push_back(b); gvec.push_back(t0); gvec.push_back(t2); gvec.push_back(t0);
		nb--; q-=3;
	}else if(pa==1 && nb==1 && p==1){
		gvec.push_back(a); gvec.push_back(b); gvec.push_back(s2);
		pa--; nb--; p--;
	}else if(qa==1 && nb==1 && q==1){
		gvec.push_back(a2); gvec.push_back(b); gvec.push_back(t0);
		qa--; nb--; q--;
	}else if(pa==1 && nb==1 && q==5){
		gvec.push_back(a); gvec.push_back(t2); gvec.push_back(t0);
		gvec.push_back(b); gvec.push_back(t0); gvec.push_back(t2); gvec.push_back(t0);
		pa--; nb--; q-=5;
	}else if(qa==1 && nb==1 && p==5){
		gvec.push_back(a2); gvec.push_back(s0); gvec.push_back(s2);
		gvec.push_back(b); gvec.push_back(s0); gvec.push_back(s2); gvec.push_back(s0);
		qa--; nb--; p-=5;
	}
	if(pa!=0 || qa!=0 || nb!=0 || p!=0 || q!=0){
		cout<<"Bad (pq,qa,nb,p,q).\n";
		cout<<"-> "<<pa<<" "<<qa<<" "<<nb<<" "<<p<<" "<<q<<"\n";
		return;
	}

	Ga3b2 prod; prod.be_identity();
	for(int i=0;i<gvec.size();i++){
		Ga3b2 gveci=gvec[i];
		prod=prod*gveci;
	}
	myassert(prod.len()==0,"prod=1");

	/******************************************************************/

	cout<<"Input t for the number of test cases.\n";
	int test;
	cin>>test;

	for(int t=1;t<=test;t++){
		Tuple<Ga3b2> g; g.init(gvec);
		int n=g.len();

		string s1="+---------------------+";
		string s2="|  Test Case "+to_string(t);
		while(s2.size()+1<s1.size()) s2=s2+" "; s2=s2+"|";
		cout<<s1<<"\n"<<s2<<"\n"<<s1<<"\n";

		cout<<"    - Step 1: generate a random tuple g=(g1,...,gn) of length n=pa+pq+nb+p+q="<<n<<"\n";
		cout<<"    - Step 2: transform (g1,...,gn) into (h1,...,hm) bullet g_non_inverse_free\n";
		cout<<"              where (h1,...,hm) is obtained from an induction whose components are short,\n";
		cout<<"                    g_non_inverse_free consists of (g,g^{-1}) and (l,l,l) with l^3=1,\n";
		cout<<"                    (h1,...,hm) is inverse-free.\n";
		cout<<"    - Step 3: normalize (h1,...,hm)\n";
		
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
		shorten_induction(g);

		cout<<"+--------+\n";
		cout<<"| Step 3 |\n";
		cout<<"+--------+\n";

	}
}