void random_locally_diagonal_conjugacy(int n,Tuple<Ga3b2>* g,int num){
	while(num--){
		int i,j;
		for(;;){
			i=rand()%n; j=i; Ga3b2 prod=g->e[i];
			while(prod.len()>0 && j+1<n) prod=prod*g->e[++j];
			if(prod.len()==0){
				int oq=rand()%3; Ga3b2 q;
				if(oq==0) q=a;
				if(oq==1) q=a2;
				if(oq==2) q=b;
				for(int k=i;k<=j;k++) g->e[k]=q.inv()*g->e[k]*q;
				break;
			}
		}
	}
}

void random_elementary_transformations(int n,Tuple<Ga3b2>* g,int max_total_l){
	for(int t=1;t<=max_total_l;t++){
		int total_l=0;
		for(int i=0;i<n;i++) total_l+=g->e[i].len();

		if(total_l>=max_total_l) continue;

		int i=rand()%(n-1)+1,epsilon=(rand()%2)*2-1;
		g->Elementary_transformation(i,epsilon);
	}
}

pair<Tuple<Ga3b2>,Tuple<Ga3b2>> random_tuple_of_short_elements(int n,bool s_prefer){
	Tuple<Ga3b2> g0;
	for(int i,pa,qa,nb,p,q;;){
		for(i=1,pa=0,qa=0,nb=0,p=0,q=0;i<=n;i++){
			int o=rand()%(5+(s_prefer?5:0));
			if(o==0) ++pa;
			if(o==1) ++qa;
			if(o==2) ++nb;
			if(o==3) ++q;
			if(o>=4) ++p; // I prefer more components conjugate to aab=s0
		}
		pair<bool,Tuple<Ga3b2>> ret=get_normal_form(pa,qa,nb,p,q);
		if(ret.first){
			cout<<"random a tuple of short elements with\n";
			cout<<"   n, pa, qa, nb, p, q = "<<n<<", "<<pa<<", "<<qa<<", "<<nb<<", "<<p<<", "<<q<<"\n";
			g0=ret.second;
			break;
		}
	}
	Tuple<Ga3b2> g=g0;
	int max_num_locally_diagonal_conjugacy=30;
	int max_total_l=900;
	random_locally_diagonal_conjugacy(n,&g,max_num_locally_diagonal_conjugacy);
	random_elementary_transformations(n,&g,max_total_l);

	return make_pair(g0,g);
}

/***************************************/

Tuple<Ga3b2> random_tuple_of_almost_short_elements(int n){
	Tuple<Ga3b2> g;
	for(int i,pa,qa,nb,p,q,n2,n3;;){
		for(i=1,pa=0,qa=0,nb=0,p=0,q=0,n2=0,n3=0;i<=n;i++){
			int o=rand()%16;
			if(o==0) ++pa;
			if(o==1) ++qa;
			if(o==2) ++nb;
			if(o==3) ++q;
			if(4<=o && o<=9) ++p;
			if(10<=o && o<=12) ++n2;
			if(13<=o && o<=15) ++n3;
		}
		pair<bool,Tuple<Ga3b2>> ret=get_normal_form(pa,qa,nb,p,q);
		if(ret.first){
			g=ret.second;
			for(int j=1;j<=n2;j++) g=bullet(g, Tuple<Ga3b2>({Ga3b2({"a^2","b","a","b"}), Ga3b2({"b","a^2","b","a"})}));
			for(int j=1;j<=n3;j++) g=bullet(g, Tuple<Ga3b2>({Ga3b2({"a^2","b","a","b"}), t0, s1}));
			break;
		}
	}
	
	int max_num_locally_diagonal_conjugacy=30;
	int max_total_l=900;
	random_locally_diagonal_conjugacy(n,&g,max_num_locally_diagonal_conjugacy);
	random_elementary_transformations(n,&g,max_total_l);

	return g;
}
