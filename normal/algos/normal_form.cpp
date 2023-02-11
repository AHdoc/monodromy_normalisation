pair<bool,Tuple<Ga3b2>> get_normal_form(int pa,int qa,int nb,int p,int q){
	Ga3b2 a({"a"}),a2({"a^2"}),b({"b"});
	Ga3b2 s0({"a^2","b"}),s1({"a","b","a"}),s2({"b","a^2"}),t0({"b","a"}),t1({"a^2","b","a^2"}),t2({"a","b"});

	Tuple<Ga3b2> g0; g0.clear();
	while(pa>=3){
		g0=bullet(g0,Tuple<Ga3b2>({a,a,a}));
		pa-=3;
	}
	while(qa>=3){
		g0=bullet(g0,Tuple<Ga3b2>({a2,a2,a2}));
		qa-=3;
	}
	while(pa>=1 && qa>=1){
		g0=bullet(g0,Tuple<Ga3b2>({a,a2}));
		pa--; qa--;
	}
	while(nb>=2){
		g0=bullet(g0,Tuple<Ga3b2>({b,b}));
		nb-=2;
	}
	while(p>=1 && q>=1){
		g0=bullet(g0,Tuple<Ga3b2>({s0,t0}));
		p--; q--;
	}
	while(p>=6){ // aab baa aab baa aab baa = 1
		g0=bullet(g0,Tuple<Ga3b2>({s0,s2,s0,s2,s0,s2}));
		p-=6;
	}
	while(q>=6){ // ba ab ba ab ba ab = 1
		g0=bullet(g0,Tuple<Ga3b2>({t0,t2,t0,t2,t0,t2}));
		q-=6;
	}
	if(qa==1 && p==2){
		g0=bullet(g0,Tuple<Ga3b2>({a2,s0,s2}));
		qa--; p-=2;
	}else if(pa==1 && q==2){
		g0=bullet(g0,Tuple<Ga3b2>({a,t2,t0}));
		pa--; q-=2;
	}else if(pa==1 && p==4){
		g0=bullet(g0,Tuple<Ga3b2>({a,s0,s0,s2,s0}));
		pa--; p-=4;
	}else if(qa==1 && q==4){
		g0=bullet(g0,Tuple<Ga3b2>({a2,t0,t2,t0,t0}));
		qa--; q-=4;
	}else if(nb==1 && p==3){ // b aab baa aab
		g0=bullet(g0,Tuple<Ga3b2>({b,s0,s2,s0}));
		nb--; p-=3;
	}else if(nb==1 && q==3){ // b ba ab ba
		g0=bullet(g0,Tuple<Ga3b2>({b,t0,t2,t0}));
		nb--; q-=3;
	}else if(pa==1 && nb==1 && p==1){
		g0=bullet(g0,Tuple<Ga3b2>({a,b,s2}));
		pa--; nb--; p--;
	}else if(qa==1 && nb==1 && q==1){
		g0=bullet(g0,Tuple<Ga3b2>({a2,b,t0}));
		qa--; nb--; q--;
	}else if(pa==1 && nb==1 && q==5){
		g0=bullet(g0,Tuple<Ga3b2>({a,t2,t0,b,t0,t2,t0}));
		pa--; nb--; q-=5;
	}else if(qa==1 && nb==1 && p==5){
		g0=bullet(g0,Tuple<Ga3b2>({a2,s0,s2,b,s0,s2,s0}));
		qa--; nb--; p-=5;
	}else if(pa==2 && p==2){
		g0=bullet(g0,Tuple<Ga3b2>({a,a,s0,s2}));
		pa-=2; p-=2;
	}else if(qa==2 && q==2){
		g0=bullet(g0,Tuple<Ga3b2>({a2,a2,t2,t0}));
		qa-=2; q-=2;
	}else if(qa==2 && p==4){
		g0=bullet(g0,Tuple<Ga3b2>({a2,a2,s0,s0,s2,s0}));
		qa-=2; p-=4;
	}else if(pa==2 && q==4){
		g0=bullet(g0,Tuple<Ga3b2>({a,a,t0,t2,t0,t0}));
		pa-=2; q-=4;
	}else if(qa==2 && nb==1 && p==1){
		g0=bullet(g0,Tuple<Ga3b2>({a2,a2,b,s2}));
		qa-=2; nb--; p--;
	}else if(pa==2 && nb==1 && q==1){
		g0=bullet(g0,Tuple<Ga3b2>({a,a,b,t0}));
		pa-=2; nb--; q--;
	}else if(qa==2 && nb==1 && q==5){
		g0=bullet(g0,Tuple<Ga3b2>({a2,a2,t2,t0,b,t0,t2,t0}));
		qa-=2; nb--; q-=5;
	}else if(pa==2 && nb==1 && p==5){
		g0=bullet(g0,Tuple<Ga3b2>({a,a,s0,s2,b,s0,s2,s0}));
		pa-=2; nb--; p-=5;
	}
	
	if(pa!=0 || qa!=0 || nb!=0 || p!=0 || q!=0){
		return make_pair(false,g0);
		//cout<<"Bad (pq,qa,nb,p,q).\n";
		//cout<<"-> "<<pa<<" "<<qa<<" "<<nb<<" "<<p<<" "<<q<<"\n";
	}

	Ga3b2 prod; prod.be_identity();
	for(Ga3b2 x:g0.e) prod=prod*x;
	myassert(prod.len()==0,"prod=1");
	
	return make_pair(true,g0);
}
