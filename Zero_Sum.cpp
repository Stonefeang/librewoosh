vi exact_sum(vi wek, int cel)
{
	int n=wek.size();
	n++;
	vi inv(n);
	inv[1]=1;
	for (int i=2; i<n; i++)
		inv[i]=inv[n%i]*(ll)(n-n/i)%n;
	int wsk=1;
	vi kto(n, -1);
	for (int i=0; i<n-1; i++)
	{
		while(kto[wsk]>=0)
			wsk++;
		int g=wsk*(ll)inv[wek[i]%n]%n;
		int bsa=1;
		int bsb=g;
		while(bsa<bsb)
		{
			int bss=(bsa+bsb)>>1;
			if (kto[(bss*(ll)wek[i])%n]==-1)
				bsb=bss;
			else
				bsa=bss+1;
		}
		kto[(bsa*(ll)wek[i])%n]=i;
	}
	vi ret;
	while(cel)
	{
		ret.push_back(kto[cel]);
		cel=(cel-wek[kto[cel]]%n+n)%n;
	}
	sort(ret.begin(), ret.end());
	return ret;
}

vi n_zero_sum(vi wek)
{
	int n=wek.size();
	n=(n+1)/2;
	if (n==1)
		return {0};
	int g=2;
	while(n%g)
		g++;
	if (g<n)
	{
		vi uzyte(2*n-1);
		vector<vi> podz;
		vi sumy;
		int wsk=0;
		vi dost;
		while(wsk<2*n-1)
		{
			while((int)dost.size()<2*g-1)
				dost.push_back(wsk++);
			vi daj;
			for (int i : dost)
				daj.push_back(wek[i]);
			vi wez=n_zero_sum(daj);
			vi nowy;
			int s=0;
			for (int i : wez)
			{
				uzyte[dost[i]]=1;
				nowy.push_back(dost[i]);
				s=(s+daj[i])%n;
			}
			podz.push_back(nowy);
			sumy.push_back(s/g);
			wez.clear();
			for (int i : dost)
				if (!uzyte[i])
					wez.push_back(i);
			dost=wez;
		}
		for (int i : n_zero_sum(sumy))
			for (int j : podz[i])
				uzyte[j]=2;
		vi ret;
		for (int i=0; i<2*n-1; i++)
			if (uzyte[i]==2)
				ret.push_back(i);
		return ret;
	}
	vector<pii> pos;
	for (int i=0; i<2*n-1; i++)
		pos.push_back({wek[i]%n, i});
	sort(pos.begin(), pos.end());
	int baza=pos.back().first;
	vi roznice;
	for (int i=0; i<n-1; i++)
	{
		if (pos[i].first==pos[i+n-1].first)
		{
			vi ret;
			for (int j=0; j<n; j++)
				ret.push_back(pos[i+j].second);
			return ret;
		}
		baza=(baza+pos[i].first)%n;
		roznice.push_back((pos[i+n-1].first-pos[i].first+n)%n);
	}
	vi ktore(n-1);
	for (int i : exact_sum(roznice, (n-baza)%n))
		ktore[i]=1;
	vi ret{pos.back().second};
	for (int i=0; i<n-1; i++)
		ret.push_back(pos[i+(n-1)*ktore[i]].second);
	sort(ret.begin(), ret.end());
	return ret;
}
