ll f(int x, int y)
{
	//Query about matrix entry (x-th row, y-th column)
}

vi smawk(const vi &rzo, vi kol)
{
	int r=rzo.size();
	int k=kol.size();
	if (r==1)
	{
		int g=-1;
		ll top=-2;
		for (int i=0; i<k; i++)
		{
			ll x=f(rzo[0], kol[i]);
			if (x>top)
			{
				top=x;
				g=i;
			}
		}
		return {g};
	}
	vi stos(k);
	iota(stos.begin(), stos.end(), 0);
	if (k>r)
	{
		stos.clear();
		for (int i=0; i<k; i++)
		{
			int s=stos.size();
			while(s && f(rzo[s-1], kol[i])>f(rzo[s-1], kol[stos.back()]))
			{
				stos.pop_back();
				s--;
			}
			if (s<r)
				stos.push_back(i);
		}
		vi pras(stos.size());
		for (int i=0; i<(int)stos.size(); i++)
			pras[i]=kol[stos[i]];
		kol=pras;
	}
	vi cod;
	for (int i=0; i<r; i+=2)
		cod.push_back(rzo[i]);
	vi wez=smawk(cod, kol);
	vi ret(r);
	for (int i=0; i<r; i+=2)
		ret[i]=wez[i>>1];
	for (int i=1; i<r; i+=2)
	{
		int a=ret[i-1];
		int b=k-1;
		if (i+1<r)
			b=ret[i+1];
		int g=-1;
		ll top=-2;
		for (int j=a; j<=b; j++)
		{
			ll x=f(rzo[i], kol[j]);
			if (x>top)
			{
				top=x;
				g=j;
			}
		}
		ret[i]=g;
	}
	for (int &i : ret)
		i=stos[i];
	return ret;
}
