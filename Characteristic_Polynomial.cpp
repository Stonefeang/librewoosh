const ll mod=998244353;

ll inv(ll v)
{
	if (v<=1)
		return v;
	return inv(mod%v)*(mod-mod/v)%mod;
}

vll characteristic(vector<vll> v)
{
	int r=v.size();
	auto zam=[&](int a, int b)
	{
		swap(v[a], v[b]);
		for (int i=0; i<r; i++)
			swap(v[i][a], v[i][b]);
	};
	auto odj=[&](int a, int b, ll x)
	{
		for (int i=0; i<r; i++)
			v[a][i]=(v[a][i]+v[b][i]*(mod-x))%mod;
		for (int i=0; i<r; i++)
			v[i][b]=(v[i][b]+v[i][a]*x)%mod;
	};
	for (int i=0; i+1<r; i++)
	{
		for (int j=i+1; j<r; j++)
		{
			if (v[j][i])
			{
				zam(i+1, j);
				break;
			}
		}
		ll dz=inv(v[i+1][i]);
		for (int j=i+2; j<r; j++)
			odj(j, i+1, v[j][i]*dz%mod);
	}
	for (int i=0; i<r; i++)
		for (int j=0; j<r; j++)
			v[i][j]=(mod-v[i][j])%mod;
	vector<vll> poly(r+1, vll(r+1));
	poly[0][0]=1;
	for (int i=0; i<r; i++)
	{
		for (int j=0; j<=i; j++)
		{
			poly[i+1][j]=(poly[i+1][j]+poly[i][j]*v[i][i])%mod;
			poly[i+1][j+1]=(poly[i+1][j+1]+poly[i][j])%mod;
		}
		ll x=1;
		for (int j=i; j+1<r; j++)
		{
			x=x*v[j+1][j]%mod;
			ll y=x*v[i][j+1]%mod;
			if ((j&1)==(i&1))
				y=(mod-y)%mod;
			for (int l=0; l<=i; l++)
				poly[j+2][l]=(poly[j+2][l]+poly[i][l]*y)%mod;
		}
	}
	return poly[r];
}
