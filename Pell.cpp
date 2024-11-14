pll pell(ll n) {
	ll s=llround(sqrtl(n))+2;
	while(s*s>n)
		s--;
	if (s*s==n)
		return {0, 0};
	ll m=0, d=1, a=s;
	__int128 num1=1, num2=a, den1=0, den2=1;
	while (num2*num2-n*den2*den2!=1)
	{
		m=d*a-m;
		d=(n-m*m)/d;
		a=(s+m)/d;
		if (num2>(1LL<<62)/a)
			return {0, 0};
		tie(num1, num2)=pll(num2, a*num2+num1);
		tie(den1, den2)=pll(den2, a*den2+den1);
	}
	return {num2, den2};
}
vector<pll> all_pell(ll n, ll limit) {
	auto [x0, y0] = pell(n);
	if (!x0)
		return {};
	vector<pll> ret;
	__int128 x=x0, y=y0;
	while (x<=limit)
	{
		ret.push_back({x, y});
		if (y0*y>(1ll<<62)/n)
			break;
		tie(x, y)=pll(x0*x+n*y0*y, x0*y+y0*x);
	}
	return ret;
}
