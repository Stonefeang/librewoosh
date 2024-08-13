ll floor_sum(ll k, ll a, ll b, ll m)//sum of (a*i+b)/m for 0<=i<k
{
	if (!k)
		return 0;
	ll ret=k*(b/m)+(k*(k-1)/2)*(a/m);
	b%=m;
	a%=m;
	return ret+floor_sum((k*a+b)/m, m, (k*a+b)%m, a);
}
