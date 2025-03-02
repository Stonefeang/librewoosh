const ll mod=998244353;

int invfast[nax];

void init_inv()//wywolac na poczatku
{
	invfast[1]=1;
	for (int i=2; i<nax; i++)
		invfast[i]=invfast[mod%i]*(mod-mod/i)%mod;
}

ll inv(ll v)
{
	if (v<nax)
		return invfast[v];
	return inv(mod%v)*(mod-mod/v)%mod;
}

struct ntt
{
	vll omega[30];
	ll gen=1;
	ntt()
	{
		while(1)
		{
			gen++;
			ll x=gen;
			for (int i=1; i<__builtin_ctzll(mod-1); i++)
				x=(x*x)%mod;
			if (x==mod-1)
				break;
		}
	}
	int lift(int v)
	{
		int ret=1;
		while(ret<v)
			ret<<=1;
		return ret;
	}
	void omegas(int v)
	{
		if (!omega[v].empty())
			return;
		int n=(1<<v);
		int should=((mod-1)&(-(mod-1)));
		ll mul=gen;
		while(n<should)
		{
			mul=(mul*mul)%mod;
			should>>=1;
		}
		omega[v].resize(n+1);
		omega[v][0]=1;
		for (int i=1; i<=n; i++)
			omega[v][i]=(omega[v][i-1]*mul)%mod;
	}
	void dft(vll &a, int dir)
	{
		int n=a.size();
		static vll b;
		b.resize(n);
		const int ch=(!dir ? 1 : -1);
		for (int i=1, w=1; i<n; i<<=1, w++)
		{
			omegas(w);
			ll *om=omega[w].data();
			b.swap(a);
			const int &d=n>>w;
			ll *pa=a.data();
			ll *pb=b.data();
			int now=(!dir ? 0 : (i<<1));
			for (int j=0; j<n; j+=d, now+=ch)
			{
				const ll &mul=om[now];
				int left=(j<<1);
				if (left>=n)
					left-=n;
				int right=(left+d);
				for (int l=0; l<d; l++)
					pa[j+l]=(pb[left+l]+pb[right+l]*mul)%mod;
			}
		}
	}
	vll multi(vll a, vll b)
	{
		if (a.empty() || b.empty())
			return {};
		int n=lift(a.size()+b.size());
		a.resize(n);
		b.resize(n);
		dft(a, 0);
		dft(b, 0);
		for (int i=0; i<n; i++)
			a[i]=(a[i]*b[i])%mod;
		dft(a, 1);
		ll div=inv(n);
		for (ll &i : a)
			i=(i*div)%mod;
		return a;
	}
};

vll norm(vll a)
{
	for (ll &i : a)
		if ((i%=mod)<0)
			i+=mod;
	while(!a.empty() && !a.back())
		a.pop_back();
	return a;
}

vll operator + (const vll &a, const vll &b)
{
	vll ret(max(a.size(), b.size()));
	for (int i=0; i<(int)a.size(); i++)
		ret[i]=a[i];
	for (int i=0; i<(int)b.size(); i++)
		ret[i]+=b[i];
	return norm(ret);
}

vll operator - (const vll &a, const vll &b)
{
	vll ret(max(a.size(), b.size()));
	for (int i=0; i<(int)a.size(); i++)
		ret[i]=a[i];
	for (int i=0; i<(int)b.size(); i++)
		ret[i]-=b[i];
	return norm(ret);
}

vll operator *(const vll &a, const vll &b)
{
	static ntt my_ntt;
	return norm(my_ntt.multi(norm(a), norm(b)));
}

vll modulo(vll a, int m)
{
	a.resize(m);
	return norm(a);
}

vll derivative(const vll &a)
{
	if (a.empty())
		return {};
	vll ret((int)a.size()-1);
	for (int i=0; i<(int)ret.size(); i++)
		ret[i]=a[i+1]*(i+1);
	return norm(ret);
}

vll integral(const vll &a)
{
	vll ret((int)a.size()+1);
	for (int i=1; i<(int)ret.size(); i++)
		ret[i]=a[i-1]*invfast[i];
	return norm(ret);
}

vll inverse(const vll &a)//a[0]!=0
{
	int r=a.size();
	if (r<=1)
		return {inv(a[0])};
	while(r&(r-1))
		r++;
	vll p=a;
	p.resize(r>>1);
	vll q=inverse(p);
	return modulo(q*(vll{2}-a*q), r);
}

vll logarithm(const vll &a)//a[0]!=0
{
	vll deri=derivative(a);
	vll inve=inverse(a);
	return modulo(integral(deri*inve), a.size());
}

vll exponent(const vll &a)//a[0]=0
{
	int r=a.size();
	if (r<=1)
		return {1};
	while(r&(r-1))
		r++;
	vll p=a;
	p.resize(r>>1);
	vll q=exponent(p);
	vll qp=q;
	qp.resize(r);
	return modulo(q*(vll{1}+a-logarithm(qp)), r);
}

pair<vll,vll> divide(const vll a, const vll &b)//{div, mod}
{
	if ((int)a.size()<(int)b.size())
		return {{}, a};
	int r=(int)a.size()-(int)b.size()+1;
	vll oa=a;
	vll ob=b;
	reverse(oa.begin(), oa.end());
	reverse(ob.begin(), ob.end());
	oa.resize(r);
	ob.resize(r);
	vll ret=oa*inverse(ob);
	ret.resize(r);
	reverse(ret.begin(), ret.end());
	ret=modulo(ret, r);
	return {ret, a-b*ret};
}

vll multipoint_evaluation(const vll &w, const vll &points)
{
	if (points.empty())
		return {};
	vector <vll> help;
	vll res;
	function<void(int, int, int)> precalc=[&](int v, int a, int b)
	{
		if ((int)help.size()<=v) help.resize(v+1);
		if (a==b)
		{
			help[v]={(mod-points[a])%mod, 1};
			return;
		}
		precalc((v<<1), a, (a+b)>>1);
		precalc((v<<1)^1, (a+b+2)>>1, b);
		help[v]=(help[v*2]*help[v*2+1]);
	};
	function<void(int, vll)> calc=[&](int v, vll wek)
	{
		wek=divide(wek, help[v]).second;
		if ((int)help[v].size()==2)
		{
			res.push_back(wek.empty() ? 0 : wek[0]);
			return;
		}
		calc((v<<1), wek);
		calc((v<<1)^1, wek);
	};
	precalc(1, 0, (int)points.size()-1);
	calc(1, norm(w));
	return res;
}

vll interpolate(const vector<pll> &points)
{
	if (points.empty())
		return {};
	vector <vll> help;
	function<void(int, int, int)> precalc=[&](int v, int a, int b)
	{
		if ((int)help.size()<=v) help.resize(v+1);
		if (a==b)
		{
			help[v]={(mod-points[a].first)%mod, 1};
			return;
		}
		precalc((v<<1), a, (a+b)>>1);
		precalc((v<<1)^1, (a+b+2)>>1, b);
		help[v]=(help[v*2]*help[v*2+1]);
	};
	function<vll(int, int, int, vll)> calc=[&](int v, int a, int b, vll wek)
	{
		wek=divide(wek, help[v]).second;
		if (a==b)
			return norm({points[a].second*inv(wek[0])});
		vll l=calc((v<<1), a, (a+b)>>1, wek*help[(v<<1)^1]);
		vll r=calc((v<<1)^1, (a+b+2)>>1, b, wek*help[(v<<1)]);
		return (l*help[(v<<1)^1])+(r*help[(v<<1)]);
	};
	precalc(1, 0, (int)points.size()-1);
	return calc(1, 0, (int)points.size()-1, {1});
}
