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
		if (min(a.size(), b.size())<=64)
		{
			vll ret(a.size()+b.size()-1);
			for (int i=0; i<(int)a.size(); i++)
				for (int j=0; j<(int)b.size(); j++)
					ret[i+j]=(ret[i+j]+a[i]*b[j])%mod;
			return ret;
		}
		int n=lift((int)a.size()+(int)b.size()-1);
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

vll integrate(const vll &a)
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
	return modulo(q*(vll{2}-a*q), a.size());
}

vll logarithm(const vll &a)//a[0]!=0
{
	vll deri=derivative(a);
	vll inve=inverse(a);
	return modulo(integrate(deri*inve), a.size());
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
	return modulo(q*(vll{1}+a-logarithm(qp)), a.size());
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

vector<vll> transpose(const vector<vll> &a)
{
	if (a.empty())
		return {};
	int ra=a.size();
	int rb=a[0].size();
	vector<vll> ret(rb, vll(ra));
	for (int i=0; i<ra; i++)
		for (int j=0; j<rb; j++)
			ret[j][i]=a[i][j];
	return ret;
}

vector<vll> operator *(const vector<vll> &a, const vector<vll> &b)
{
	if (a.empty() || b.empty())
		return {};
	int oa=(int)a.size()+(int)b.size()-1;
	int ob=(int)a[0].size()+(int)b[0].size()-1;
	vll pa(oa*ob);
	vll pb(oa*ob);
	for (int i=0; i<(int)a.size(); i++)
		for (int j=0; j<(int)a[i].size(); j++)
			pa[i*ob+j]=a[i][j];
	for (int i=0; i<(int)b.size(); i++)
		for (int j=0; j<(int)b[i].size(); j++)
			pb[i*ob+j]=b[i][j];
	vll pc=pa*pb;
	pc.resize(oa*ob);
	vector<vll> ret(oa, vll(ob));
	for (int i=0; i<oa; i++)
		for (int j=0; j<ob; j++)
			ret[i][j]=pc[i*ob+j];
	return ret;
}

vll power_projection(const vll &a, const vll &b, int n)
{
	if (a.empty() || b.empty() || !n)
		return vll(n);
	vector<vll> ta{a};
	vector<vll> tb{{1}, b};
	reverse(ta[0].begin(), ta[0].end());
	tb[0].resize(a.size());
	tb[1].resize(a.size());
	for (int i=0; i<(int)a.size(); i++)
		tb[1][i]=(mod-tb[1][i])%mod;
	while((int)ta[0].size()>1)
	{
		int r=ta[0].size();
		auto rb=tb;
		for (int i=0; i<(int)rb.size(); i++)
			for (int j=1; j<(int)rb[0].size(); j+=2)
				rb[i][j]=(mod-rb[i][j])%mod;
		auto pa=ta*rb;
		auto pb=tb*rb;
		ta.resize(pa.size());
		tb.resize(pb.size());
		for (int i=0; i<(int)ta.size(); i++)
		{
			ta[i].resize((r+1)/2);
			for (int j=0; j<(r+1)/2; j++)
				ta[i][j]=pa[i][2*j+1-(r&1)];
		}
		for (int i=0; i<(int)tb.size(); i++)
		{
			tb[i].resize((r+1)/2);
			for (int j=0; j<(r+1)/2; j++)
				tb[i][j]=pb[i][2*j];
		}
	}
	vll pa=transpose(ta)[0];
	vll pb=transpose(tb)[0];
	pa.resize(n);
	pb.resize(n);
	vll ret=pa*inverse(pb);
	ret.resize(n);
	return ret;
}

vll composition(const vll &a, const vll &b)
{
	if (a.empty())
		return {};
	function<vector<vll>(vector<vll>)> reverse_2d=[&](vector<vll> a)
	{
		reverse(a.begin(), a.end());
		for (int i=0; i<(int)a.size(); i++)
			reverse(a[i].begin(), a[i].end());
		return a;
	};
	function<vector<vll>(vector<vll>, vector<vll>)> multiply_by_reversed=[&](const vector<vll> &a, const vector<vll> &b)
	{
		if (a.empty() || b.empty())
			return a;
		auto d=a*reverse_2d(b);
		auto ret=a;
		for (int i=0; i<(int)a.size(); i++)
			for (int j=0; j<(int)a[0].size(); j++)
				ret[i][j]=d[i+(int)b.size()-1][j+(int)b[0].size()-1];
		return ret;
	};
	vector<vll> tb{{1}, b};
	tb[0].resize(a.size());
	tb[1].resize(a.size());
	for (int i=0; i<(int)a.size(); i++)
		tb[1][i]=(mod-tb[1][i])%mod;
	vector<vector<vll>> rem;
	vi bits;
	while((int)tb[0].size()>1)
	{
		int r=tb[0].size();
		auto rb=tb;
		for (int i=0; i<(int)rb.size(); i++)
			for (int j=1; j<(int)rb[0].size(); j+=2)
				rb[i][j]=(mod-rb[i][j])%mod;
		auto pb=tb*rb;
		tb.resize(pb.size());
		rem.push_back(rb);
		bits.push_back(1-(r&1));
		for (int i=0; i<(int)tb.size(); i++)
		{
			tb[i].resize((r+1)/2);
			for (int j=0; j<(r+1)/2; j++)
				tb[i][j]=pb[i][2*j];
		}
	}
	vll pb=transpose(tb)[0];
	pb=inverse(pb);
	pb.resize(a.size());
	tb=transpose({pb});
	auto ta=transpose({a});
	ta=multiply_by_reversed(ta, tb);
	int lowering=0;
	for (auto &i : rem)
		lowering+=i.size()-1;
	for (int h=(int)bits.size()-1; h>=0; h--)
	{
		auto &mul=rem[h];
		int bit=bits[h];
		vector<vll> pa(ta.size(), vll((int)ta[0].size()*2-1+bit));
		for (int i=0; i<(int)ta.size(); i++)
			for (int j=0; j<(int)ta[i].size(); j++)
				pa[i][2*j+bit]=ta[i][j];
		ta=multiply_by_reversed(pa, mul);
		lowering-=mul.size()-1;
		ta.resize(lowering+1);
	}
	ta[0].resize(a.size());
	reverse(ta[0].begin(), ta[0].end());
	return norm(ta[0]);
}

vll sums_of_root_powers(const vll &a, int n)
{
	if (a.empty() || !n)
		return vll(n, 0);
	vll b=a;
	reverse(b.begin(), b.end());
	if (n>(int)a.size())
		b.resize(n-1);
	vll ret=vll{0}-derivative(b)*inverse(b);
	ret.resize(n-1);
	ret.insert(ret.begin(), (int)a.size()-1);
	return ret;
}

vll compositional_inverse(vll a)
{
	if ((int)a.size()<=1)
		return vll(a.size());
	int n=a.size();
	ll c=inv(a[1]);
	for (int i=0; i<n; i++)
		a[i]=a[i]*c%mod;
	vll vals(n);
	vals[n-1]=1;
	vll b=power_projection(vals, a, n);
	for (int i=1; i<n; i++)
		b[i]=b[i]*(n-1)%mod*inv(i)%mod;
	reverse(b.begin(), b.end());
	b=logarithm(b);
	b.resize(n);
	for (int i=0; i<n; i++)
		b[i]=b[i]*(mod-inv(n-1))%mod;
	b=exponent(b);
	b.insert(b.begin(), 0);
	ll w=1;
	for (ll &i : b)
	{
		i=i*w%mod;
		w=w*c%mod;
	}
	return modulo(b, n);
}

vector<vector<vll>> half_gcd(vll a, vll b)
{
	vll oa=a;
	vll ob=b;
	vector<vll> mula{{1}, {}};
	vector<vll> mulb{{}, {1}};
	int goal=a.size()/2;
	int max_step=((int)a.size()-goal+1)/2;
	while(1)
	{
		if ((int)a.size()<=goal)
			return {mulb, mula};
		{
			auto take=divide(b, a);
			b=take.second;
			mulb[0]=mulb[0]-mula[0]*take.first;
			mulb[1]=mulb[1]-mula[1]*take.first;
		}
		if ((int)b.size()<=goal)
			return {mula, mulb};
		int r=a.size();
		int goal_now=max(goal, r-max_step);
		int sta=max(2*goal_now-r, 0);
		vll na, nb;
		for (int i=sta; i<(int)a.size(); i++)
			na.push_back(a[i]);
		for (int i=sta; i<(int)b.size(); i++)
			nb.push_back(b[i]);
		{
			auto take=half_gcd(na, nb);
			na=take[0][0]*a+take[0][1]*b;
			nb=take[1][0]*a+take[1][1]*b;
			a=nb;
			b=na;
			vector<vll> nmula={take[0][0]*mula[0]+take[0][1]*mulb[0], take[0][0]*mula[1]+take[0][1]*mulb[1]};
			vector<vll> nmulb={take[1][0]*mula[0]+take[1][1]*mulb[0], take[1][0]*mula[1]+take[1][1]*mulb[1]};
			mula=nmulb;
			mulb=nmula;
		}
	}
};

pair<vll,pair<vll,vll>> gcd(const vll &a, const vll &b)//first=a*second.first+b*second.second
{
	if (a.empty())
	{
		if (b.empty())
			return {{}, {{}, {}}};
		ll x=inv(b.back());
		return {norm(b*vll{x}), {{}, norm({x})}};
	}
	auto take1=half_gcd(a, b);
	vll na=a*take1[0][0]+b*take1[0][1];
	vll nb=a*take1[1][0]+b*take1[1][1];
	auto take2=gcd(nb, na);
	return {take2.first, {take1[0][0]*take2.second.second+take1[1][0]*take2.second.first, take1[0][1]*take2.second.second+take1[1][1]*take2.second.first}};
}

vll power_modulo(vll a, ll expo, vll m)
{
	vll ret=divide({1}, m).second;
	while(expo)
	{
		if (expo&1)
			ret=divide(ret*a, m).second;
		expo>>=1;
		if (expo)
			a=divide(a*a, m).second;
	}
	return ret;
}

vll root_finding(const vll &a, bool repeated)
{
	vll ret;
	function<void(vll)> search=[&](vll b)
	{
		if ((int)b.size()<=1)
			return;
		if ((int)b.size()==2)
		{
			ret.push_back((mod-b[0])*inv(b[1])%mod);
			return;
		}
		vll c;
		for (int i=1; i<(int)b.size(); i++)
			c.push_back(rand()%mod);
		c=power_modulo(c, (mod-1)/2, b)-vll{1};
		if (c.empty())
		{
			search(b);
			return;
		}
		c=gcd(c, b).first;
		search(c);
		search(divide(b, c).first);
	};
	function<void(vll, vll, int)> extract=[&](vll b, vll prod, int expo)
	{
		prod=divide(prod, b).second;
		vector<vll> powers{prod};
		for (int i=1; i<expo; i++)
		{
			vll c=divide(powers.back()*powers.back(), b).second;
			powers.push_back(c);
		}
		for (int i=expo; i>=0; i--)
		{
			vll g=gcd(b, powers[max(0, i-1)]).first;
			if (!i)
			{
				search(g);
			}
			else
			{
				extract(divide(b, g).first, prod, i-1);
				b=g;
			}
		}
	};
	int expo=0;
	while((1<<expo)<(int)a.size()-1)
		expo++;
	extract(a, power_modulo({0, 1}, mod, a)-vll{0, 1}, repeated ? expo : 0);
	sort(ret.begin(), ret.end());
	return ret;
}

vll product(const vector<vll> &seq)
{
	function<vll(int, int)> dc=[&](int a, int b)
	{
		if (a>b)
			return vll{1};
		if (a==b)
			return seq[a];
		int s=(a+b)>>1;
		return dc(a, s)*dc(s+1, b);
	};
	return dc(0, (int)seq.size()-1);
}
