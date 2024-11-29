//zaczyna z x = 0 i ans
//n razy wykonuje x = (x + a) % m
//jeśli x się nie przekręci:
//  jeśli x < c, to do ans przykłada AX
//  jeśli x >= c, to do ans przykłada AY
//jeśli x się przekręci:
//  jeśli x < c, to do ans przykłada BX
//  jeśli x >= c, to do ans przykłada BY
//w powyższych porównaniach patrzymy na wartość x przed dodaniem do niej a
//struct ans musi mieć przeciążony operator mnożenia (składania monoidu)

template<class F>
int find_max_true(int lw,int up,F f){
	while(up-lw>1){
		const int mid=(lw+up)/2;
		if(f(mid))lw=mid;
		else up=mid;
	}
	return lw;
}

pll abnum(ll a,ll b,ll n){
	assert(1<=a);
	assert(1<=b);
	ll x=find_max_true(0,n+1,[&](ll z){
		return z+(a*z)/b<=n;
	});
	ll y=(a*x)/b;
	ll cur=x*a-y*b;
	n-=x+y;
	if(n){
		cur+=a;
		x++;
		n--;
		assert(cur/b>=n);
		y+=n;
		n=0;
	}
	return pll(x,y);
}

ll fld(ll a, ll b) { // floored division
	return a / b - ((a ^ b) < 0 && a % b); }
ll cld(ll a, ll b) { // ceiled division
	return a / b + ((a ^ b) > 0 && a % b); }

template<class N>
N gauss_sum_monoid_super(ll a,ll m,ll c,ll n,N ax,N ay,N bx,N by,N ans){
	ll b=m-a;
	assert(1<=a);
	assert(1<=b);
	assert(gcd(a,b)==1);
	assert(0<=c && c<=a+b);
	assert(0<=n);
	
	auto x_pow=[&](N&x,N v,ll k){
		while(k){
			if(k&1)x=x*v;
			if(k>>=1)v=v*v;
		}
	};
	auto pow_x=[&](N v,ll k,N&x){
		while(k){
			if(k&1)x=v*x;
			if(k>>=1)v=v*v;
		}
	};
	
	auto [an,bn]=abnum(a,b,n);
	while(an||bn){
		assert(an>=0);
		assert(bn>=0);
		assert(0<=c && c<=a+b);
		if(a==0){
			assert(b==1);
			assert(bn==0);
			pow_x(0<c?ax:ay,an,ans);
			break;
		}
		if(b==0){
			assert(a==1);
			assert(an==0);
			pow_x(0<c?bx:by,bn,ans);
			break;
		}
		assert(a>0);
		assert(b>0);
		if(a<b){
			{
				ll fin=(a*an-b*bn);
				ll p=fin/a,q=fin%a;
				ll u=min(p,cld(c-q,a));
				pow_x(ay,p-u,ans);
				pow_x(ax,u,ans);
				an-=p;
				fin-=a*p;
				assert(an>=0);
			}
			ll p=b/a,q=b%a;
			if(c<=q){
				pow_x(ay,p,bx);
				pow_x(ay,p,by);
			}else{
				auto sub=[&](ll z){
					ll u=min(p,cld(c-z,a));
					N res=z+p*a<c?bx:by;
					pow_x(ay,p-u,res);
					pow_x(ax,u,res);
					return res;
				};
				ll cnx=q+(c-q)%a;
				N bxnx=sub(q);
				N bynx=sub(cnx);
				bx=bxnx;
				by=bynx;
				c=cnx;
			}
			b-=a*p;
			an-=bn*p;
		}else{
			{
				ll fin=(a*an-b*bn);
				ll p=(a+b-1-fin)/b;
				ll u=min(p,cld(fin+b*p+1-c,b));
				pow_x(bx,p-u,ans);
				pow_x(by,u,ans);
				bn-=p;
				fin+=b*p;
				assert(bn>=0);
				ans=(fin-a<c?ax:ay)*ans;
				an--;
				fin-=a;
				assert(an>=0);
			}
			ll p=a/b,q=a%b;
			if(c<=b+q){
				x_pow(ax,by,p);
				x_pow(ay,by,p);
			}else{
				auto sub=[&](ll z){
					N res=z<c?ax:ay;
					ll u=min(p,cld(z+a+1-c,b));
					x_pow(res,by,u);
					x_pow(res,bx,p-u);
					return res;
				};
				ll cnx=(c-q)%b;
				N axnx=sub(0);
				N aynx=sub(cnx);
				ax=axnx;
				ay=aynx;
				by=bx;
				c=cnx;
			}
			a-=b*p;
			bn-=an*p;
		}
	}
	return ans;
}
