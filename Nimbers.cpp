using ull=unsigned long long;
const int d=64;
 
ull rem[d][d];
 
ull mno(ull a, ull b)
{
	ull ret=0;
	for (int i=0; i<d; i++)
		for (int j=0; j<d; j++)
			if ((a&(1ULL<<i)) && (b&(1ULL<<j)))
				ret^=rem[i][j];
	return ret;
}
 
void prepro()
{
	for (int i=0; i<d; i++)
	{
		for (int j=0; j<=i; j++)
		{
			if (!j)
			{
				rem[j][i]=rem[i][j]=1ULL<<i;
				continue;
			}
			int p=0;
			while((1<<(p+1))<=i)
				p++;
			if (j&(1<<p))
				rem[j][i]=rem[i][j]=rem[i^(1<<p)][j]^mno(rem[i^(1<<p)][j^(1<<p)], 1ULL<<((1<<p)-1));
			else
				rem[j][i]=rem[i][j]=rem[i^(1<<p)][j]<<(1<<p);
		}
	}
}
