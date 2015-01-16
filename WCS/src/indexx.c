void indexx(n,arrin,indx)
int n,indx[];
float arrin[];
{
	int l,j,ir,indxt,i;
	float *temp,q;
	void free_vector();
	float *vector();

	temp = vector(1,n);
	

	for (j=1;j<=n;j++) indx[j]=j;
	l=(n >> 1) + 1;
	ir=n;
	for (;;) {
		if (l > 1)
			q=arrin[(indxt=indx[--l])];
		else {
			q=arrin[(indxt=indx[ir])];
			indx[ir]=indx[1];
			if (--ir == 1) {
				indx[1]=indxt;
				for (j=1;j<=n;j++) temp[j]=arrin[j];
				for (j=1;j<=n;j++) arrin[j] = temp[indx[j]];
				free_vector(temp,1,n);
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
			if (q < arrin[indx[j]]) {
				indx[i]=indx[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		indx[i]=indxt;
	}
}


