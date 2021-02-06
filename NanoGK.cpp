#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<iomanip>

#define filename "Data-160"
#define N0 7.5
#define K a 
#define D b
#define ode (K*D/(K*y+D))*(1.0*0.001028-0.001028*0.6778*N0*y*y*y)	//define your ode as y'=f(x,y)
#define func RK4
#define m 7										//define the starting point
#define n 10									//define the half-maximum point
#define l 23 									//define the end point
#define Tmax 100							//define the starting temp.
#define Tmin 0.1								//define the ending temp.
#define rate_cool 0.99							//define the cooling coefficient
#define R 0.01									//define the Eng/RT
#define Engf LSR(a,b,c,expm)

double LSR(double a, double b, double c, double expm[l][2]){
	double x,y,SSE=0.0;
	double thr[l][2]={0};
	for(int i=0;i<l;i++){
		thr[i][0]=expm[i][0];
		}
	double xx=expm[n][0],yy=pow(expm[n][1]/(4.1888*N0),1.0/3.0),k1,k2,k3,k4,h=0.0001;
	thr[n][0]=expm[n][0], thr[n][1]=expm[n][1];
	int k=n+1;
	while(xx<40){
			y=yy;x=xx;
	        k1=ode;
	        y=yy+0.5*h*k1;x=xx+0.5*h;
	        k2=ode;
	        y=yy+0.5*h*k2;x=xx+0.5*h;
	        k3=ode;
	        y=yy+h*k3;x=xx+h;
	        k4=ode;
	        yy=yy+h*(k1+2*k2+2*k3+k4)/6.0;xx=xx+h;
	        if(fabs(xx-expm[k][0])<0.000000001){
				thr[k][1]=4.1888*N0*pow(yy,3.0);
				//printf("%lf\t%lf\n",thr[k][0],thr[k][1]);
				k=k+1;
				}
			else;
			}					//RK4-forward
	xx=expm[n][0], yy=pow(expm[n][1]/(4.1888*N0),1.0/3.0), k=n-1;
	while(xx>=expm[m][0]&&yy>0){
			y=yy;x=xx;
	        k1=ode;
	        y=yy-0.5*h*k1;x=xx-0.5*h;
	        k2=ode;
	        y=yy-0.5*h*k2;x=xx-0.5*h;
	        k3=ode;
	        y=yy-h*k3;x=xx-h;
	        k4=ode;
	        yy=yy-h*(k1+2*k2+2*k3+k4)/6.0;xx=xx-h;
	        if(fabs(xx-expm[k][0])<0.000000001){
	        	thr[k][0]=xx;
				thr[k][1]=4.1888*N0*pow(yy,3.0);
				//printf("%lf\t%lf\n",thr[k][0],thr[k][1]);
				k=k-1;
				}
			else;
			}					//RK4-backward
	for (int i=m;i<l;i++)
		{
		//printf("%lf\t%lf\n",thr[i][0],thr[i][1]);
		double slope=thr[i][1]-thr[i-1][1];
		SSE=SSE+pow((thr[i][1]-expm[i][1]),2)*(1/(slope*slope+1));
		}

return (SSE);
}								//calc. of 1-r^2, "The Coefficient of Determination"

double pLSR(double a, double b, double c, double expm[l][2]){
	double x,y,SSE=0.0,EX2=0.0,EX=0.0,SST=0.0;
	double thr[l][2]={0};
	for(int i=0;i<l;i++){
		thr[i][0]=expm[i][0];
		}
	double xx=expm[n][0],yy=pow(expm[n][1]/(4.1888*N0),1.0/3.0),k1,k2,k3,k4,h=0.0001;
	thr[n][0]=expm[n][0], thr[n][1]=expm[n][1];
	int k=n+1;
	while(xx<expm[l-1][0]){
			y=yy;x=xx;
	        k1=ode;
	        y=yy+0.5*h*k1;x=xx+0.5*h;
	        k2=ode;
	        y=yy+0.5*h*k2;x=xx+0.5*h;
	        k3=ode;
	        y=yy+h*k3;x=xx+h;
	        k4=ode;
	        yy=yy+h*(k1+2*k2+2*k3+k4)/6.0;xx=xx+h;
	        if(fabs(xx-expm[k][0])<0.000001){
				thr[k][1]=4.1888*N0*pow(yy,3.0);
				//printf("%lf\t%lf\n",thr[k][0],thr[k][1]);
				k=k+1;
				}
			else;
			}					//RK4-forward
	xx=expm[n][0], yy=pow(expm[n][1]/(4.1888*N0),1.0/3.0), k=n-1;
	while(xx>=expm[m][0]&&yy>0){
			y=yy;x=xx;
	        k1=ode;
	        y=yy-0.5*h*k1;x=xx-0.5*h;
	        k2=ode;
	        y=yy-0.5*h*k2;x=xx-0.5*h;
	        k3=ode;
	        y=yy-h*k3;x=xx-h;
	        k4=ode;
	        yy=yy-h*(k1+2*k2+2*k3+k4)/6.0;xx=xx-h;
	        if(fabs(xx-expm[k][0])<0.000001){
	        	thr[k][0]=xx;
				thr[k][1]=4.1888*N0*pow(yy,3.0);
				//printf("%lf\t%lf\n",thr[k][0],thr[k][1]);
				k=k-1;
				}
			else;
			}					//RK4-backward
	
	FILE *fp=NULL;
	char fname[64];
	strcpy (fname,filename);
	strcat (fname,"-Check-FB.txt");
	fp=fopen(fname,"w");   //name your file
	for (int i=0;i<l;i++)
		{
		printf("%lf\t%lf\n",thr[i][0],thr[i][1]);
		fprintf(fp,"%lf\t%lf\n",thr[i][0],thr[i][1]);
		}
return (0);
}	

int main()
{
char fname[64];
strcpy (fname,filename);
strcat (fname,".txt");
double expm[l][2]={0};  
double a=1,b=1,c=1,aa,bb,cc,aaa,bbb,ccc,T=Tmax,S=0.0,r=rate_cool,Eng,Eng00,Eng0;
Eng00=65535; Eng0=Eng00; 
srand((unsigned)time(0));  

double data[l][2]={0};  
    FILE *fpRead=fopen(fname,"r");  						//input data file
    for(int i=0;i<l;i++){  
        fscanf(fpRead,"%lf\t%lf\n",&data[i][0],&data[i][1]);   
    }
	fclose(fpRead);     									//Read expm. values
	for(int i=0;i<l;i++){
	expm[i][0]=data[i][0];
	expm[i][1]=data[i][1];
	printf("%lf\t%lf\n",expm[i][0],expm[i][1]);
	}
	
FILE *fp=NULL;
strcpy (fname,filename);
strcat (fname,"-Result-FB.txt");
fp=fopen(fname,"w");   										//name the file
fprintf(fp,"a\t b\t R2\n");
printf("T\tK\tD\tSSR\n");
//////SimulateD Annealing///////

for(int j=0;j<200;j++){
T=Tmax; a=1; aa=0; b=1; bb=0; c=1; cc=1; ccc=0; Eng0=Engf;
while(T>Tmin)
	{
	a=fabs((T-2*T*rand()/double(RAND_MAX))+aa);
	b=fabs((T-2*T*rand()/double(RAND_MAX))+bb);
	c=fabs((T-2*T*rand()/double(RAND_MAX))+cc);
	Eng=Engf;
	//printf("%d\t%.2f\t%f\t%f\t%f\t%f\n",j,T,a,b,c,Eng);
	if((Eng<=Eng0)){
		T=T*r;
		Eng0=Eng;
		aa=a;
		bb=b;
		cc=c;
		printf("%d\t%.2f\t%f\t%f\t%f\t%f\n",j,T,aa,bb,cc,Eng0);
		}
	else{if((exp(-(Eng-Eng0)/(R*T))>rand()/double(RAND_MAX))){
			Eng0=Eng;
			aa=a;
			bb=b;
			cc=c;
			T=T*r;
			}
		}
	}
	printf("%d\t%.2f\t%f\t%f\t%f\t%f\n",j,T,aa,bb,cc,Eng0);
	fprintf(fp,"%f\t%f\t%f\t%.4f\n",aa,bb,cc,Eng0);
	if(Eng0<Eng00){
		Eng00=Eng;
		aaa=aa;
		bbb=bb;
		ccc=cc;
		}
	else;
}
printf("Best try:\t%f\t%f\t%f\t%f\n",aaa,bbb,ccc,Eng00);
pLSR(aaa,bbb,ccc,expm);

fclose(fp);

return (0);

}
