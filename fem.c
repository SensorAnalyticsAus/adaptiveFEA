#define T800         0
#define Less         1 /*ok[] is not used*/
#define Timer        0 
#define DEBUG        0 /* prints ok[] and solution vectors */
#define DEBUGT800    0 /*prints the reflected data */
#define Maxproc      20
#define RHSsub       100   /* 50  nodes/proc and 50 band for debug */
#define MaxOk        RHSsub*50 /* part matrix */
#define CLK          15625

/*This a FEM program for constant strain triangular elements with
  an adaptivity module */


  /* Master main() */

#include <stdio.h>

#include <stdlib.h>

#include <string.h>

#include <math.h>

#include <time.h>

#include <ctype.h>

#if T800

#include <timer.h>

#include <chan.h>

#include <alt.h>

#if Timer
#include <timer.h>
#endif

#endif

#define NE 200

#define MM 25

#define LM  5

#define BUFF 256

FILE *fp2,*fp3;

char buf[BUFF];



int nn,ne,nb,nm,nl,ndime,nstr,nbe,kmax,nmax,iln,nq;

int *nod, nco[MM][4], mat[5];

int nsb;

double dt,damp;

double *cord;

double *ste;

double amat[5][10];

double est[6][6], *s;

int nn2, arg_def;

float *ok;

double *f;

int iter, opt_proc;

unsigned gauss;

time_t t1, t2;

struct {int nc; double aload[3];} NCLOAD[LM];

void fst(int n);
void /*FUNCTION*/ setup(int nn, int ne, double *cord, int *nod);
void routput(char *file_sav);
void print_orig_data(void);
void rinput(char *file_sav);
void openfile_fp2(char *file_sav);
void openfile_fp3(char *file_sav);
void semi_band(void);
void elem(void);
void zero_ok(void);
void aload_f(void);
void form_ok(void);
void bc_ok(void);
void gauss_el(void);
void fpr_defl(void);
void calc_stress(void);
void form_D(double v, double e, double d[][3]);
char* trimwhitespace(char *str);
void chkf(char *str);

#if T800
void preliminary_data(void);
void distribute_rows_per_processor(void);
void ok_transmitter(int nrow, int ncol, int ja, int jb, int flag);
void load_vector_distributer(void);
void receive_and_print_reflected_data(void);
void receive_iter(void);
void  receive_data_f(void);


#if Timer
void rec_pr_tim(void);
#endif

#else
extern void  jcclib(unsigned msize, unsigned bsize, float  *A, double *r);
#endif

#if T800
struct {unsigned nproc, nn2, bsize, row_proc[Maxproc]; int opt_proc;} PRE;
struct {unsigned proc_no, nrow, ncol, rstart, rend, bcflag, end_of_mess;
	float ok;} OK;
#if DEBUGT800
struct {unsigned rstart, rend; flaot A[MaxOk]; double r[RHSsub];} REF;
#endif
/*struct {unsigned destin, end_of_mess;} LOAD;*/
struct {unsigned iter, end_of_iter_flag; double pAp, conv_test;} ITER;
unsigned nproc, rows_per_proc[Maxproc], bcflag;

CHAN **in_p, **out_p;

#if Timer
float proc_tim[2*Maxproc];
#endif

main(argc, argv, envp, in_ports, ins, out_ports, outs)

int argc, ins, outs;
char *argv[], *envp[];
CHAN *in_ports[], *out_ports[];

#else
int main(int argc, char **argv)
#endif

{
char buf[90];
char file_sav[30];
unsigned yes_no;

#if T800
in_p = in_ports;
out_p = out_ports;
#endif
/*-------------------- Read in the dat file ----------------------*/

arg_def = 0;               /* FALSE */
if(argc == 4) { 
arg_def = 1;               /* TRUE */
strcpy(file_sav,argv[1]);
}

rinput(file_sav);

if(arg_def == 0){
printf("Guass [1]  JCC [0]\n");
chkf(fgets(buf,90,stdin)); sscanf(buf,"%d",&gauss);
}else{
sscanf(argv[2],"%d",&gauss);
}

 /* Opening the data output file  for FEM and ADAPTIVITY*/

openfile_fp2(file_sav);

ste = (double *)calloc(ne*3, sizeof(double));
if(ste == NULL)
{printf(" out of core for ste\n");exit(1);}



/*----------------------- F. E. Program --------------------------- */

  openfile_fp3(file_sav); /* open wr file stress.dat */

  semi_band();       /* calc nsb */

  printf("\nSemi-Band-Width before = %d\n",nsb);
  
  if(arg_def == 0){
  printf("\n Semi-Band-Width reduction Yes[1] No[0]\n");
  chkf(fgets(buf,90,stdin)); sscanf(buf,"%d",&yes_no);
  } else {
  sscanf(argv[3],"%d",&yes_no);
  }

time(&t1);

  if(yes_no)
  setup(nn, ne, cord, nod);

#if T800
#else
#endif


  semi_band();       /* calc nsb */
  printf("\nSemi-Band-Width after = %d\n",nsb);

#if T800
  printf("\n*Optimum number of processors for this problem <= %d\n",
	 2*nn/nsb);
printf("Enter number of processors in the pipeline\n");
gets(buf); sscanf(buf,"%d",&nproc);

/*
if(nproc <= 2*nn/nsb) opt_proc = 1;
else opt_proc = 0;

while(4*nn/nproc < nsb) nproc--;
printf("problem requires %d processors\n", nproc);
*/

if(nproc < 3){
printf("This program needs atleast 3 processors but your\n");
exit(1);}

if(opt_proc == 1){ printf("\n\n***communication optimiztion is on***\n\n");
printf("if you want optimization to be 'on' enter [1] else [0]\n");
gets(buf);sscanf(buf,"%d",&opt_proc);
}
#endif



  nn2 = 2*nn;

#if Less
  ok = (float *)calloc(nn2*nsb, sizeof(float));
  if(ok == NULL)
  {printf(" out of core for ok\n");exit(1);}
#endif

  f = (double *)calloc(nn2, sizeof(double));
  if(f == NULL)
  {printf(" out of core for f\n");exit(1);}

  s = (double *)calloc(ne*6*3, sizeof(double));
  if(s == NULL)
  {printf(" out of core for s\n");exit(1);}

#if DEBUG
  print_orig_data();
#endif

  elem();

#if T800
#if Timer
  rec_pr_tim();
#endif
#endif

#if T800
#else
time(&t2);
#endif

fprintf(fp2,"\n\n*** Time for analysis = %f min or %ld sec ***\n",
       ((float)t2-(float)t1)/60., t2-t1);
fprintf(fp2,"*** no of iterations = %d ***\n", iter);
fprintf(stdout,"\n\n*** Time for analysis = %f min or %ld sec ***\n",
       ((float)t2-(float)t1)/60., t2-t1);
fprintf(stdout,"*** no of iterations = %d ***\n", iter);

printf("*** The stresses are stored in %s.str ***\n",file_sav);
printf("*** The results of FE  are stored in %s.fem ***\n", file_sav);

routput(file_sav);

fclose(fp2);
fclose(fp3);

return 0;

}



/*                           AO2.C                           */



/* Function openfile */
void openfile_fp2(char *file_sav)
{
char file[30];
strcpy(file,file_sav); strcat(file,".fem");
fp2=fopen(file,"w");
	if(fp2==(FILE *)NULL)
	{
	printf("cannot open %s\n",file);
	exit(1);
	}
return;
}

void openfile_fp3(char *file_sav)
{
char file[30];
strcpy(file, file_sav); strcat(file,".str");
fp3=fopen(file,"w");
	if(fp3==(FILE *)NULL)
	{
	printf("cannot open %s\n", file);
	exit(1);
	}
return;
}

/*** Finite Element function  ***/


void elem(void)
{
#if DEBUG
FILE *debug;
if((debug=fopen("debug.dat","w")) == NULL)
{printf("cannot open debug.dat\n"); exit(1);}
#endif

printf("\n F E analysis started with node = %d and elements = %d\n",nn,ne);

   /* zero s, ok and f */

#if Less
   zero_ok();
#endif

   /* put aload[] into f[] */
   aload_f();


#if T800
  /* send prelimenary data */

  distribute_rows_per_processor();

  preliminary_data();

#endif

  /* form stiff mat ok */
  printf("\n \n Assembling stiffness matrix \n");
  form_ok();


  /* apply boundary conds to the ok[][] */
  printf("applying boundary conditions\n");
  bc_ok();


#if T800
  printf("releasing worker tasks \n");
  ok_transmitter(0,0,0,0,1); /* release the workers */
#endif

#if T800
  printf("distributimg load vector\n");
  load_vector_distributer();
#endif

#if T800

#if DEBUGT800
  receive_and_print_reflected_data();
#endif

  t1 = ((float)timer_now())/((float)CLK);

  printf("receiving the iteration\n");
  receive_iter();
#endif



#if Less
#if DEBUG
  fprintf(debug,"%d %d\n",nsb-1, nn2);
  for(i = 0; i < nn2; i++)
  {
     fprintf(debug," %5d ",i+1);
     for(j = 0; j < nsb; j++)
     fprintf(debug," %12.5e",ok[nsb*i + j]);
     fprintf(debug," 0.0");                    /* solution if known*/
     fprintf(debug," %f", f[i]);               /* R.H.S. of the system */
     fprintf(debug,"\n");
  }
printf("###### DEBUG.dat generated #######\n");
fclose(debug);
#endif
#endif

  /* solution by gaussian elimination */
  printf("\n \n Solving stiffness matrix \n \n \n");


#if T800
#else
#endif

#if Less
  if(gauss)
  gauss_el();
  else
#endif
#if T800
  receive_data_f();
#else
  jcclib(nn2, nsb, ok, f);
#endif

#if T800
  t2 =  ((float)timer_now())/((float)CLK);
#else
#endif



  /* print the nodal displacements */
  fpr_defl();

  /* calculate and fprint the stresses */
  calc_stress();

return;
}

/* Function to zero back subs matrix s, overall stiffness matrix ok
 and load vector f */
#if Less
void zero_ok(void)
{
int i,j,k,l;
	for(i=0; i< ne; i++)
	{
		for(k=0; k< 3; k++)
		{
			for(l=0; l< 6; l++) s[18*i+6*k+l]=0.;
		}
	}

	for(i=0; i< 2*nn; i++)
	{
		for(j=0; j< nsb; j++) ok[nsb*i + j]=0.;
	f[i]=0.;
	}
return;
}
#endif

/* function to read in applied aload()  into f(array), must have load
card for last node even if no applied load */
void aload_f(void)
{
int i,ii,jj,k;

	if(iln >0 )
	{
	 printf("iln = %d\n",iln);
		for(i=0; i < iln; i++)
		{
		if(NCLOAD[i].nc > nn) break;
		ii=2*NCLOAD[i].nc-2;
		jj=2*NCLOAD[i].nc-1;
			for(k=0 ; k < 2; k++)
			{
			if(k==0) f[ii]=NCLOAD[i].aload[k];
			else
			f[jj]=NCLOAD[i].aload[k];
			}
		}
	}
return;
}


/* function to assemble the stiffness mat by adding the indivivual
   element matrices formed by function fst() */
void form_ok(void)
{
int i,n,ni,j,nj,ja,jb,k,l, nrow, ncol;
/* loop for each element */

	for(n=0; n < ne; n++)
	{
	fst(n);
		for(i=0; i< 3; i++)
		{
		ni=2*(nod[4*n + i]-1);
			for(j=0; j< 3; j++)
			{
			nj=2*(nod[4*n + j]-1);
				for(k=0; k< 2; k++)
				{
					for(l=0; l< 2; l++)
					{
					ja = i*2 + k;
					jb = j*2 + l;

					nrow = ni + k;
					ncol = nj + l;
					if(ncol < nrow) continue;
					ncol = ncol - nrow;
#if Less
			  ok[nsb*nrow + ncol] += est[ja][jb];
#endif

#if T800
			  bcflag = 0; /* this is  not boundary data */
			  ok_transmitter(nrow, ncol, ja, jb,  0);
#endif
						}
				}
			}
		}
	}

return;
}


/* function applies boundry conditions to stiffness matrix ok and check
each node in turn */
void bc_ok(void)
{
int ip,il,nrow,nri;
	for(ip=0; ip< nb; ip++)
	{
	il=nco[ip][0];
	nrow=(il*2)-2;
	/* check in  x direction */
	if(nco[ip][1] < 1) goto lab1;
#if Less
	ok[nsb*nrow + 0]=pow(10.,33.);
#endif
#if T800
	bcflag = 1; ok_transmitter(nrow,0,0,0,0);
#endif
	f[nrow]=0.;
lab1:
	/* check in y direction  */
	if(nco[ip][2] < 1) continue;
	nri=nrow+1;
#if Less
	ok[nsb*nri + 0]=pow(10.0,33.);
#endif
#if T800
	bcflag = 1; ok_transmitter(nri,0,0,0,0);
#endif
	f[nri]=0.;
	}
return;
}


/* function for the solution of ok[][] by gaussian elimination */
#if Less
void gauss_el(void)
{
int m1,mm,i,ii,j,k,jk, nba1,l;
double fac,sum;

	mm = nn2;
	m1 = mm - 1;
	for(i = 0 ; i < m1 ; i++)
	{
	ii = i + 1;
	nba1 = nsb + i ;
	  for(j = ii; j < nba1; j++)
	  {

	  if(j > mm - 1) break;
	  fac=ok[nsb*i + (j - i)]/ok[nsb*i + 0];

	      for(k = 0;k < nsb; k++)
	      {
	      l = k + (j - i);
	      if(l > nsb - 1) continue;
	      ok[nsb*j + k] += -fac*ok[nsb*i + l];
	      }
	  f[j]=f[j]-fac*f[i];
	  }
	}

/* Back subsitution  */

	f[mm-1]=f[mm-1]/ok[nsb*(mm-1) + 0];
	for(i = 0 ; i < m1; i++)
	{
	j = mm - i - 2;
	sum = 0.;
		for(jk = 1; jk < nsb; jk++)
		{
		if((j + jk) > mm - 1)  goto skip_2;
		sum +=ok[nsb*j + jk]*f[j + jk];
		}
	skip_2:
	f[j] =(f[j] -sum)/ok[nsb*j + 0];
	}
return;
}
#endif


/* function to print the nodal displacements */
void fpr_defl(void)
{
int i,l,m;

	m = 2*nn;

	fprintf(fp2,"\n \n        Results \
	   \n \n Deflections\n");
	fprintf(fp2,"   Node   x-def          y-def\n");

	for(i=0; i< m; i=i+2)
	{
	l=(i+2)/2;
	fprintf(fp2,"%5d   %15.8e %15.8e\n", l, f[i], f[i+1]);
	}
return;
}


/* function to calculate and fprint stresses */
void calc_stress(void)
{
int n,l,i,j,k;
double def[6], cc,aa,smax,smin,ang;
	fprintf(fp2," \n \n Stresses \n");

	fprintf(fp2,"Elm  Sx           Sy           Sxy \
	  Smax         Smin         Angle\n");

	for(n=0; n < ne; n++)
	{
		for(l=0 ;l < 3; l++)
		ste[3*n + l]=0. ;
		
		for(i=0; i < 3; i++)
		{
		def[i*2]=f[nod[4*n + i]*2-2];
		def[i*2+1]=f[nod[4*n + i]*2-1];
		}

		/* Multiply the stress matrix by def */

		for(j=0; j < 3; j++)
		{
			for(k=0; k< 6; k++)
			ste[3*n + j] += s[18*n+6*j+k]*def[k];
		}

		cc=(ste[3*n + 0]+ste[3*n + 1])/2. ;
		aa=sqrt(pow((ste[3*n + 0]-ste[3*n + 1]),2.)/4.+
		pow(ste[3*n + 2],2.));
		smax=cc+aa;
		smin=cc-aa;
		if(fabs(ste[3*n + 1] - smin) <= pow(10,-10.))
		ang = 90.0;
		else
		ang=57.29578*atan(ste[3*n + 2]/(ste[3*n + 1]-smin));
		fprintf(fp2,"%3d %12.5e %12.5e %12.5e %12.5e\
  %12.5e %11.4e\n",n+1,ste[3*n + 0],ste[3*n + 1],ste[3*n + 2],smax,smin,ang);
		fprintf(fp3,"%12.5e %12.5e %12.5e %12.5e %12.5e \n",
		ste[3*n + 0],ste[3*n + 1],ste[3*n + 2],smax,smin);
	
	}
return;
}                
 



			/* FUNCTION FST */

void fst(int n)
{
int n1,n2,n3,ll,i,j,k;
double v,e,xn1,xn2,xn3,yn1,yn2,yn3,comm,b1,b2,b3,c1,c2,c3,area,thic;
double b[3][6], d[3][3];        
/* This function works out the element stiffnesses of triangles */

	/* get nodes out of nod */
	
	n1=nod[4*n + 0]-1;
	n2=nod[4*n + 1]-1;
	n3=nod[4*n + 2]-1;
	
	/* get material properties and thichnesses */
	
	ll=nod[4*n + 3]-1;
	thic=amat[ll][2];
	v=amat[ll][0];
	e=amat[ll][1];

	/* define coordinates */

	xn1=cord[2*n1 + 0];
	xn2=cord[2*n2 + 0];
	xn3=cord[2*n3 + 0];

	yn1=cord[2*n1 + 1];
	yn2=cord[2*n2 + 1];
	yn3=cord[2*n3 + 1];

	/* Form elasticity matrix D */
	form_D(v,e,d);

	/* form strain displacement matrix b */

	b1=yn2-yn3;
	b2=yn3-yn1;
	b3=yn1-yn2;

	c1=xn3-xn2;
	c2=xn1-xn3;
	c3=xn2-xn1;

	area=(-(c2*b2)+(2.*(c1*b3))+(c3*b3)+(c1*b1))/2.;

	b[0][0]=b1/(2.*area);
	b[0][1]=0.;
	b[0][2]=b2/(2.*area);
	b[0][3]=0.;
	b[0][4]=b3/(2.*area);
	b[0][5]=0.;

	b[1][0]=0.;
	b[1][1]=c1/(2.*area);
	b[1][2]=0.;
	b[1][3]=c2/(2.*area);
	b[1][4]=0.;
	b[1][5]=c3/(2.*area);

	b[2][0]=c1/(2.*area);
	b[2][1]=b1/(2.*area);
	b[2][2]=c2/(2.*area);
	b[2][3]=b2/(2.*area);
	b[2][4]=c3/(2.*area);
	b[2][5]=b3/(2.*area);

	/* form back subsitution matrix */

	for(j=0; j< 6; j++)
	{
		for(i=0; i< 3; i++)
		{
			for(k=0; k< 3; k++)
			s[18*n+6*i+j] += d[i][k]*b[k][j];
		}
	}

	/* form element stiffness */

	for(i=0; i < 6; i++)
	{
		for(j = 0; j < 6; j++)
		{
		est[i][j]=0.;
			for(k=0; k < 3; k++)
			est[i][j]+=b[k][j]*s[18*n+6*k+i]*area*thic;
		}
	}
}
	

/* function to form elsticty matrix D */
void form_D(double v, double e, double d[][3])
{
double comm;
	comm=e/(1.-pow(v,2.));
	d[0][0]=comm;
	d[0][1]=comm*v;
	d[0][2]=0.0;

	d[1][0]=comm*v;
	d[1][1]=comm;
	d[1][2]=0.;

	d[2][0]=0.;
	d[2][1]=0.;
	d[2][2]=comm*(1.0-v)/2.;
return;
}

/* function to read in the data from the primary data file */

void rinput(char *file_sav)
{
char file[30];
FILE *fp;
int i,j,k,n;
double x1,y1,z1;
  
if(arg_def == 0){ 
   printf("Enter the input file (.dat assumed)\n");
   chkf(fgets(file,30,stdin)); 
   strcpy(file,trimwhitespace(file)); 
   strcpy(file_sav, file); 
   strcat(file,".dat");
   }else{
   strcpy(file,file_sav);
   strcat(file,".dat");
   }


   if((fp = fopen(file, "r")) == (FILE *) NULL)
   {printf("cannot open %s\n",file); exit(1);}
   chkf(fgets(buf,256,fp));
   fputs(buf,stdout);
   chkf(fgets(buf,256,fp));
   sscanf(buf,"%lf%lf",&dt,&damp);
   printf("%f %f \n",dt,damp);
   chkf(fgets(buf,256,fp));
   sscanf(buf,"%d %d %d %d %d %d %d %d %d %d %d",&nn,&ne,&nb,&nm,&nl,
   &ndime, &nstr, &nbe, &kmax,&nmax, &iln);

   chkf(fgets(buf,256,fp));
   sscanf(buf,"%d",&nq);


   if((cord = (double *)calloc(2*nn, sizeof(double))) == NULL)
   {printf(" out of core for cord\n");exit(1);}



   if((nod = (int *)calloc(4*ne, sizeof(int))) == NULL)
   {printf(" out of core for nod\n");exit(1);}


   for(i = 0; i < nn; i++)
   {
   chkf(fgets(buf,256,fp));
   sscanf(buf,"%d %lf %lf %lf",&n, &x1, &y1,&z1);
   cord[2*i] = x1;
   cord[2*i + 1] = y1;
   }



   if(ne > 0)
   {
      for(i = 0; i < ne; i++)
      {
      chkf(fgets(buf,256,fp));
      sscanf(buf,"%d %d %d %d %d",&n,&nod[4*i + 0],&nod[4*i + 1],
      &nod[4*i + 2],&nod[4*i + 3]);
      }
   }



   if(nl > 0)
   {
      for(i = 0; i < nl;i++)
      {
      chkf(fgets(buf,256,fp));
      }
   }



   if(nq > 0)
   {
      for(i = 0; i < nq; i++)
      {
      chkf(fgets(buf,256,fp));
      }
   }



   for(i = 0; i < nm; i++)
   {
   chkf(fgets(buf,256,fp));
   sscanf(buf,"%d %d %lf %lf %lf %lf %lf",&n,&mat[i],&amat[i][0],&amat[i][1],
   &amat[i][2],&amat[i][3],&amat[i][4]);
   }



   for(i = 0; i < nb; i++)
   {
   chkf(fgets(buf,256,fp));
   sscanf(buf,"%d %d %d %d",&nco[i][0],&nco[i][1],&nco[i][2],&nco[i][3]);
   }



   if(iln > 0)
   {
   iln = 0;
more:
   iln++;
   i = iln - 1;

   chkf(fgets(buf,256,fp));
   sscanf(buf,"%d %lf %lf %lf",&NCLOAD[i].nc,&NCLOAD[i].aload[0],&NCLOAD[i].aload[1],&NCLOAD[i].aload[2]);

   if(NCLOAD[i].nc < nn)
   goto more;
   }

   if(nstr > 0)
   {
   for(i = 0; i < nstr; i++)
   {
   chkf(fgets(buf,256,fp));
   chkf(fgets(buf,256,fp));
   }
   }
fclose(fp);
return;
}


/* function calculates the semi band width by finding the max diff
   between the node numbers */

   void semi_band(void)
   {
   int k,l,n1,n2,ndiff;

   nsb = 0;
   for(k = 0; k < ne; k++)
   {
      n1 = nod[4*k + 0];
      for(l = 2; l >= 0; l--)
      {
      n2 = nod[4*k + l];
      ndiff = (int)fabs((double)n1 - (double)n2);
      if(ndiff > nsb) nsb = ndiff;
      n1 = n2;
      }
   }
   nsb = (nsb + 1)*2;
   return;
   }



void routput(char *file_sav)
{
char filename[30];
size_t len;
int i,j,k,n;
double tL[3],rL;
FILE *fp0;

   strcpy(filename, file_sav);
   len = strlen(filename);
   filename[len-1] = '_';
   strcat(filename,".dat");

   if((fp0 = fopen(filename,"w")) == (FILE *) NULL)
   {printf("cannot open %s\n",filename); exit(1);}

   fprintf(fp0,"out.dat\n");

   fprintf(fp0,"%10.5f%10.5f\n", dt, damp);

   fprintf(fp0,"%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d\n",nn,ne,nb,nm,nl,
   ndime, nstr, nbe, kmax,nmax, iln);

   fprintf(fp0,"%5d\n",nq);


   for(i = 0; i < nn; i++)
   {
   fprintf(fp0,"%5d%15.5f%15.5f\n",i+1, cord[2*i+0], cord[2*i+1]);
   }



   if(ne > 0)
   {
      for(i = 0; i < ne; i++)
      {
      fprintf(fp0,"%5d%5d%5d%5d%5d\n",i+1,
      nod[4*i+0],nod[4*i+1],nod[4*i+2],nod[4*i+3]);
      }
   }






   for(i = 0; i < nm; i++)
   {
   fprintf(fp0,"%5d%5d%12.5e%12.5e%12.5e%12.5e%12.5e\n",i+1,mat[i],amat[i][0],
   amat[i][1],amat[i][2],amat[i][3],amat[i][4]);
   }



   for(i = 0; i < nb; i++)
   {
   fprintf(fp0,"%5d%5d%5d%5d\n",nco[i][0],nco[i][1],nco[i][2],nco[i][3]);
   }



   if(iln > 0)
   {
   iln = 0;
more:
   iln++;
   i = iln - 1;

   fprintf(fp0,"%5d%15.5f%15.5f%15.5f\n",NCLOAD[i].nc,NCLOAD[i].aload[0],NCLOAD[i].aload[1],
   NCLOAD[i].aload[2]);

   if(NCLOAD[i].nc < nn)
   goto more;
   }


printf("*** renumbered input data file %s generated ***\n", filename);
fclose(fp0);
free(ste);
free(ok);
free(f);
free(s);
free(cord);
free(nod);
return;
}

char* trimwhitespace(char *str)
{
  char *end;

  /* Trim leading space */
  while(isspace((unsigned char)*str)) str++;

  if(*str == 0)  /* All spaces? */
    return str;

  /* Trim trailing space */
  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)*end)) end--;

  /* Write new null terminator character */
  end[1] = '\0';

  return str;
}

void chkf(char *str){
   if( str == NULL){
     printf("line not read\n");
     exit(1);
   }
}

#if Less
void print_orig_data(void)
{
FILE *orig;
register int i, j;

if((orig = fopen("orig.dat","w")) == NULL)
{printf("cannot open orig.dat\n"); exit(1);}

fprintf(orig,"FILE *** orig.dat ***\n");

for(i = 0; i < nn2; i++)
{
   fprintf(orig,"%5d", i+1);
   for(j = 0; j < nsb; j++)
   fprintf(orig,"%12.5e", ok[nsb*i+j]);
   fprintf(orig,"%12.5e", f[i]);
   fprintf(orig,"\n");
}
fclose(orig);
printf("### orig.dat generated ###\n");
}
#endif



#if T800
void preliminary_data(void)
{
unsigned msglen;

msglen = sizeof(PRE);
PRE.nn2 = nn2;
PRE.nproc = nproc;
PRE.bsize = nsb;
PRE.opt_proc = opt_proc;
memcpy(PRE.row_proc, rows_per_proc, sizeof(rows_per_proc));
chan_out_message(msglen, &PRE, out_p[2]);
}


void ok_transmitter(int nrow, int ncol, int ja, int jb, int flag)
{
register int i;
unsigned target_proc, sum, msglen, rstart, rend, vali;
float   val;

msglen = sizeof(OK);



sum = 0;
for(target_proc = 0, i = 0; i < nproc; i++, target_proc++)
{
sum += rows_per_proc[i];

if(nrow < sum)
{
vali=OK.rstart = sum - rows_per_proc[i];        /* global start no of rows */
vali=OK.rend   = sum;                           /* global end no of the rows */
vali=OK.proc_no = target_proc;
vali=OK.nrow = nrow;
vali=OK.ncol = ncol;

if(bcflag == 1)                                 /* bcflag [1] indicates boundary condition */
val =OK.ok = 1.e+33;
else
val =OK.ok = est[ja][jb];

vali=OK.end_of_mess = flag;                     /* 1 marks end of communication */
chan_out_message(msglen, &OK, out_p[2]);
break;
}
}

}


#if DEBUGT800
void receive_and_print_reflected_data(void)
{
FILE *reflect;
register int knt, i,ip,j;
unsigned  msglen, rstart, rend;

msglen = sizeof(REF);

if((reflect = fopen("reflect.dat","w")) == NULL)
{printf("cannot open reflect.dat\n"); exit(1);}

fprintf(reflect," FILE *** reflect.dat ***\n");

for(knt = 0; knt < nproc; knt++)
{
alt_wait(1, in_p[2]);
chan_in_message(msglen,&REF, in_p[2]);




rstart = REF.rstart;
rend   = REF.rend;


fprintf(reflect,"reflected data received from processor no. %d\n", knt);
fprintf(reflect,"from rows no. %d to row no. %d\n", rstart, rend);
/*fprintf(stdout, "reflected data received from processor no. %d\n", knt);*/

for(ip = 0, i = rstart; i < rend; i++, ip++)
{
  fprintf(reflect,"%5d", i+1);

  for(j = 0; j < nsb; j++)
  fprintf(reflect,"%12.5e", REF.A[nsb*ip+j]);
  fprintf(reflect,"%12.5e", REF.r[ip]);
  fprintf(reflect,"\n");
}
}

printf("### reflect.dat generated ###\n");
fclose(reflect);
}
#endif


void distribute_rows_per_processor(void)
{
register int i;
unsigned knt;
knt = 0;
for(i = 1; i < nproc; i++)
{
rows_per_proc[i] = nn2 / nproc;
knt += rows_per_proc[i];
}
rows_per_proc[0] += nn2 - knt;
}

void load_vector_distributer(void)
{
int i;
unsigned sum, index, msglen;
int sizeof_f;

sum = 0;
for( i = 0; i < nproc; i++)
{

sum += rows_per_proc[i];


index = sum - rows_per_proc[i];
sizeof_f = sizeof(double)*rows_per_proc[i];
chan_out_word(i, out_p[2]);
chan_out_word(sizeof_f, out_p[2]);
chan_out_message(sizeof_f, (char *) &f[index], out_p[2]);
}
}


void receive_iter(void)
{
register int knt;
int msglen, dum;

msglen = sizeof(ITER);

do{
/*for(knt = 0; knt < nproc; knt++)
{*/
alt_wait(1, in_p[2]);
/*chan_in_word(&dum, in_p[2]); printf(" *** \n");*/
chan_in_message(msglen,&ITER, in_p[2]);
iter = ITER.iter;
/*if(knt == nproc-1)*/
printf("%d %12.5e %12.5e\n", ITER.iter, ITER.pAp, ITER.conv_test);
/*} */
}
while(ITER.end_of_iter_flag == 0);

}

void  receive_data_f(void)
{
int msglen, dum;

alt_wait(1, in_p[2]);
chan_in_word(&dum,    in_p[2]);
chan_in_word(&msglen, in_p[2]);
chan_in_message(msglen,(char *) f, in_p[2]);
}

#if T800
#if Timer

void rec_pr_tim(void)
{
register int i;
float tt3;
int msglen;
float tt1, tt2, t1_sum, t2_sum;

t1_sum = 0;
t2_sum = 0;

msglen = 2*nproc*sizeof(float);

alt_wait(1, in_p[2]);
chan_in_message(msglen, (char *) proc_tim, in_p[2]);

for(i = 0; i < nproc; i++)
{

tt1=proc_tim[2*i+0];  tt2=proc_tim[2*i+1]; tt3 = tt1/tt2;

fprintf(stdout,"proc_no %3d t_cal=%7.2f t_com=%7.2f t_cal/t_com = %7.2f\n",
	i, tt1, tt2, tt3);
fprintf(fp2   ,"proc_no %3d t_cal=%7.2f t_com=%7.2f t_cal/t_com = %7.2f\n",
	i, tt1, tt2, tt3);
t1_sum += tt1;   t2_sum += tt2;
}
fprintf(stdout,"t_calc_avg (s) =%7.2f t_comm_avg =%7.2f t_calc_avg/t_comm_avg=%7.2f\n",
		t1_sum/((float)nproc), t2_sum/((float)nproc), t1_sum/t2_sum);
fprintf(fp2,"t_calc_avg (s) =%7.2f t_comm_avg =%7.2f t_calc_avg/t_comm_avg =%7.2f\n",
		t1_sum/((float)nproc), t2_sum/((float)nproc), t1_sum/t2_sum);

}
#endif
#endif
#endif
