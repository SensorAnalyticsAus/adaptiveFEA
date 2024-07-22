#define II 15000

#define JJ 20000

#define LL 50

#define MM 25

#define NN 30

#define LIL 10

#define BUFF 256


#include<stdio.h>
#include<math.h>

void rinput();
void routput();
void tri_len(int i,double tL[]);  
void cabl_len(int i, double *prL); 
void seglen(double x1,double y1, double x2, double y2, double *pL);   
void chkf(char *str);

FILE *fp1;

 char buf[BUFF];

 char titel[32];

 int nn,ne,nb,nm,nl,ndime,nstr,nbe,kmax,nmax,iln,nq;

 int nod[JJ][4],nodl[LIL][3],nodq[LIL][5], nco[MM][4],nsnod[LIL][50],
     sort[II][20], mat[20],nc[30],ngs[LIL];

 double para1, para2;
 double dt,damp;

 double cord[3*II];

 double amat[20][10],aload[30][3];



main(int argc, char *argv[])
{
register i,ii;

if(argc == 1) {
fprintf(stdout,"usage: filename.dat  nb bn1 ... bnn iln ln1 Fx1 Fy1 lnn ...\
Fxn Fyn\n");
exit(1);
}

if((fp1=fopen(argv[1],"r"))==(FILE *)NULL)
{printf("cannot open %s\n",argv[1]);exit(1);}   

rinput();

/* read the boundary nodes and fix in x,y, and z directions */
sscanf(argv[2],"%d",&nb);

if(nb > 0)
 if(nb < nn)
 {
 for(i = 0; i < nb; i++)
   {
   sscanf(argv[i+3],"%d",&nco[i][0]);
   nco[i][1] = 1; nco[i][2] = 1; nco[i][3] = 1;
   }
 }else{
 fprintf(stdout,"bad boundary node number argument %d\n",nb);
 exit(1);
 }

/* read no of loaded nodes, node Fx, Fy, Fz=0, add last node, if not
   already specified as zero loaded */
sscanf(argv[3+nb],"%d",&iln);
for(i = 0; i < iln; i++)
 {
  ii = (4+nb)+i*3;
  sscanf(argv[ii],"%d",&nc[i]);
  sscanf(argv[ii+1],"%lf",&aload[i][0]); 
  sscanf(argv[ii+2],"%lf",&aload[i][1]); 
  aload[i][2] = 0.;
 }

if(nc[iln-1] < nn) /* if last loaded node is the last node */
{
nc[iln] = nn;
aload[iln][0] = 0.; aload[iln][1]=0.; aload[iln][2]=0.;
iln += 1; 
}

routput();
}

/* function to read in the data from the primary data file */

void rinput()
{
int i,n;
double val,x1,y1,z1;


   chkf(fgets(titel,BUFF,fp1));
   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%lf%lf",&para1,&para2);
   sscanf(buf,"%lf%lf",&dt,&damp);
   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%d %d %d %d %d %d %d %d %d %d %d",&nn,&ne,&nb,&nm,&nl,
   &ndime, &nstr, &nbe, &kmax,&nmax, &iln);

   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%d",&nq);



   for(i = 0; i < nn; i++)
   {
   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%d %lf %lf %lf",&n, &x1, &y1,&z1);
   val = cord[3*i] = x1;
   val = cord[3*i + 1] = y1;
   val = cord[3*i + 2] = z1;
   }



   if(ne > 0)
   {
      for(i = 0; i < ne; i++)
      {
      chkf(fgets(buf,BUFF,fp1));
      sscanf(buf,"%d %d %d %d %d",&n,&nod[i][0],&nod[i][1],
      &nod[i][2],&nod[i][3]);
      }
   }



   if(nl > 0)
   {
      for(i = 0; i < nl;i++)
      {
      chkf(fgets(buf,BUFF,fp1));
      sscanf(buf,"%d %d %d %d",&n,&nodl[i][0],&nodl[i][1],
      &nodl[i][2]);
      }
   }



   if(nq > 0)
   {
      for(i = 0; i < nq; i++)
      {
      chkf(fgets(buf,BUFF,fp1));
      sscanf(buf,"%d %d %d %d %d %d",&n,&nodq[i][0],&nodq[i][1],
      &nodq[i][2],&nodq[i][3],&nodq[i][4]);
      }
   }



   for(i = 0; i < nm; i++)
   {
   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%d %d %lf %lf %lf %lf %lf",&n,&mat[i],&amat[i][0],&amat[i][1],
   &amat[i][2],&amat[i][3],&amat[i][4]);
   }



   for(i = 0; i < nb; i++)
   {
   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%d %d %d %d",&nco[i][0],&nco[i][1],&nco[i][2],&nco[i][3]);
   }



   if(iln > 0)
   {
   iln = 0;
more:
   iln++;
   i = iln - 1;

   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%d %lf %lf %lf",&nc[i],&aload[i][0],&aload[i][1],&aload[i][2]);

   if(nc[i] < nn)
   goto more;
   }

   if(nstr > 0)
   {
   for(i = 0; i < nstr; i++)
   {
   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%d",&ngs[i]);
   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%d %d %d %d %d %d %d %d",&nsnod[i][0],&nsnod[i][1],
   &nsnod[i][2],&nsnod[i][3],&nsnod[i][4],&nsnod[i][5],&nsnod[i][6],
   &nsnod[i][7]);
   }
   }
return;
}

/* function to write in the data from the primary data file */

void routput()
{
int i;
double tL[3],rL;

   fprintf(stdout,"%s",titel);
   fprintf(stdout,"%10.5f%10.5f\n",para1,para2);
   fprintf(stdout,"%6d%6d%5d%5d%5d%5d%5d%5d%5d%5d%5d\n",nn,ne,nb,nm,nl,
   ndime, nstr, nbe, kmax,nmax, iln);

   fprintf(stdout,"%5d\n",nq);



   for(i = 0; i < nn; i++)
   {
   fprintf(stdout,"%5d%15.5f%15.5f%15.5f\n",i+1, cord[3*i + 0], cord[3*i + 1],
   cord[3*i + 2]);
   }



   if(ne > 0)
   {
      for(i = 0; i < ne; i++)
      {
      tri_len(i, tL);
      fprintf(stdout,"%5d%5d%5d%5d%5d%15.8f%15.8f%15.8f\n",i+1,nod[i][0],nod[i][1],
      nod[i][2],nod[i][3],tL[0],tL[1],tL[2]);
      }
   }



   if(nl > 0)
   {
      for(i = 0; i < nl;i++)
      {
      cabl_len(i, &rL);
      fprintf(stdout,"%5d%5d%5d%5d%15.8f\n",i+1,nodl[i][0],nodl[i][1],
      nodl[i][2],rL);
      }
   }



   if(nq > 0)
   {
      for(i = 0; i < nq; i++)
      {
      fprintf(stdout,"%5d%5d%5d%5d%5d%5d",i+1,nodq[i][0],nodq[i][1],
      nodq[i][2],nodq[i][3],nodq[i][4]);
      }
   }



   for(i = 0; i < nm; i++)
   {
   fprintf(stdout,"%5d%5d%12.5e%12.5e%12.5e%12.5e%12.5e\n",i+1,mat[i],amat[i][0],
   amat[i][1],amat[i][2],amat[i][3],amat[i][4]);
   }



   for(i = 0; i < nb; i++)
   {
   fprintf(stdout,"%5d%5d%5d%5d\n",nco[i][0],nco[i][1],nco[i][2],nco[i][3]);
   }



   if(iln > 0)
   {
   iln = 0;
more:
   iln++;
   i = iln - 1;

   fprintf(stdout,"%5d%15.5f%15.5f%15.5f\n",nc[i],aload[i][0],aload[i][1],
   aload[i][2]);

   if(nc[i] < nn)
   goto more;
   }

   if(nstr > 0)
   {
   for(i = 0; i < nstr; i++)
   {
   fprintf(stdout,"%5d\n",ngs[i]);
   fprintf(stdout,"%5d%5d%5d%5d%5d%5d%5d%5d\n",nsnod[i][0],nsnod[i][1],
   nsnod[i][2],nsnod[i][3],nsnod[i][4],nsnod[i][5],nsnod[i][6],
   nsnod[i][7]);
   }
   }
fclose(stdout);
return;
}

/* function to calculate the lengths of the three sides of each element
for storage in tlen[] */
void tri_len(int i,double tL[])
{
int j,k;
double x1,y1,x2,y2,L;

   k = 0;
      for(j = 0; j < 3; j++)
      {
      k = j + 1;
      if(k > 2)
      k = 0;
      x1 = cord[3*(nod[i][j] - 1) + 0];
      y1 = cord[3*(nod[i][j] - 1) + 1];
      x2 = cord[3*(nod[i][k] - 1) + 0];
      y2 = cord[3*(nod[i][k] - 1) + 1];
      seglen(x1,y1,x2,y2,&L);
      tL[j] = L;
      k++;
      }
return;
}                                

/* function to calculate the length of the cable segments */
void cabl_len(int i, double *prL)
{
double x1,y1,x2,y2,L;

   x1 = cord[3*(nodl[i][0] - 1) + 0];
   y1 = cord[3*(nodl[i][0] - 1) + 1];
   x2 = cord[3*(nodl[i][1] - 1) + 0];
   y2 = cord[3*(nodl[i][1] - 1) + 1];
   seglen(x1,y1,x2,y2,&L);
   *prL = L;
return;
}

/* Function to calculate the length of a segment */
void seglen(double x1,double y1, double x2, double y2, double *pL)
{
double L;
L=sqrt(pow(x2-x1,2.)+pow(y2-y1,2.));
*pL = L;
return;
}

void chkf(char *str){
   if( str == NULL){
     printf("line not read\n");
     exit(1);
   }
}                   
