			 /* aom.c*/

  /* Master main() */


#define SUBDIV 0

#include <stdio.h>

#include <stdlib.h>

#include <string.h>

#include <math.h>

#include <time.h>

#include <ctype.h>

#include "incl.c"


IN ab;
OUT cd;




void aow(int i_i);
void load_info(int *pi);
void rinput();
void r_node_mesh(void);
void node_para(void);
void un_load_st(void);
void end_load(void);
void tri_len(int i,double tL[]);
void cabl_len(int i, double *prL);
void openfile1(void);
void openw_info(char *s);
void search_node(void);
void routput(char *file_out);
void stack(void);
int  comp_st(int node);
void revise_nod(int trace[], int *pjj);
void nb_iln(int n_iln);
void min_he(double he[],double *phi_min, int *ppos);
void max_he(double he[],double *phi_max, int *ppos);
void smooth1(void);
int  check_bound(int i);
void load_cord(double x[], double y[] ,int i, int k);
void interdiag(void);
void nodecord(int node, double *px, double *py);
void elemcord(int n,int nn[], double x[], double y[]);
void sidelen(double x[], double y[], double L[]);
void maxlen(double L[], int n, double *pL, int *pn);
void lenod(int no, int n[]);
int  comp(int n1[], int n[]);
void swap(int n1[], int n[], int e1, int e);
void minang(int n1[], int n[], double *pmA);
void calcang(double x1[],double y1[], double x[], double y[], double A[]);
void dot(double x1, double y1, double x2, double y2,
     double x3, double y3, double x4, double y4, double *pD);
void msm(int node, int no, int nen[]);
void elimcom(int node, int *pkk, int nrn[]);
void avgcom(int node,int kk, int nrn[]);
int  interseg(double x1,double y1,double x2,double y2,
     double x3,double y3,double x4,double y4, int *ptch, int *pwth, 
     int *near_parr);
int  parallel(double x1, double y1, double x2, double y2,
	     double x3, double y3, double x4, double y4);
void seglen(double x1,double y1, double x2, double y2, double *pL);
void store_bg(int *pi);
void scan_div();
int  cmp_side_hhi(int i);
void p_stack(int ni);
int ck_small_hhi(int i, double *phhi_avg);
#if SUBDIV
void subdivide(int i, double hhi_avg);
#endif
void ex_bg();
void minang_1(int n1[], double *pmA);
void calcang_1(double x1[],double y1[], double A[]);
void  bound_send(void);
void  magnify(void);   /* for small triangle sizes */
void  demagnify(void); /* 5 April 96 */
char* trimwhitespace(char *str); /*https://stackoverflow.com/a/122721/4106464*/
void chkf(char *str);

int main(int argc, char **argv)
{
char file_sav[BUFF], file_out[BUFF], buf[LINEBUFF]; 
int  pos, meshsel, bc[LL];
int i;

/*buf = (char *) malloc(sizeof(char) * LINEBUFF);
  if(buf == NULL) 
    { 
        printf("Memory allocation failed!"); 
        exit(1); 
    } */

if(argc > 1 && argc != 4){
printf("usage filename 1 1 \n");
exit(1);
}

/*-------------------- Read in the dat file ----------------------*/

/*if((cord = (double *)calloc(3*II, sizeof(double))) == NULL)
{printf(" out of core for cord\n");exit(1);}*/

printf("Max of %d nodes and %d elements can be meshed\n\n", nn = II, ne = JJ);

printf("***Node para are going to be used*** \n"); 

if(argc != 4){
printf("Enter the project file name (without extension)\n");
chkf(fgets(buf,LINEBUFF,stdin));
strcpy(buf,trimwhitespace(buf));
} else {
strcpy(buf,argv[1]);
}
strcpy(file_sav, buf); 
strcat(buf,".dat");
if((fp1=fopen(buf,"r"))==(FILE *)NULL)
{printf("cannot open %s\n",buf);exit(1);}

strcpy(file_out, file_sav); strcat(file_out, "o.dat");
if((fp0=fopen(file_out,"w"))==(FILE *)NULL)
{printf("cannot open %s\n",file_out);exit(1);}
rinput();
fclose(fp1);
fclose(fp0);

/*printf("Enter the ratio of avg_parameter/avg_side_length\n");
gets(buf);sscanf(buf,"%lf",&Fac);*/

Fac = 0.2;

strcpy(buf, file_sav); strcat(buf,".men");
if((hhi = (double *)calloc(II,sizeof(double))) == NULL)
{printf(" out of core for hhi\n");exit(1);}

printf("opening %s for reading \n", buf);
if((fp5=fopen(buf,"r"))==(FILE *)NULL)
{
 printf("cannot open %s  opening %s.me for read\n",buf, file_sav); 
 strcpy(buf, file_sav); strcat(buf, ".me");
 if((fp5=fopen(buf,"r"))==(FILE *)NULL)
 {printf("cannot open %s \n",buf); exit(1);}
 
 if((hi = (double *)calloc(JJ,sizeof(double))) == NULL)
 {printf(" out of core for hi\n");exit(1);}
 node_para();
}
else
{
r_node_mesh();
}

magnify();

fclose(fp5);


/*-----------------------------------------------------------------*/
if(argc != 4){
printf("Enter [1] if coarse background mesh is to be refined\n");
/*gets(buf);*/
chkf(fgets(buf,LINEBUFF,stdin));
sscanf(buf,"%d",&meshsel);
} else {
meshsel = 1;
}

if(argc != 4) {
printf("Enter [1] if mesh post-processing required\n");
chkf(fgets(buf,LINEBUFF,stdin));
sscanf(buf,"%d",&pos); 
} else {
pos =  1;
}

if(meshsel == 1)
{
meshsel = 0; /* zeroed for latter use */

/* ---------------- dynamic allocator block for Mesh Generator---------- */

if((nod_st = (int *)calloc(4,sizeof(int))) == NULL)
{printf(" out of core for nod_st\n");exit(1);}

if((cord_st = (double *)calloc(3*3, sizeof(double))) == NULL)
{printf(" out of core for cord_st\n");exit(1);}



/* ----------------- Mesh Generator Block ------------------------- */


       printf("The program is to be run without debug\n");

       time(&tima);             /* begin time count for pre-processing */


       nodeb = 1;
       
       nn1 = 0;
       ne1 = 0;
       nb1 = 0;


       store_bg(&n_iln); /* store bg mesh */

       if((ne_ii = (int *)calloc(ne_bg, sizeof(int))) == NULL)
       {printf(" out of core for ne_ii\n");exit(1);}

       if((id = (int *)calloc(ne_bg, sizeof(int))) == NULL)
       {printf(" out of core for id\n");exit(1);}


       if((ratio_ii = (double *)calloc(ne_bg, sizeof(double))) == NULL)
       {printf(" out of core for ratio_ii\n");exit(1);}



       /*ex_bg();  debug examination */

       printf("\n The No. of p_stack elems = %d\n",knt_bg);

#if SUBDIV
       printf("\n The No. of subdiv elems = %d\n", knt_div);
#endif

     openw_info(file_sav);
     bound_send();
     fclose(fp7);

     ne = ne1;
     nn = nn1;

     if(pos == 1){  
     printf("\n over all diagonal exchange started \n");
     interdiag();

     search_node();

     printf("\n overall smoothing started\n");
     for(i = 0; i < 5; i++)
     smooth1(); 
     }

printf("iln = %d\n",iln);
nb_iln(iln);

time(&timb);


printf("\n Time taken for mesh generation = %10.4f min\n", ((float)timb - (float)tima)/60.);
printf("\n\n\n  *** Mesh Compiled with %d Nodes and %d Element *** \n",nn , ne);



free((char *)nod_st);

free((char *)cord_st);


free((char *)hhi);


}


/* ---------------------- input.dat file generation ---------------- */

end_load();
demagnify();
routput(file_out);
printf("*** File %s with %d nodes & %d elements has been generated ***\n"
      ,file_out, nn, ne);

printf("*** Neural net training data saved in %s.inf ***\n", file_sav);
/*------------------------------------------------------------------- */

return(0);
}






/*
void   un_load_st(void)
   {
   nn = nn1;
   ne = ne1;


     nb_iln(n_iln);
}
*/


/* function to add a max node value to nc[] to terminate
   load reading condition for FE analysis */
void end_load(void)
{
int kt;
   /* count the no. of entries in nc[] */
   for(kt=0; nc[kt] != 0; kt++);
   
   nc[kt] = II;  /* take last node as II */

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





/* Function openfile1 */
void openfile1(void)
{
if(name != 1){
printf("Enter the mesh file\n");
chkf(fgets(buf,LINEBUFF,stdin));
sscanf(buf,"%s",nout1);
fp5=fopen(nout1,"r");
	if(fp5==(FILE *)NULL)
	{
	printf("cannot open %s\n",nout1);
	exit(1);
	}
}else
{
fp5=fopen("me.dat","r");
	if(fp5==(FILE *)NULL)
	{
	printf("cannot open %s\n",nout1);
	exit(1);
	}
}
return;
}



void openw_info(char *s)
{
char fils[30];
strcpy(fils, s);
strcat(fils, ".inf");
fp7=fopen(fils,"w");
	if(fp7==(FILE *)NULL)
	{
	printf("cannot open %s\n",fils);
	exit(1);
	}
return;
}




/*** Function to identify common nodes in the elements ***/

void search_node(void)
{
   int i,j,k,l,m;
   for(i=0; i < nn; i++)
      for(j=0; j < 20; j++)
       sort[i][j]=0;

/* Begin searching and sorting the elements connecting to a node,
at a time */
   for(i=0; i < nn; i++)
   {
    l=0;
      for(j=0; j < ne; j++)
      {
	  for(k=0; k < 3; k++)
	  {
	     if((m=i+1) == nod[j][k])
	     {
	      sort[i][l]= j+1;
	      l++;
	      break;
	      }
	  }
      }
   }              

return;
}



/* function to read in the data from the primary data file */

void rinput()
{
int i,n;
double val,x1,y1,z1;


   chkf(fgets(buf,BUFF,fp1));
   fputs(buf,stdout);
   fputs(buf,fp0);
   chkf(fgets(buf,BUFF,fp1));
   sscanf(buf,"%lf%lf",&dt,&damp);
   printf("%f %f \n",dt,damp);
   fprintf(fp0,"%10.5f%10.5f\n",dt,damp);
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

void routput(char *file_out)
{
int i;
double tL[3],rL;

if((fp0=fopen(file_out,"a"))==(FILE *)NULL)
{printf("cannot open %s \n", file_out);exit(1);}

   fprintf(fp0,"%6d%6d%5d%5d%5d%5d%5d%5d%5d%5d%5d\n",nn,ne,nb,nm,nl,
   ndime, nstr, nbe, kmax,nmax, iln);

   fprintf(fp0,"%5d\n",nq);



   for(i = 0; i < nn; i++)
   {
   fprintf(fp0,"%5d%15.5f%15.5f%15.5f\n",i+1, cord[3*i + 0], cord[3*i + 1],
   cord[3*i + 2]);
   }



   if(ne > 0)
   {
      for(i = 0; i < ne; i++)
      {
      tri_len(i, tL);
      fprintf(fp0,"%5d%5d%5d%5d%5d%15.8f%15.8f%15.8f\n",i+1,nod[i][0],nod[i][1],
      nod[i][2],nod[i][3],tL[0],tL[1],tL[2]);
      }
   }



   if(nl > 0)
   {
      for(i = 0; i < nl;i++)
      {
      cabl_len(i, &rL);
      fprintf(fp0,"%5d%5d%5d%5d%15.8f\n",i+1,nodl[i][0],nodl[i][1],
      nodl[i][2],rL);
      }
   }



   if(nq > 0)
   {
      for(i = 0; i < nq; i++)
      {
      fprintf(fp0,"%5d%5d%5d%5d%5d%5d",i+1,nodq[i][0],nodq[i][1],
      nodq[i][2],nodq[i][3],nodq[i][4]);
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

   fprintf(fp0,"%5d%15.5f%15.5f%15.5f\n",nc[i],aload[i][0],aload[i][1],
   aload[i][2]);

   if(nc[i] < nn)
   goto more;
   }

   if(nstr > 0)
   {
   for(i = 0; i < nstr; i++)
   {
   fprintf(fp0,"%5d\n",ngs[i]);
   fprintf(fp0,"%5d%5d%5d%5d%5d%5d%5d%5d\n",nsnod[i][0],nsnod[i][1],
   nsnod[i][2],nsnod[i][3],nsnod[i][4],nsnod[i][5],nsnod[i][6],
   nsnod[i][7]);
   }
   }
fclose(fp0);
return;
}



/* function arranges the output nod[] and cord[] from each trian call
   into stack arrays st_nod[] and st_cord[] */
void stack(void)
{
int i,ii,jj,test1;
int trace[II],vali;
double val;
   ii = nn1;
   jj = ne1;

   for(i = 0; i < 3; i++)
   {

     test1 = 0;
     test1 = comp_st(i); /* test1 returns the compared node_st */  

     if(test1 == -1)
     {
     val = cord[3*ii + 0] = cord_st[3*i + 0];
     val = cord[3*ii + 1] = cord_st[3*i + 1];
     vali = trace[i] = ii + 1;
     ii++;
     continue;
     }
     vali = trace[i] = test1;
  }

   /* load in the new topology */

   revise_nod(trace,&jj);
     

  nn1 = ii;
  ne1 = jj; 

return;
}

/* function to compare the packet nodes with the stack nodes */
int comp_st(int node)
{
int i;
double x,y,x_st,y_st;
double tol = pow(10.,-4.);

   x = cord_st[3*node + 0];
   y = cord_st[3*node + 1];

   for(i = 0; i < nn1; i++)
   {
   x_st = cord[3*i + 0];
   y_st = cord[3*i + 1];
      if((fabs(x - x_st) <= tol) && (fabs(y - y_st) <= tol))
      return(i + 1);
   }
return(-1);
}

/* function to revise the node topology */

void revise_nod(int trace[], int *pjj)
{
int vali;
int jj = *pjj;
   vali = nod[jj][0] = trace[0];
   vali = nod[jj][1] = trace[1];
   vali = nod[jj][2] = trace[2];
   vali = nod[jj][3] = nod_st[3];
   jj++;

*pjj = jj;
return;
}

/* function to identify the restrained and loaded nodes from the
   compiled mesh */

void nb_iln(int n_iln)
{
int i,j,k,l,m;
double tol = pow(10.,-6.);
double cx_b,cy_b,cx_l,cy_l;
   k = 0;
   m = 0;
   nb = 0;

   if(nc[iln-1] > nn) nc[iln-1] = II;


   for(i = 0; i < nn; i++)
   {

      for(j = 0; j  < nb_bg; j++)
      {
       cx_b = fabs(cord[3*i + 0] - cord_bg[nco_bg[j] - 1][0]);
       cy_b = fabs(cord[3*i + 1] - cord_bg[nco_bg[j] - 1][1]);
       if(cx_b <= tol && cy_b<= tol)
       {
	nco[k][0] = i + 1;
	nco[k][1] = 1;
	nco[k][2] = 1;
	nco[k][3] = 1;
	k++;
	break;
       }
      }

       for(l = 0;l < n_iln; l++)
       {
	if(nc_bg[l] > nn ) break;
	cx_l = fabs(cord[3*i + 0] - cord_bg[nc_bg[l] - 1][0]);
	cy_l = fabs(cord[3*i + 1] - cord_bg[nc_bg[l] - 1][1]);
	if(cx_l <= tol &&  cy_l <= tol)
	{
	 nc[m] = i + 1;
	 m++;
	 break;
	}
       }           
   }

   nb = k;
return;
}
 


/* function to find the min of hi[] */
void min_he(double he[],double *phi_min, int *ppos)
{
int i;

   *ppos = 0;

   *phi_min = pow(10.,30.);

   for(i = 0; i < ne; i++)
   {

   if(*phi_min > he[i])
   {
   *phi_min = he[i];

   *ppos = i;
   }

   }
return;
}


/* function to find the max of hi[] */
void max_he(double he[],double *phi_max, int *ppos)
{
int i;

   *ppos = 0;

   *phi_max = pow(10.,-30.);

   for(i = 0; i < ne; i++)
   {

   if(*phi_max < he[i])
   {
   *phi_max = he[i];

   *ppos = i;
   }

   }
return;
}
 



/* function to identify the boundary nodes and then smooth the
   compiled mesh nodes */

void smooth1(void)
{

int i,j,k,test1, nen[10];
   for(i=0; i < nn; i++)
   {
       test1 = -1;
       test1 = check_bound(i);
       if(test1 == -1)
       continue;
       for(j=0; (k=sort[i][j]) != 0; j++)
       nen[j]=k;
       msm(i+1,j,nen);
   }
return;
}


/* function checks the compiled nodes if they are located at
 the boundary by adding the appex angles of all the elements
 common to the ref node and ensuring that their sum < 180 degrees
 to qualify for a boundary node */

int check_bound(int i)
{
int j,k;
double x[3],y[3],angle,L1,L2,DD;

   angle = 0.;
      for(j = 0;(k = sort[i][j]) != 0; j++)
      {
      load_cord(x,y,i,k);
      seglen(x[1],y[1],x[0],y[0],&L1);
      seglen(x[2],y[2],x[0],y[0],&L2);
      dot(x[0],y[0],x[1],y[1],x[0],y[0],x[2],y[2],&DD);
      angle += acos(DD/(L1*L2))*57.29577951;
      }
   angle = fabs(angle - 360);

   if(angle <= 3.)
   return(1);

return(-1);
}
      
/* function to load the cord of the ref node and the remaining nodes
   of each common element to the ref node */

void load_cord(double x[], double y[] ,int i, int k)
{
int l,m,n;
	 l = 0;
	 m = 0;
	 n = 0;
	 do
	 {
	 m = l + 1;
	 if(m > 2) m = 0;
	 n = m + 1;
	 if(n > 2) n = 0;       
	    if(nod[k-1][l] == i + 1)
	    {
	   x[0] = cord[3*(nod[k-1][l] - 1) + 0];
	   y[0] = cord[3*(nod[k-1][l] - 1) + 1];
	   x[1] = cord[3*(nod[k-1][m] - 1) + 0];
	   y[1] = cord[3*(nod[k-1][m] - 1) + 1];
	   x[2] = cord[3*(nod[k-1][n] - 1) + 0];
	   y[2] = cord[3*(nod[k-1][n] - 1) + 1];
	   }
	 l++;
	 }
	 while(nod[k-1][l-1] != i + 1);
return;
}        



/* Function to exchange the diagonal of irregular elements */
void interdiag(void)
{
int i,j,no,n[3],n1[3], test;
double x[3],y[3],L[3],mL;

    for(i = 0; i < ne; i++)
    {
    elemcord(i,n,x,y);
    sidelen(x,y,L);
    
    mL= 0.;
    no= -1;
    maxlen(L,3,&mL,&no);
    lenod(no,n);
  
       for(j=0; j < ne; j++)
      {
      if(j==i)
      continue;
      
      elemcord(j,n1,x,y);
      test = 1000;
      test = comp(n1,n);
	 if(test == 1)
	 {
	swap(n1,n,j,i);
	goto next;
	}
      continue;
      }
   next:
   continue;
   }
return;
}

/*function returns cords of a node */
void nodecord(int node, double *px, double *py)
{
   *px = cord[3*(node-1) + 0];
   *py = cord[3*(node-1) + 1];
return;
}

/* function elemcord returns nodes & cords of elements */
void elemcord(int n,int nn[], double x[], double y[])
{
   nn[0] = nod[n][0];
   nn[1] = nod[n][1];
   nn[2] = nod[n][2];

   nodecord(nn[0],&x[0], &y[0]);
   nodecord(nn[1],&x[1], &y[1]);
   nodecord(nn[2],&x[2], &y[2]);
return;
}

/* function to calc. the side lengths of an element */
void sidelen(double x[], double y[], double L[])
{
int i,j;
   for(i=0; i < 3; i++)
   {
   j=i+1;
   if(j > 2)
   j = 0;
   seglen(x[i],y[i],x[j],y[j],&L[i]);
   }
return;
}

/* function returns the max length or magnitude from an array of 
n, no of lengths  or magnitudes and also indicates its point
 of occurance in the array  */

void maxlen(double L[], int n, double *pL, int *pn)
{
int i;
   for(i=0; i < n; i++)
   {
   if(*pL < L[i])
   {
   *pL = L[i];
   *pn = i;
   }
   }
return;
}

/* function lenod takes the no. which refers to the longest side
of an element and alligns the n[] ie n[0] is the 1st node of the
largest side, n[1] the 2nd node of the same */
void lenod(int no, int n[])
{
int i,j,nb[3];
   j=no;
   for(i=0; i < 3; i++)
   {
   nb[i]=n[j];
   j++;
   if(j > 2)
   j=0;
  }

  for(i=0; i < 3; i++)
  n[i] = nb[i];

return;
}

/* function to compare n[0] ,n[1] with n1[0] and n1[1] and if there
is a match then return 1 else -1 */
int comp(int n1[], int n[])
{
int i,j,k,l, test, tch,wth,nb[3], npar, test1;
double x1,y1,x2,y2,x3,y3,x4,y4;
   j=0;
   for(i=0; i < 3; i++)
   {
   j=i+1;
   if(j > 2)
   j = 0;
   if( (n[0] == n1[i] && n[1] == n1[j]) ||
       (n[1] == n1[i] && n[0] == n1[j]) )
  {
  k=j+1;
  if(k > 2)
  k=0;
  nodecord(n[0],&x1,&y1);
  nodecord(n[1],&x2,&y2);
  nodecord(n[2],&x3,&y3);
  nodecord(n1[k],&x4,&y4);
  test = interseg(x1,y1,x2,y2,x3,y3,x4,y4,&tch,&wth, &npar);
 
 /* allign the n1 array ie common diag nodes are in the
  beginning */
  nb[0]=n1[i];
  nb[1]=n1[j];
  nb[2]=n1[k];
  for(l=0; l < 3; l++)
  n1[l] = nb[l];

    test1 = -1;



  if(test == 1 && test1 == -1)
  return(1);

  return(-1);
  }
  }
return(-1);
}


/* function interchanges the diagonal of the selected elements
by interchanging the last entries of n1[] and n[] and then
feeding n1[] and n[] to nod[j][], nod[i][] */
void swap(int n1[], int n[], int e1, int e)
{
int bin1,bin, i;
double  minA, minA1;

   /* calculate min angle of elems i & j before swap */
   minA = pow(10.,10.);
   minang(n1,n,&minA);

   /* swap diagonals */
   bin1 = n1[2];
   bin = n[2];

   n1[0] = bin;
   n[0] =  bin1;

   /* calculate min angle for the new configuration */
   minA1 = pow(10.,10.);
   minang(n1,n,&minA1);

   if(minA < minA1)   
   for(i=0; i < 3; i++)
   {
   nod[e1][i] = n1[i];
   nod[e][i] = n[i];
   }
return;
}

/* function to calculate minimum value of angle out two vectors */
void minang(int n1[], int n[], double *pmA)
{
int i;
double x1[3],y1[3],x[3],y[3],A[6], min_A;

   /* load the cords */
   for(i=0; i < 3; i++)
   {
   nodecord(n1[i], &x1[i], &y1[i]);
   nodecord(n[i], &x[i], &y[i]);
   }

   calcang(x1,y1,x,y,A);

   /* calculate the min angle from array A */
   min_A = pow(10.,10.);
   for(i=0; i < 6; i++)
   {
   min_A = (min_A > A[i]) ? A[i] : min_A;
   }
   *pmA = min_A;

return;
}



/* function calculates the angles of two arrays n1 & n and 
stores  in A[] */
void calcang(double x1[],double y1[], double x[], double y[], double A[])
{
int i,ii,j,k;
double Aa,Bb,D;
   j = 0;
   k = 0;
   for(i=0; i < 3; i++)
   {

   ii = 2*i;

   j = i+1;
   if(j > 2)
   j=0;

   k = j+1;
   if(k > 2)
   k = 0;

   seglen(x1[i],y1[i],x1[j],y1[j],&Aa);
   seglen(x1[j],y1[j],x1[k],y1[k],&Bb);
   dot(x1[j],y1[j],x1[i],y1[i],
       x1[j],y1[j],x1[k],y1[k], &D);
   A[ii]=acos(D/(Aa*Bb));

   seglen(x[i],y[i],x[j],y[j],&Aa);
   seglen(x[j],y[j],x[k],y[k],&Bb);
   dot(x[j],y[j],x[i],y[i],
       x[j],y[j],x[k],y[k], &D);
   A[ii+1]=acos(D/(Aa*Bb) );
   }
return;
}

/* function to calculate dot product of two vectors */
void dot(double x1, double y1, double x2, double y2,
    double x3, double y3, double x4, double y4, double *pD)
{
   *pD = (x2-x1)*(x4-x3)+(y2-y1)*(y4-y3);
return;
}


/* function to average the neighbouring node cord of the
   internal node */
void msm(int node, int no, int nen[])
{
int i,nrn[30],ii,jj,kk;
   for(i=0; i < no; i++)
   {
   ii = 3*i;
   jj = ii+1;
   kk = jj+1;
   nrn[ii] = nod[nen[i]-1][0];
   nrn[jj] = nod[nen[i]-1][1];
   nrn[kk] = nod[nen[i]-1][2];
   }
   kk++;
   elimcom(node,&kk, nrn);
   avgcom(node,kk,nrn);
return;
}

/*function eliminates the common and the central nodes from the
  array nrn */
void elimcom(int node, int *pkk, int nrn[])
{
int i,j,no,nno, nrnb[30];
   no = *pkk;
   nno = 0;
   for(i=0; i < no; i++)
   {
      for(j=0; j < nno; j++)
      {
      if(nrn[i] == nrnb[j])
      goto skip;
      }
      if(nrn[i] == node)
      continue;
   nrnb[nno] = nrn[i];
   nno++;
   skip:
   continue;
   }

   for(i=0; i < nno; i++)
   nrn[i]=nrnb[i];
   *pkk = nno;

return;
}   
/* the function averages the neighbouring node cord to
   obtain the suitable value for the node under reference */
void avgcom(int node,int kk, int nrn[])
{
int i;
double x,y;
   x = 0;
   y = 0.;
   for(i=0; i < kk; i++)
   {
   x += cord[3*(nrn[i]-1) + 0];
   y += cord[3*(nrn[i]-1) + 1];
   }
   x = x/kk;
   y = y/kk;
    
   cord[3*(node-1) + 0] = x;
   cord[3*(node-1) + 1] = y;
return;
}



/* Function for calculating the intersection of two line segments */
int interseg(double x1,double y1,double x2,double y2,
double x3,double y3,double x4,double y4, int *ptch, int *pwth, 
int *near_parr)
{
double a11=x2-x1;
double a12=x3-x4;
double a21=y2-y1;
double a22=y3-y4;
double c11=x3-x1;
double c21=y3-y1;
double D=a11*a22-a12*a21;
double t,s,tt,ss,tol;
double x,y;

   /* check if the line is near parallel to any segment */
   if(fabs(D) <= .0009)
   *near_parr = 1;
   else
   *near_parr = -1;

   if(fabs(D) > pow(10.,-4.))
   {
   t=(a22*c11-a12*c21)/D;
   s=(a11*c21-a21*c11)/D;
   x=x1+t*(x2-x1);
   y=y1+t*(y2-y1);
   
   tol=pow(10.,-11.);
   tt=fabs(t-1.);
   ss=fabs(s-1.);

   if(t < .999 && t > 0. && s < .999 && s > 0.)
   return(1);

   
   *ptch = -1;
   s=fabs(s);
   
   if(s < tol && t >= 0. && t < .999)
   *ptch = 0;
   
   if(ss < tol && t >= 0. && t < .999)
   *ptch = 1;
   
   return(0);
   }
   else
   {
   *pwth = parallel(x3,y3,x4,y4,x1,y1,x2,y2);
   return(-1);
   }
}




 /* Function parallel check that parallel line (x1,y1)(x2,y2)
    is within line (x3,y3)(x4,y4) and if such is the 
    case then the function returns 1 else -1 */
int parallel(double x1, double y1, double x2, double y2,
	     double x3, double y3, double x4, double y4)
{
double tol, L,  L1, L2, t1, t2, xrhs1, yrhs1, xrhs2, yrhs2;
double xdiff1, ydiff1, xdiff2, ydiff2;
double tot, t, L12;
   seglen(x3,y3,x1,y1, &L1);
   seglen(x3,y3,x2,y2, &L2);
   seglen(x3,y3,x4,y4, &L);
   seglen(x1,y1,x2,y2, &L12);

   t1 = L1/L;
   t2 = L2/L;
   t  = L12/L;

   tol = pow(10.,-4.)*L;
   xrhs1 = x3+t1*(x4-x3);
   yrhs1 = y3+t1*(y4-y3);
   xrhs2 = x3+t2*(x4-x3);
   yrhs2 = y3+t2*(y4-y3);

   xdiff1=fabs(xrhs1-x1);
   ydiff1=fabs(yrhs1-y1);
   xdiff2=fabs(xrhs2-x2);
   ydiff2=fabs(yrhs2-y2);

   tot = xdiff1+ydiff1+xdiff2+ydiff2;

   if(fabs(t1 - 1.) <= tol)
   t1 = 1.;
   if(fabs(t2 - 1.) <= tol)
   t2 = 1.;
   
   if(tot < tol &&  t < .999 && t1 <= 1. && t2 <= 1. )
   return(1);
   else
   return(-1);
}

/* Function to calculate the length of a segment */
void seglen(double x1,double y1, double x2, double y2, double *pL)
{
double L;
L=sqrt(pow(x2-x1,2.)+pow(y2-y1,2.));
*pL = L;
return;
}

/* function stores bg mesh cords and topology. The Function also loads
and stacks the elements which have their hhi_min > hi[]
 */
void store_bg(int *pi)
{
int i,j,vali;
double val, hhi_avg;

       inc_nn = nn_bg = nn;
       inc_ne = ne_bg = ne;
       nb_bg = nb;

       for(i=0; i < nn_bg; i++)
       {
	  for(j=0; j < 2; j++)
	  val = cord_bg[i][j] = cord[3*i + j];
       }

       knt_bg = 0;
       knt_div = 0;
       for(i=0; i < ne_bg; i++)
       {
       vec[i] = 0;
#if SUBDIV 
       vali = ck_small_hhi(i, &hhi_avg);
       
       if(vali == -1){
       subdivide(i,hhi_avg);
       continue;} 
#endif
	  for(j=0; j < 4; j++)
	  nod_bg[i][j] = nod[i][j];

       }


       for(i = 0; i < nb_bg; i++)
	nco_bg[i] = nco[i][0];

   if(iln > 0)
   {
   iln = 0;
more_1:
   iln++;
   i = iln - 1;

   nc_bg[i] = nc[i];

   if(nc[i] < nn)
   goto more_1;
   }

*pi = iln;


       scan_div(); /* mark no-div elements */


       nn_bg = inc_nn;     /* refresh nn_bg and ne_bg */
       ne_bg = inc_ne;


return;
}





/* function to scan over the subdiv bg mesh and mark the elements
   which do not require refinement */

   void scan_div()
   {
   int i, vali;

   for(i = 0; i < ne_bg; i++)
   {
   vali = cmp_side_hhi(i);
   if(vali == -1){ vec[i] = 1;/* if vali= -1 then stack the large element & mark vec=1*/
   p_stack(i);knt_bg++;}
   }

return;
}

/* function scans over each side of the element and finds the min val
   of hhi[] at either of its ends.  Then it checks whether this hh_min
   is larger than the length of the segment  if its true a conter is added
   up, this procedure is repeated for all the 3 sides and then in the
   end if the counter has not been able to reach a val of 3 a 1 val is
   returned else -1 */


int cmp_side_hhi(int i)
{
int j,k;
int nn_n[3], count, vali_1, vali_2;
double  x[3],y[3],L[3],hh_1, hh_2, hhi_min, val;

    elemcord(i, nn_n, x, y);

       k = 0;
       count = 0;
       for(j = 0; j < 3; j++)
       {    
       k = j + 1; if(k > 2) k = 0;
       vali_1 = nod[i][j] - 1;
       vali_2 = nod[i][k] - 1;
       hh_1 = hhi[vali_1];                  /* hhi at 1st end  */
       hh_2 = hhi[vali_2];                  /* hhi at 2nd end  */
       seglen(x[j],y[j],x[k],y[k], &L[j]);
       val = L[j];
       hhi_min = (hh_1 < hh_2) ? hh_1 : hh_2;
       if(hhi_min > .5*val) count++;
       }

      if(count == 3) return(-1);

return(1);

}


/* function to pre-stack the elements which exhibit hhi_min > the
relevent element size hi[] */

void p_stack(int ni)
{
int j,k,l, vali;
double val;

      nn = 3;    ne = 1;

      vali = nod_st[0] =  1;
      vali = nod_st[1] =  2;
      vali = nod_st[2] =  3;
      vali = nod_st[3] = nod[ni][3];

      j = nod[ni][0]; k = nod[ni][1]; l = nod[ni][2];

      val = cord_st[0] = cord_bg[j - 1][0];
      val = cord_st[1] = cord_bg[j - 1][1];

      val = cord_st[3] = cord_bg[k - 1][0];
      val = cord_st[4] = cord_bg[k - 1][1];

      val = cord_st[6] = cord_bg[l - 1][0];
      val = cord_st[7] = cord_bg[l - 1][1];


      stack();

return;
}

/* function to average the hi[] for an element and then average the side
   lengths of the element and then take their ratio.  If the ratio is
   less than .125 say then return -1 else 1.   */

int ck_small_hhi(int i, double *phhi_avg)
{
int j,k;
int nn_n[3];
double  x[3],y[3],L[3],ratio,hhi_avg,L_avg, Min_ang, Fac_1;

    hhi_avg = 0.;                                 /* compute the avg hhi */
    for(j = 0; j < 3; j++)
    hhi_avg += hhi[nod[i][j] - 1];
    hhi_avg = hhi_avg/3.;
    *phhi_avg = hhi_avg;

    elemcord(i, nn_n, x, y);

       k = 0;                                    /* compute L_avg */
       L_avg = 0;
       for(j = 0; j < 3; j++)
       {    
       k = j + 1; if(k > 2) k = 0;
       seglen(x[j],y[j],x[k],y[k], &L[j]);
       L_avg += L[j];
       }
       L_avg = L_avg/3.;

       ratio = hhi_avg/L_avg;

       minang_1(nn_n, &Min_ang);

       Fac_1 = Fac*(Min_ang/30. - 1);

      if(ratio <= Fac_1) return(-1);

return(1);

}

/* function subdivides an elemnent with ck_small_hhi = -1 into 3 elements
   and stores the 1st sub element in the originals place and the remaining
   two at the end of the list */

#if SUBDIV
void subdivide(int i, double hhi_avg)
{
int j,k,l,vali, nn_n[3];
double val, xc,yc,c1,c2,a1,b1,a2,b2, x[3],y[3];

   elemcord(i, nn_n,x,y); /* calc the coords of elem nodes */
   c1 = x[1]*y[2] - x[2]*y[1] - x[2]*y[0] + x[0]*y[2];
   c2 = x[2]*y[0] - x[0]*y[2] - x[0]*y[1] + x[1]*y[0];
   a1 = y[0] + y[1] - 2.*y[2];
   b1 = -x[0] - x[1] + 2.*x[2];
   a2 = -2.*y[0] + y[1] + y[2];
   b2 = 2.*x[0] -x[1] - x[2];
   xc = -c1/a1 - (a1*c2 - a2*c1)/(a2*b1 - a1*b2)*b1/a1;/* centriod of elem */
   yc = (a1*c2 - a2*c1)/(a2*b1 - a1*b2);

   for(j = 0; j < 2 ; j++){  /* load 2  elems at end of nod_bg list */
      k = j + 1; if(k > 2) k = 0;l = j + 2; if(l > 2) l = 0;
      vali = nod_bg[inc_ne + j][0] = nod[i][k];
      vali = nod_bg[inc_ne + j][1] = nod[i][l];
      vali = nod_bg[inc_ne + j][2] = inc_nn + 1;
      vali = nod_bg[inc_ne + j][3] = nod[i][3];
      vec[inc_ne] = 0;}            /* load the vec for transmission enbld*/
      inc_ne += 2;  /* inc the addl 2 elems created at end of list of list */
      knt_div += 2; /* count the elements created */

      vali = nod_bg[i][0] = nod[i][0];
      vali = nod_bg[i][1] = nod[i][1];
      vali = nod_bg[i][2] = inc_nn + 1; /* replace original elem with new topology */

      cord_bg[inc_nn][0] = xc;   /* load the cords for the new node */
      cord_bg[inc_nn][1] = yc;

      val = hhi[inc_nn] = hhi_avg; /* load the hhi val for gen node */
      inc_nn++;



return;
}

#endif

/* function to calculate minimum value of angle out one vectors */
void minang_1(int n1[], double *pmA)
{
int i;
double x1[3],y1[3],A[3], min_A,val;

   /* load the cords */
   for(i=0; i < 3; i++)
   {
   nodecord(n1[i], &x1[i], &y1[i]);
   }

   calcang_1(x1,y1,A);

   /* calculate the min angle from array A */
   min_A = pow(10.,10.);
   for(i=0; i < 3; i++)
   {
   val = A[i];
   min_A = (min_A > val) ? val : min_A;
   }
   *pmA = min_A;

return;
}



/* function calculates the angles of two arrays n1 & n and 
stores  in A[] */
void calcang_1(double x1[],double y1[], double A[])
{
int i,j,k;
double Aa,Bb,D,val;
   j = 0;
   k = 0;
   for(i=0; i < 3; i++)
   {


   j = i+1;
   if(j > 2)
   j=0;

   k = j+1;
   if(k > 2)
   k = 0;

   seglen(x1[i],y1[i],x1[j],y1[j],&Aa);
   seglen(x1[j],y1[j],x1[k],y1[k],&Bb);
   dot(x1[j],y1[j],x1[i],y1[i],
       x1[j],y1[j],x1[k],y1[k], &D);
   val = A[i]=acos(D/(Aa*Bb))*57.29577951;

   }
return;
}


void r_node_mesh(void)
{
register int i;
int dum;

   for(i = 0; i < nn; i++)
   {
   chkf(fgets(buf,BUFF,fp5));
   if(buf[0] == '!') continue;
   sscanf(buf,"%d  %lf",&dum, &hhi[i]);
   printf("%3d    %f\n",i+1, hhi[i]);
   }
}



void node_para(void)
{
int i,j,k;

  
   for(i = 0; i < ne; i++)      /* read in *.me */
   {chkf(fgets(buf,BUFF,fp5)); sscanf(buf,"%d %lf", &j, &hi[i]);}

   search_node();

   for(i=0; i < nn; hhi[i++]=0.);

   /* avg the hi values w.r.t. the connecting elements */
   for(i = 0; i < nn; i++)
   {
       for(j = 0; (k = sort[i][j]) != 0 ; j++)
       hhi[i] += hi[k-1];

       if(j == 0)
       continue;

       hhi[i] = hhi[i]/(double) j;

   }


  }


/* thread to send the information to worker task */

void  bound_send(void)
{

   int i,j,k;
   int ns1,ns2;


   

   for(i = 0; i < ne_bg; i++)
   {
   if(vec[i] == 1) continue;
   if( nn == II && i < ne_bg) exit(1);
     k = 0;
      for(j = 0; j < 3; j++)
      {
      k = j + 1;      if(k > 2)      k = 0;
      ns1 = nod_bg[i][j]; ns2 = nod_bg[i][k];
      ab.hhis[j] = hhi[ns1 - 1];               
      ab.cords[j][0] = cord_bg[ns1 - 1][0];
      ab.cords[j][1] = cord_bg[ns1 - 1][1];
      }

      ab.mats = nod_bg[i][3];
      /*ab.dels = hi[i];*/

      ab.id_no = i + 1; 

      aow(i);     /* call aow which will stack the elements as they are generated*/

    }

}






/* function to load the information on the ratio(hhi_avg/L_avg), ne generated
   per processor and the processor time */

void load_info(int *pi)
{
int i, vali;
double val;
	  i = *pi;
	  vali = ne_ii[i] = cd.ne_i;
	  vali = id[i] = cd.id_no;
	  val = ratio_ii[i] = cd.ratio_i;
	  i++;
	  *pi = i;
return;
}


void magnify()
{
register int i, j;

for(i=0; i < nn; i++)
  for(j = 0; j < 2; j++);
  cord[3*i+j] = cord[3*i+j]*MAG; 
return;
}

void demagnify()
{
register int i, j;

for(i=0; i < nn; i++)
  for(j = 0; j < 2; j++);
  cord[3*i+j] = cord[3*i+j]/MAG; 
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
