/*This a adaptive module for constant strain triangular elements with
  an adaptivity module */


  /* Master main() */

#include <stdio.h>

#include <stdlib.h>

#include <string.h>

#include <math.h>

#include <ctype.h>

#define MM 25

#define BUFF 256

FILE *fp2,*fp3,*fp4,*fp5,*fp6;

char buf[BUFF];

char buf_s[BUFF];

int nn,ne,nb,nm,nl,ndime,nstr,nbe,kmax,nmax,iln,nq;

int *nod, nco[MM][4],*sort, mat[20],nc[30];

double dt,damp;

double *cord;

double *cord_st, *hhi, *steAN, *avg_ste, *avg_err, *err_norm,
	      *disp_norm, *ste;

double amat[5][10],aload[10][3];

double *hi, hmin, hmax; 
int hmin_def;




void rinput(void);
void openfile_fp2(void);
void openfile2(void);
void elem_size_calc(void);
void elemcord(int n,int nn[], double x[], double y[]);
void nodecord(int node, double *px, double *py);
void seglen(double x1,double y1, double x2, double y2, double *pL);
void seglen(double x1,double y1, double x2, double y2, double *pL);
void opena_fp6(void);
void r_stress(void);
void search_node(void);
void avg(void);
int err(void);
void estm(void);
void pr_avg_ste(void);
void pr_avg_err(void);
void pr_D(double d[][3]);
void energy_disp__norm(void);
void energy_norm_g(double *pe_norm );
void u_norm_g(double *pu_norm);
void err_alloc_glob(double energy_g,double disp_g);
void min_he(double he[], double *phi_min, int *ppos);
void max_he(double he[],double *phi_max, int *ppos);
void pr_sort(void);
int inv(int n, double A[][3]);
void elem_area(int el_no,double *pA);
void form_D(double v, double e, double d[][3]);
char* trimwhitespace(char *str);
void chkf(char *str);

void main(int argc, char **argv)
{

hmin_def = 0;
if( argc == 3){
hmin_def = 1; /* TRUE */
strcpy(buf_s,argv[1]);
sscanf(argv[2],"%lf", &hmin);
}

/*-------------------- Read in the dat file ----------------------*/

rinput();

  if((sort = (int *) calloc(nn*20, sizeof(int))) == NULL)
{printf(" out of core for sort\n");exit(1);}



 /* Opening the data output file  for FEM and ADAPTIVITY*/
 openfile_fp2();


  if((ste = (double *) calloc(ne*3, sizeof(double))) == NULL)
{printf(" out of core for ste\n");exit(1);}

/* ----------------Error estimation and adaptivity----------------- */

if((steAN = (double *) calloc(nn*3, sizeof(double))) == NULL)
{printf(" out of core for steAN\n");exit(1);}

if((avg_ste = (double *) calloc(ne*3, sizeof(double))) == NULL)
{printf(" out of core for avg_ste\n");exit(1);}

if((avg_err = (double *) calloc(ne*3, sizeof(double))) == NULL)
{printf(" out of core for avg_err\n");exit(1);}

if((err_norm = (double *) calloc(ne, sizeof(double))) == NULL)
{printf(" out of core for err_norm\n");exit(1);}

if((disp_norm = (double *) calloc(ne, sizeof(double))) == NULL)
{printf(" out of core for disp_norm\n");exit(1);}

if((hi = (double *) calloc(ne, sizeof(double))) == NULL)
{printf(" out of core for hi\n");exit(1);}


elem_size_calc();
openfile2();
opena_fp6();

r_stress();
search_node();
err();

printf("*** Element errors are stored in %s.err ***\n", buf_s);
printf("*** Mesh parameters are in %s.me \n", buf_s);
printf("*** Adaptivity results are in %s.adp\n", buf_s);


fclose(fp2);
fclose(fp3);
fclose(fp4);
fclose(fp6);

}



/*                           AO2.C                           */



/* Function openfile */
void openfile_fp2(void)
{
strcpy(buf,buf_s); strcat(buf,".adp");
fp2=fopen(buf,"w");
	if(fp2==(FILE *)NULL)
	{
	printf("cannot open %s\n",buf);
	exit(1);
	}
}

void openfile2(void)
{
strcpy(buf,buf_s); strcat(buf,".str");
fp3=fopen(buf,"r");
	if(fp3==(FILE *)NULL)
	{
	printf("cannot open %s\n",buf);
	exit(1);
	}
printf("opening %s for reading the stresses\n",buf);
strcpy(buf, buf_s); strcat(buf,".err");
fp4=fopen(buf,"w");
	if(fp4==(FILE *)NULL)
	{
	printf("cannot open %s\n", buf);
	exit(1);
	}
}

void opena_fp6(void)
{
strcpy(buf, buf_s); strcat(buf,".me");
fp6=fopen(buf,"w");
	if(fp6==(FILE *)NULL)
	{
	printf("cannot open %s\n",buf);
	exit(1);
	}
}


/* function to form elsticity matrix D */
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
}



/*** Function to identify common nodes in the elements ***/

void search_node(void)
{
   int i,j,k,l;
   for(i=0; i < nn; i++)
      for(j=0; j < 20; j++)
       sort[20*i + j]=0;

/* Begin searching and sorting the elements connecting to a node,
at a time */
   for(i=0; i < nn; i++)
   {
    l=0;
      for(j=0; j < ne; j++)
      {
          for(k=0; k < 3; k++)
          {
	     if((i+1) == nod[4*j + k])
             {
	      sort[20*i + l]= j+1;
              l++;
              break;
	      }
          }
      }
   }              

}

/* fuction to read the opened stress.dat file for stresses */
void r_stress(void)
{
int i;
double smax,smin;
for(i = 0; i < ne; i++)
{
chkf(fgets(buf,256,fp3));
sscanf(buf,"%lf%lf%lf%lf%lf",&ste[3*i + 0],&ste[3*i + 1],&ste[3*i + 2],&smax,&smin);
printf("%12.5e %12.5e %12.5e\n",ste[3*i + 0],ste[3*i + 1],ste[3*i + 2]);
}
}






/* function to read in the data from the primary data file */

void rinput(void)
{
FILE *fp;
int i,n;
double x1,y1,z1;

if(hmin_def == 0){
   printf("Enter the input file (.dat assumed)\n");   
   chkf(fgets(buf_s,BUFF,stdin)); 
   strcpy(buf_s,trimwhitespace(buf_s));
}

strcpy(buf, buf_s); strcat(buf,".dat");
 
  if((fp=fopen(buf,"r")) == (FILE *) NULL)
   {printf("cannot open %s\n",buf); exit(1);}

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


   if((cord = (double *)calloc(3*nn, sizeof(double))) == NULL)
   {printf(" out of core for cord\n");exit(1);}

   if((nod = (int *)calloc(4*ne, sizeof(int))) == NULL)
   {printf(" out of core for nod\n");exit(1);}


   for(i = 0; i < nn; i++)
   {
   chkf(fgets(buf,256,fp));
   sscanf(buf,"%d %lf %lf %lf",&n, &x1, &y1,&z1);
   cord[3*i] = x1;
   cord[3*i + 1] = y1;
   cord[3*i + 2] = z1;
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
   sscanf(buf,"%d %lf %lf %lf",&nc[i],&aload[i][0],&aload[i][1],&aload[i][2]);

   if(nc[i] < nn)
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
}




/*                        AO7.C                            */




/*** Function to estimate error at each node ***/
int err(void)
{

/* call function to calculate the average nodal stresses */
avg();

/* Function to estimate error in each element by 
the mean of the average nodal stresses from the element stress */
estm();
return(1);
}    

/*** Function avg() calculates the average stress on each node ***/

void avg(void)
{
   int i,j,k,l;
   FILE *gnuplot;
   

/* zero steA */
   for(i=0; i< nn; i++)
   {
      for(j=0; j < 3; j++)
       steAN[3*i + j]= 0.;
   }

/* average the nodal stresses w.r.t the connecting elements */
   for(i=0; i < nn; i++)
   {
      for(j=0; (k=sort[20*i + j]) != 0; j++)
      {
         for(l=0; l < 3; l++)
	  steAN[3*i + l] += fabs(ste[3*(k-1) + l]); /* absolute val */
      }
      for(l=0; l < 3; l++)
      steAN[3*i + l]=steAN[3*i + l]/(double) j; /* avg absolute val */
   }

/* output steAN  */
   fprintf(fp2,"\n\nThe average nodal  stresses steAN\n");
   fprintf(fp2,"Node   Avg.Sx       Avg.Sy     Avg.Sxy \n");
   for(i = 0; i < nn; i++)
   {
   fprintf(fp2,"%3d ", i+1);
      for(j=0; j < 3; j++)
      {
       fprintf(fp2,"%12.5e ", steAN[3*i + j]);
      }
   fprintf(fp2,"\n");
   }           

   /* create data file for gnuplot*/
   strcpy(buf, buf_s);
   strcat(buf,".gst");
   if((gnuplot=fopen(buf,"w")) == (FILE *)NULL){
   printf("cannot open %s\n",buf);
   exit(1);
   }

   for(i = 0; i < nn; i++)
   {
   fprintf(gnuplot,"%12.5f %12.5f ", cord[3*i+0], cord[3*i+1]);
      for(j=0; j < 3; j++)
      {
       fprintf(gnuplot,"%12.5e ", steAN[3*i + j]);
      }
   fprintf(gnuplot,"\n");
   }           
   fclose(gnuplot);

}

/* Function to calculate the error in each element estm() */
void estm(void)
{
   int i,j,k;
   double energy_g, disp_g;
   for(i = 0; i < ne; i++)
   {
   for(j= 0; j < 3; j++)
   avg_ste[3*i + j] = 0.;
   }

   for(i=0; i < ne; i++)
   {
      for(j=0; j < 3; j++) /* sx,sy,sxy */
      {
         for(k=0; k < 3; k++) /* n1,n2,n3 */
	 avg_ste[3*i  + j] += steAN[3*(nod[4*i + k]-1) + j];
	 avg_ste[3*i + j] = avg_ste[3*i + j]/3.;

	 avg_err[3*i + j] = fabs(fabs(ste[3*i + j]) - avg_ste[3*i + j]);
      }
   }

   free((char *)steAN);
   free((char *) ste);

/* print the averaged elem stresses from steAN */
pr_avg_ste();
pr_avg_err(); /* | s* - s^ |  */

/* calc error energy norm || e ||**2 and total energy  
norm || u ||**2 (disp. norm) per element */
energy_disp__norm();


/* calc global values of ||e|| and ||u|| */
energy_norm_g(&energy_g);
u_norm_g(&disp_g);
     
/* translate the error to a range of element sizes and generate
   mesh.dat file at stream fp6 */

     err_alloc_glob(energy_g, disp_g);

}

/* function to print the averaged steAN on each elem as avg_ste */
void pr_avg_ste(void)
{
int i;
   fprintf(fp2,"\n \n elem    avg_ste_x       avg_ste_y       avg_ste_xy \n");
   for(i = 0; i < ne; i++)
   fprintf(fp2,"%5d   %12.5e   %12.5e   %12.5e \n",i+1,
   avg_ste[3*i + 0],avg_ste[3*i + 1],avg_ste[3*i + 2]);
}

/* function to print the avg_err[][] on each elem */
void pr_avg_err(void)
{
int i;
   fprintf(fp2,"\n \n elem    avg_err_x       avg_err_y       avg_err_xy \n");
   for(i = 0; i < ne; i++)
   fprintf(fp2,"%5d   %12.5e   %12.5e   %12.5e \n",i+1,
   avg_err[3*i + 0],avg_err[3*i + 1],avg_err[3*i + 2]);
}

/* function to print D matrix */
void pr_D(double d[][3])
{
int i,j;
   fprintf(fp2,"\n \n D Matrix \n");
   for(i = 0; i < 3; i++)
   {
      for(j = 0; j < 3; j++)
      fprintf(fp2,"%12.5e   ",d[i][j]);
   fprintf(fp2,"\n");
   }
}

/* function to calculate the energy norm */
void energy_disp__norm(void)
{
int i,j,k,test;
double d[3][3];
double  C[3];
double  U[3];
double area;
   /* this is valid for single material type otherwise all nod[][3]
   will have to be scanned to find indivivual d matrices for each elem */

   form_D(amat[0][0], amat[0][1], d);
   
   pr_D(d);

   test = inv(3, d);

   if( test == -1)
   pr_D(d);
   else  
   fprintf(fp2," inversion stopped \n");



   /* zero err_norm & disp_norm */
   for(i = 0; i < ne ; i++)
   {
    err_norm[i] = 0.;
    disp_norm[i] = 0.;
   }

   fprintf(fp2,"\n \n elem   ||e||**2       ||u||**2 \n");

   for(i = 0; i < ne ; i++)
   {

   /* form and invert d matrix for each element of the mesh */
   form_D(amat[nod[4*i + 3] - 1][0], amat[nod[4*i + 3] - 1][1], d);
   test = inv(3, d);
   if(test != -1)
   printf("inversion stopped at element no = %d \n",i+1);


      /* zero C U*/
      for(j = 0; j < 3; j++)
      {
      C[j] = 0.;
      U[j] = 0.;
      }

      /* mult row of avg_err  and avg_ste with col of D(inv) */
      for(j = 0; j < 3; j++)
      {
         for(k = 0; k < 3; k++)
	 {
	 C[j] += avg_err[3*i + k]*d[k][j];
	 U[j] += avg_ste[3*i + k]*d[k][j];
	 }
      }

      /* mult row C[j] with avg_err_T */
         for(j = 0; j < 3; j ++)
         {
	 err_norm[i] += C[j]*avg_err[3*i + j];
	 disp_norm[i] += U[j]*avg_ste[3*i + j];
	 }
   elem_area(i, &area);
   err_norm[i] = err_norm[i]*area;
   disp_norm[i] = disp_norm[i]*area;
   fprintf(fp2,"%5d   %12.5e   %12.5e \n",i+1, err_norm[i], disp_norm[i]);
   }

   free((char *)avg_ste);

}    


/* global energy norm */
void energy_norm_g(double *pe_norm )
{
int i;
double val = 0.;

   *pe_norm = 0.;
   for(i = 0; i < ne ; i++)
   val  += err_norm[i];

*pe_norm = val;

}


/* global u_norm */
void u_norm_g(double *pu_norm)
{
int i;
double val = 0.;


   for(i = 0; i < ne ; i++)
   val += disp_norm[i];

*pu_norm = val;

}
   


/* function allocates the elem size hi[] */
void err_alloc_glob(double energy_g,double disp_g)
{
int i,pos;
double em, *he,he_min,he_max, exi;
double nita, nita_i, mag;

if((he = (double *)calloc(ne, sizeof(double))) == NULL)
{printf(" out of core for he\n");exit(1);}


   nita = sqrt(energy_g/(disp_g));

/* print the actual values of ||e||g**2, ||u||g**2, em, nita
   hi_act[] and the applied val hi[] */

   fprintf(fp2,"\n \n ||e||g**2 = %12.5e \n ||u||g**2 = %12.5e\n",
   energy_g, disp_g);

   fprintf(fp2," nita = %f \n", nita);
   printf("\n nita (percent) from the current mesh = %f\n",nita*100.);

   em = .05*sqrt((disp_g)/(double) ne);

   fprintf(fp2," em = %f \n", em); 

   fprintf(fp2,"\n \n elem   ||e||i**2      ||u||i**2          exi      \
      hi_bg          nita_i %% \n");
   
   for(i =0 ; i < ne; i++)
   {
  


     nita_i = (sqrt(err_norm[i])/(em/.05))*100.;

   exi = sqrt(err_norm[i])/em;

   fprintf(fp2,"%5d   %12.5e   %12.5e   %12.5e   %12.5e   %4.1f \n",i+1, 


   err_norm[i],disp_norm[i],exi, hi[i], nita_i);

   fprintf(fp4,"%12.5f \n",nita_i);

   he[i] = hi[i] /exi;
   }

   /* the min val of he[] is found. The ratio of hmin / he_min
   gives the magnification factor  'mag' , this is used to scale
   up the he[] values provided  the the magnified vals do not exceed
   the hmax */




   
   min_he(he, &he_min, &pos);   
   max_he(he, &he_max, &pos);



   printf(" he_min = %f  he_max = %f\n",he_min, he_max);

if(hmin_def == 0) {
  printf("Enter your own hmin \n");
  chkf(fgets(buf,BUFF,stdin));
  sscanf(buf,"%lf", &hmin);
} /* else hmin is defined in main() */

   mag = hmin/he_min;
   printf("Magnification = %5.1f (hmin = %5.1f)\n", mag, hmin);

   fprintf(fp2,"\n he_min = %f el.no = %d\n", he_min, pos);
   fprintf(fp2,"\n user defined hmin = %f\n", hmin);
   fprintf(fp2,"\n magnification factor mag (hmin/he_min) = %f \n",mag);
   fprintf(fp2,"\n elem    hi[]_actual     hi[]  \n");

   for(i = 0; i < ne ; i++)
   {
  
   hi[i] = he[i]*mag;

   fprintf(fp2," %5d   %12.5e   %12.5e \n",i+1,he[i], hi[i]);

   fprintf(fp6," %d   %f \n",i+1, hi[i]);

   }
}       

/* function to find the min of hi[] */
void min_he(double he[], double *phi_min, int *ppos)
{
int i, ps;
double min;

   min = pow(10.,30.);

   for(i = 0; i < ne; i++)
   {

     if(min > he[i])
     {
     min = he[i];

     ps = i + 1;
     }

   }
   *phi_min = min;
   *ppos = ps;
}


/* function to find the max of hi[] */
void max_he(double he[],double *phi_max, int *ppos)
{
int i, ps;
double max;

   max = pow(10.,-30.);

   for(i = 0; i < ne; i++)
   {

     if(max < he[i])
     {
     max = he[i];
     ps = i + 1;
     }

   }
   *phi_max = max;
   *ppos = ps;
}
 

/* function to output the sort matrix */
void pr_sort(void)
{
int i, j;
   fprintf(fp2," Elements relating to each node \n \n");
   fprintf(fp2,"  Node         Elements (Listed row-wise)\n");
   for(i=0; i < nn; i++)
   {
    fprintf(fp2,"%3d   ", i+1);
      for(j=0; j < 20; j++)
      {
       fprintf(fp2," %3d",sort[20*i + j]);
      }     
    fprintf(fp2,"\n");
   }
}


/*function to invert a matrix */
int inv(int n, double A[][3])
{
int i, j, k, jz, jf;
double p,r, B[3][6];
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < 2*n; j++)
      B[i][j] = 0.;
   }


   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      B[i][j] = A[i][j];
   }

   for(i = 0; i < n; i++)
   {
      jz = n;
      jf = 2*n ;

      for(j = jz; j < jf; j++)
      {
      if(i == j - n )
      goto d3;
      B[i][j] = 0.;
      continue;
d3:   B[i][j] = 1.;
      }
  }
  
   for(k = 0; k < n; k++)
   {
   p = B[k][k];

   if(p == 0.)
   return(1);

   for(j = 0; j < jf; j++)
   B[k][j] = B[k][j]/p;

   for(i = 0; i < n; i++)
   {
   if(i == k)
   continue;

   r = B[i][k];

      for(j = 0; j < jf; j++)
      B[i][j] = B[i][j] - r*B[k][j];

   }
   }

   for(i = 0; i < n; i++)
   {
      k = 0;
      for(j = jz; j < jf; j++)
      {
      A[i][k] =  B[i][j];
      k++;
      }
   }

return(-1);
}



/*------------ bin.c ------------*/
/* function calculates the area of the triangular element */
void elem_area(int el_no,double *pA)
{
double x1,y1,x2,y2,x3,y3;
   
   x1 = cord[3*(nod[4*el_no + 0] - 1) + 0];
   y1 = cord[3*(nod[4*el_no + 0] - 1) + 1];
   x2 = cord[3*(nod[4*el_no + 1] - 1) + 0];
   y2 = cord[3*(nod[4*el_no + 1] - 1) + 1];
   x3 = cord[3*(nod[4*el_no + 2] - 1) + 0];
   y3 = cord[3*(nod[4*el_no + 2] - 1) + 1];

   *pA = .5*((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1));
}

/* function to scan over the elements and work out their sizes
   for storatge in hi[] */
void elem_size_calc(void)
{
int i,j,k;
int nn_n[3];
double  x[3],y[3],L[3], L_avg;

    for(i = 0; i < ne; i++)
    {
    elemcord(i, nn_n, x, y);

       k = 0;
       for(j = 0; j < 3; j++)
       {    
       k = j + 1;
       if(k > 2)
       k = 0;
       seglen(x[j],y[j],x[k],y[k], &L[j]);
       }

   L_avg = (L[0] + L[1] + L[2])/3. ;
   hi[i] = L_avg ;

   }
}

/* function elemcord returns nodes & cords of elements */
void elemcord(int n,int nn[], double x[], double y[])
{
   nn[0] = nod[4*n + 0];
   nn[1] = nod[4*n + 1];
   nn[2] = nod[4*n + 2];

   nodecord(nn[0],&x[0], &y[0]);
   nodecord(nn[1],&x[1], &y[1]);
   nodecord(nn[2],&x[2], &y[2]);
}


/*function returns cords of a node */
void nodecord(int node, double *px, double *py)
{
   *px = cord[3*(node-1) + 0];
   *py = cord[3*(node-1) + 1];
}

/* Function to calculate the length of a segment */
void seglen(double x1,double y1, double x2, double y2, double *pL)
{
double L;
L=sqrt(pow(x2-x1,2.)+pow(y2-y1,2.));
*pL = L;
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
