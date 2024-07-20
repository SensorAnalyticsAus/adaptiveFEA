
/* aow.c */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "inclx.c"


void boundnode(int i_i, int bc[],double cord_bgx[][3],int nod_bgx[][4], int ei);
int  fprop(double xo1, double yo1, double xo2, double yo2,
int  nod1, int nod2, int bci, double cord_bgx[][3]);
void node_gen(int *pi, double L, double *pL1, 
     double *pcc, double *pss, double *pxo1, double *pyo1, double *pxo2,
     double *pyo2,double Ss[], int *prod, double *pdel, int nod1, int nod2,
     double cordgen[][2]);
int  dist_err(int *pchk, double del, int Ki, int Kf, double ex, 
     double ey, double Lerr, 
     double L1, double Ss[], 
int  bci, double cordgen[][2],double cord_bgx[][3]); 
void reorder(int l, int Ki);
void afront(int nod1, int nod2, int Ki, int Kf, int chk);
void trian(int i_i, double cord_bgx[][3], int nod_bgx[][4], int ei);
int  search_elbg1( double xa, double ya, 
     double *phi, double cord_bgx[][3], int nod_bgx[][4], int ei);
void seldel(double del1, double L, double *pseldel1);
int  selnode(double PER,double cri, int iseg,
	    double del1,int nod1, int nod2,
	    double ccc[][2], int node[], int test_sel);
int  pt_dist(double xa,double ya, double xb,double yb, 
     double ccc[][2], double xp, double yp);
int  selnode_drs( int iseg,int nod1, int nod2,
     double ccc[][2], int node[]);
int  inter(int iseg, double xa, double ya, 
     double xb, double yb, double xn,
     double yn);
int  dirnode(int nod1, int nod2, double x3, double y3);
void common_seg(int iseg,int selnod, int *pseg1, int *pseg2);
void nearnode_ult(int nod1,int nod2, double x, double y, int node[]);
void nearnode(int nod1,int nod2, double x, double y, int node[], 
	 double del, int *ptest);
void seglen(double x1,double y1, double x2, double y2, double *pL);
void perpcord(double x1,double y1,double x2,double y2,double del1,
double ccc[][2]);
void renum(void);
int  interseg(double x1,double y1,double x2,double y2,
     double x3,double y3,double x4,double y4, int *ptch, int *pwth, 
     int *near_parr);
int  parallel(double x1, double y1, double x2, double y2,
	     double x3, double y3, double x4, double y4);
void prune(int selnod,int iseg,int en_i, int en_f);
void addsegment(int iseg, int selnod);
int  common1(void);
int  common2(void);
void areacord(double x1, double y1, double x2, double y2,
     double x3, double y3, double x, double y, double *pl1,
     double *pl2,double *pl3);
void incord(double x1,double y1,double x2,double y2,
     double *px, double *py);
void outcord(double x1,double y1,double x2,double y2,
     double *px, double *py);
void ck_small_hhi(int i, double *pratio);
void elemcord(int n,int nn[], double x[], double y[]);
void nodecord(int node, double *px, double *py);
XX  void stack(void);
void write_info(double cord_bgx[][3], int en);

/* ------------------- worker task -------------------------------- */
/* function to generate the boundary segments and then call on trian
   to fill the boundary with the triangular elements */


XX IN ab;
XX OUT cd;
XX int *nod_st;
XX double *cord_st;


int mats;

void aow(int i_i)
{
int i,j,vali;
int bc[KK];
int nod_bgx[3][4];

double cord_bgx[3][3], val;

/*----------------------------- worker waits here -----------------*/

      vali = dbg = JJ;
      vali = bsgt = 3;   /* load the information received */
      vali = bdn = 3;
      vali = mats = ab.mats;

      j = 0;
      for(i = 0; i < 3; i++)
      {
      j = i + 1; if(j > 2) j = 0;
      vali = nodx[0][i] = nod_bgx[0][i] = i + 1;
      vali = seg[i][0] = i + 1;
      vali = seg[i][1] = j + 1;
      val  = cordx[3*i+0] = cord_bgx[i][0] = ab.cords[i][0];
      val  = cordx[3*i+1] = cord_bgx[i][1] = ab.cords[i][1];
      val = hhix[i] = ab.hhis[i];
      val = dels = ab.dels;
      }
      vali = nod_bgx[0][3] = ab.mats;

      ck_small_hhi(0,&ratio);

      val = cd.ratio_i =  ratio;

      vali = cd.id_no = ab.id_no;

      boundnode(i_i, bc,cord_bgx,nod_bgx,i);
}





void boundnode(int i_i, int bc[],double cord_bgx[][3],int nod_bgx[][4], int ei)
{
int i;

int nod1,nod2;

double xo1,yo1,xo2,yo2;

   /* Generate nodes on the boundary segment */

   bdn1=bdn;
   nsgt = 0;
   for(i=0; i < bsgt; i++)
   {
   nod1= seg[i][0];  /* identify nodes */
   nod2= seg[i][1];
   xo1=cord_bgx[nod1-1][0]; /* assign x, y cords to identified nodes */
   yo1=cord_bgx[nod1-1][1];
   xo2=cord_bgx[nod2-1][0];
   yo2=cord_bgx[nod2-1][1];

/* generate boundary nodes proportionally */

   fprop(xo1,yo1,xo2,yo2,nod1,nod2,bc[i],cord_bgx);
 }
   /* The final number of boundary nodes after generation are
      stored in extrenal variable bdn2 */

   bdn = bdn1;

   /* write the discretised boundary seg on fp6 mesh.dat */


   /* call triangle generation function */

trian(i_i,cord_bgx, nod_bgx, ei); 

return;
}

/* Function to generate nodes on a boundary segment propotional 
   to the back ground mesh values of hi */
/* Function returns the no. of boundry nodes generated */
int fprop(double xo1, double yo1, double xo2, double yo2,
int nod1, int nod2, int bci, double cord_bgx[][3])
{

   int i,l;
   int roder = 1;
   int Ki = bdn1, Kf; /* the initial boundary nodes and final no.*/
   int chk=0;   
   double L, L1, Lerr;
   double ss,cc;
   double ex,ey;
   double del;

   seglen(xo1,yo1,xo2,yo2,&L);


   i=bdn1;


   cc=(xo2-xo1)/L; /* calculate the direction cosines */
   ss=(yo2-yo1)/L; /* do */

   L1=L; /* load initially the actual length of the segment */

   /* generate nodes proportionally along seg */
   node_gen(&i, L, &L1, &cc, &ss, &xo1, &yo1, &xo2, &yo2, Ss, 
   &roder, &del, nod1,nod2,cordgen); 


   /* distribute the final error  ie final val of L1 to 
      the generated nodes */
   Kf=i;
   ex=L1*cc;
   ey=L1*ss;
   Lerr=L-L1;

   l = dist_err(&chk, del, Ki, Kf, ex, ey, Lerr, L1, 
	    Ss, bci, cordgen,cord_bgx);

   if(roder == -1)
   reorder(l, Ki);


    /* assemble the initial front */
    afront(nod1,nod2,Ki,Kf,chk);



return 0;
}

   /* Function node_gen  distributes nodes along seg */
void  node_gen(int *pi, double L, double *pL1, 
   double *pcc, double *pss, double *pxo1, double *pyo1, double *pxo2,
   double *pyo2,double Ss[], int *prod, double *pdel, int nod1, int nod2,
   double cordgen[][2])
   {
   int i = *pi;
   double x, y, sx, sy, del, L1 = *pL1;
   double del1,del2;
   double cc,ss;
   double  Ltest;
   double dummy1,dummy2, dum;

   cc = *pcc;
   ss = *pss;

   del1 = hhix[nod1 - 1];


   del2 = hhix[nod2 - 1];

   if(del2 > del1) /* interchange node1 with node2 */
   {
   dummy1 = *pxo1;
   dummy2 = *pyo1;
   *pxo1 = *pxo2;
   *pyo1 = *pyo2;
   *pxo2 = dummy1;
   *pyo2 = dummy2;

   dum = del1;
   del1 = del2;
   del2 = dum;

   *pcc = (*pxo2 - *pxo1)/L;
   *pss = (*pyo2 - *pyo1)/L;
   cc = *pcc;
   ss = *pss;
   *prod = -1;
   }

   x = *pxo1; /* load 1st and 2nd ends of ref seg */
   y = *pyo1;

   del = del1;

   do  
   {
   seglen(x,y, *pxo2, *pyo2, &Ltest);

   Ss[i] = Ltest*del/(Ltest + del - del2);

   sx=Ss[i]*cc;
   sy=Ss[i]*ss;
   cordgen[i][0]=x+sx;
   cordgen[i][1]=y+sy;
   x=cordgen[i][0];
   y=cordgen[i][1];

   seglen(x,y,*pxo1,*pyo1,&L1);
   L1=(L-L1);
   *pL1 = L1;

   del = (del1*(L-L1) + del2*L1)/L;

   i++; 
   *pi=i;

   }
   while(L1  > del2);

*pdel = del2;

return;
}

   /* Function dist_err distributes the error at dn2 end
   of segment over the inter nodes gen */

int dist_err(int *pchk, double del, int Ki, int Kf, double ex, 
double ey, double Lerr, 
double L1, double Ss[], 
int bci, double cordgen[][2],double cord_bgx[][3]) 
{
int j ,chk = *pchk;
int l, k;
double sumx, sumy;
double Dfx, Dfy, tol;
      if(fabs(del-L1) <= pow(10.,-7.))
      {
      ex = ey = 0.;
      chk = 1;
      }
   l=bdn1;
   sumx=0.;
   sumy=0.;
   for(j=Ki; j < Kf; j++) 
   {
   sumx +=Ss[j]/Lerr*ex;
   sumy +=Ss[j]/Lerr*ey;
      for(k=0 ; k < bdn; k++)
      {
      Dfx = fabs(cordgen[j][0]+sumx-cord_bgx[k][0]);
      Dfy = fabs(cordgen[j][1]+sumy-cord_bgx[k][1]);
      tol=pow(10.,-9.);
      if(Dfx <= tol && Dfy <= tol)
      goto below;
      }
      cordx[3*j+0] =sumx+cordgen[j][0];
      cordx[3*j+1] =sumy+cordgen[j][1];

      /* load nco[] if bci = 1 x & y freedoms supressed
		       bci = 2   x   freedom supressed
		       bci = 3   y   freedom supressed */
      switch(bci)
      {
      case 1:
      nco[nb][0] = j+1;
      nco[nb][1] = 1;
      nco[nb][2] = 1;
      nb++;
      break;
      case 2:
      nco[nb][0] = j+1;
      nco[nb][1] = 1;
      nco[nb][2] = 0;
      nb++;
      break;
      case 3:
      nco[nb][0] = j+1;
      nco[nb][1] = 0;
      nco[nb][2] = 1;
      nb++;
      break;
      }
      l++;
below:
   continue;
   }

   bdn1 = l;
   *pchk = chk;
return(l);
}

/* Function roder reorders the generated nodes on a
   segment, in ascending order wrt the actually defined
   1st node of the segment */
void reorder(int l, int Ki)
{
int i,j;
 j=0;
   
    for(i=l-1; i >= Ki; i--)
    {
    dum[i][0]=cordx[3*(Ki+j)+0];
    dum[i][1]=cordx[3*(Ki+j)+1];
    j++;
    }

    for(i=l-1; i >= Ki; i--)
    {
    cordx[3*i+0]=dum[i][0];
    cordx[3*i+1]=dum[i][1];
    }
 return;
}

/* function to assemble initial front */
void afront(int nod1, int nod2, int Ki, int Kf, int chk)
{
int i,j;
   /* Assemble the Front */
   j=nsgt;
      if(chk != 0)
      Kf++;
   for(i=Ki; i < Kf; i++)
   {

   if(i==Ki && i==Kf-1)
   {
   bseg[j][0]=nod1;
   bseg[j][1]=nod2;
   j++;
   continue;
   }

   if(i == Ki)
   {
   bseg[j][0]=nod1;
   bseg[j][1]=Ki+1;
   j++;
   continue;
   }

   if(i==Kf-1)
   {
   bseg[j][0]=Kf-1;
   bseg[j][1]=nod2;
   j++;
   continue;
   }
  
   bseg[j][0]=i;
   bseg[j][1]=i+1;
   j++;
   }   
   nsgt=j;
   return;
}








/* Function for triangle generation */

void trian(int i_i, double cord_bgx[][3], int nod_bgx[][4], int ei)
{

int i,j, elem_no;
int nod1,nod2;
int node[20];
int iseg;
int selnod;
int test_dis;
int test_sel;
int kt;
int en_i, en_f, vali;
double x1,y1,x2,y2,xm,ym,L;
double ccc[5][2];
double d1, seldel1, val;
double cri;

/* debug vars */
double dx,dy;
double dL1,dL2, PER; 

       ei = 0;

       if(dbg == 0)
       goto debuga; 

       PER = 0.; /* %age tolerance over selection criteria */
      

       en = 0; /* No of elem generated is set = 0 */

       /* triangs will be generated till nsgt > 0 */

      renum(); /* call renum to place the smallest segment in bseg[][]
	       as the first segment in seg[][] */

       en_i=0;
       en_f=0;

       kt=0;

     while(nsgt > 0)
     { 
		   
       if(dbg == 0)
       break;

       en_i=en_f;

       nod1=seg[0][0];
       nod2=seg[0][1];
       iseg = 0;
		 
	x1=cordx[3*(nod1-1)+0];
	y1=cordx[3*(nod1-1)+1];
	x2=cordx[3*(nod2-1)+0];
	y2=cordx[3*(nod2-1)+1];

	seglen(x1,y1,x2,y2,&L);
	xm = (x1+x2)/2.;
	ym = (y1+y2)/2.;

	d1 = dels;

	elem_no = search_elbg1(xm,ym,&d1, cord_bgx, nod_bgx, ei);

	seldel(d1, L, &seldel1);     

	perpcord(x1,y1,x2,y2,seldel1,ccc);

	for(i=0; i < 20; i++)
	node[i] = 0;
	nearnode(nod1,nod2,ccc[0][0], ccc[0][1],node,seldel1, &test_sel);

	test_dis = pt_dist(x1,y1,x2,y2,ccc,ccc[4][0],ccc[4][1]);

	/* test_sel = 1 indicates that no near node has been selected,
	test_dis  = 1 indicates that the node from ccc is too close
	to a segment of the front hence a near node is required relaxing
	the condition of node being within a specified radius */
	if(test_sel == 1 && test_dis == 1)
	{
	for(i=0; i < 20; i++)
	node[i] = 0;
	nearnode_ult(nod1,nod2,ccc[0][0], ccc[0][1],node);
	test_sel = -1;
	}     

	cri = 1.5;
	selnod=selnode(PER,cri,iseg,seldel1,nod1,nod2,ccc,node,test_sel);

	if(selnod == -1)
	{
	for(i=0; i < 20; i++)
	node[i] = 0;
	nearnode_ult(nod1,nod2,ccc[0][0], ccc[0][1],node);
	selnod=selnode_drs(iseg,nod1,nod2,ccc,node);
	}
  
	if(bdn1 >= IIX || en >= JJX){
	printf("\n ran out of specified memory at nn=%d ne=%d increase IIX or/and JJX \n",bdn1,en);
	exit(1);
	}
 
	if(selnod == -1)
	{
	printf("non-selection of apex node for bg elem no.%d curr nn %d curr ne %d\n",i_i+1, bdn1, en);
/*	break; */
        exit(1);
	}  

	nodx[en][0] = nod1;
	nodx[en][1] = nod2;
	nodx[en][2] = selnod;
	nodx[en][3] = mats;

	vali = cd.mats = mats;

	val = cd.cords[0][0] = cordx[3*(nod1 - 1)+0];
	val = cd.cords[0][1] = cordx[3*(nod1 - 1)+1];
	val = cd.cords[1][0] = cordx[3*(nod2 - 1)+0];
	val = cd.cords[1][1] = cordx[3*(nod2 - 1)+1];
	val = cd.cords[2][0] = cordx[3*(selnod - 1)+0];
	val = cd.cords[2][1] = cordx[3*(selnod - 1)+1];

 for(j = 0; j < 3; j++)
      {
      vali = nod_st[j] = j + 1;
      val = cord_st[3*j + 0] = cd.cords[j][0];
      val = cord_st[3*j + 1] = cd.cords[j][1];
      }
      vali = nod_st[3] = cd.mats;

      stack();

/*
      stat += 1;

      id_num = cd.id_no;

      ck_fail = cd.sign;


      if(ck_fail == 1)
      load_info(&knt);             

      signal += cd.sign;

*/
	en++;
	en_f=en;

       
	prune(selnod,iseg,en_i, en_f);

	vali = cd.sign = 0;
	if(nsgt == 0)
	{vali = cd.sign = 1;     /* load in processor's statistics */
	vali = cd.ne_i = en;}

	if(nsgt == 0)
	break;
   
	renum();
   
	/* debug area */
       debuga:
       dx=cordx[3*(selnod-1)+0];
       dy=cordx[3*(selnod-1)+1];
       seglen(x1,y1,dx,dy, &dL1);
       seglen(x2,y2,dx,dy, &dL2);
       kt++;
       if(kt==dbg)
       {
       if(nodeb == 1)
       break;
       }


   }

/*   nn = bdn1;
   ne = en;      */

printf("Coarse BG mesh elem no.%d resulted in --> nn = %d and ne = %d\n",i_i+1, bdn1,en);

write_info(cord_bgx, en);

return;
}      

		       
/* function searches the element/elements of bg mesh, which contain the
  reference coordinates under reference and returns the inter-
polated value of the element size */
int search_elbg1( double xa, double ya, 
double *phi, double cord_bgx[][3], int nod_bgx[][4], int ei)
{
int m,n1[3];
double x[3], y[3], l[3], sum;
   
      for(m=0; m < 3; m++)
      {
       n1[m] = nod_bgx[ei][m];
       x[m] = cord_bgx[n1[m]-1][0];
       y[m] = cord_bgx[n1[m]-1][1];
      }
   
   areacord(x[0],y[0], x[1], y[1], x[2], y[2], xa,ya,
   &l[0], &l[1], &l[2]);
      sum = 0.;
      for(m=0; m < 3; m++) 
      sum += hhix[n1[m]-1]*l[m];
      *phi = sum;

   return(ei + 1);
}


       /* Function to return an appropiate value of del1 */
void seldel(double del1, double L, double *pseldel1)
{

       if(del1 < 0.55*L)
       *pseldel1 = 0.55*L;

       if(0.55*L <= del1 && del1 <= 2.0*L)
       *pseldel1 = del1;

       if(2.0*L < del1)
       *pseldel1 = 2.0*L;
return;
}

/*Function to finally select the node from node[] and C (x,y) */
int selnode(double PER,double cri, int iseg,
	    double del1,int nod1, int nod2,
	    double ccc[][2], int node[], int test_sel)
{
int i;
int selnod;
double xa,ya,xb,yb, L1,L2,xn,yn;
double tol1, tol2;
int test;
int test2;
int result,result2,stat[20];
int selection = -1;
int tr;

   xa=cordx[3*(nod1-1)+0];
   ya=cordx[3*(nod1-1)+1];
   xb=cordx[3*(nod2-1)+0];
   yb=cordx[3*(nod2-1)+1];

   cknode = 0;
   ckn = 0;
   for(i=0; i < 10; i++)
   stat[i] = 0;
   result = -1;

     if(test_sel != 1)
     {
     for(i=0; node[i] != 0 ; i++)
     {
     cknode = 0; /* at bgin of each iter cknode = 0 */
     ckn=0;
     xn=cordx[3*(node[i]-1)+0];
     yn=cordx[3*(node[i]-1)+1];
     seglen(xa,ya,xn,yn,&L1);
     seglen(xb,yb,xn,yn,&L2);
     test=dirnode(nod1,nod2,xn,yn);
     test2=inter(iseg,xa,ya,xb,yb,xn,yn);

     /* stat[] keeps trace of node/nodes passing both the tests */
     if(test == -1 && test2 == -1)
     {
     result = 1;
     stat[i] = 1;
     }

      tol1=(L1-cri*del1)/L1*100.;
      tol2=(L2-cri*del1)/L2*100.;

      if(test == -1 && test2 == -1 && ((L1 < cri*del1 && L2 < cri*del1) 
      || (tol1 < PER && tol2 < PER)))
      {
      selnod=node[i];
      cknode=1;
      selection = 1;
      goto out;
      }
      }
      }

      for(i=0; i < 5; i++)
      {
      xn = ccc[i][0];
      yn = ccc[i][1];
      test2=inter(iseg,xa,ya,xb,yb,xn,yn);
	 if(test2 == -1) 
	 {
	 bdn1++ ;       
	 selnod=bdn1;
	 selection = 1;
	 ckn=1;
	 tr = i;
	 goto out;
	 }
      }  
   
 
   out:
      /* if a selection is made from ccc and gives result2 = 1
      then the inner ccc[i][] shall be selected till i = 4.
      if still result = 1 then a selection will be made from
      node[] passing both the results */

      if(ckn==1 && selection == 1)
      {

       do
       {
       result2 = pt_dist(xa,ya,xb,yb,ccc,xn,yn);   
       if(result2 == 1 && tr < 3)
       tr++;
       xn = ccc[tr][0];
       yn = ccc[tr][1];
       }
       while(tr < 3 && result2 == 1);

       if(result2 == 1 && result == 1)
	 {
	    bdn1--; /* de-select the selnod */
	       for(i=0;node[i] != 0; i++)
	       {
		  if(stat[i] == 1)
		  return(node[i]);
	       }
	 }

      cordx[3*(bdn1-1)+0]=xn;  /* Load the cords of the selnod as a new
			      node has been gen  */
      cordx[3*(bdn1-1)+1]=yn;
      }

if(ckn != 1 && cknode != 1)
return(-1);

return(selnod);
}


/* function to scan over the segments and find the nearest distance
from the selnod from ccc[][] and if this distance is < distance
between the ref segment mid pt. and the last entry in ccc[][] then
the function returns a value of [1] else [-1] */

int pt_dist(double xa,double ya, double xb,double yb, 
double ccc[][2], double xp, double yp)
{
int i;
int res,tch,wth,parr;
double xm,ym;
double x1,y1,x2,y2, s,t, a11,a12,a21,a22,c1,c2,D,x,y;
double Lc, Lmin, L1;

   /* calculate the distance from the mid pt of ref seg to min ccc */
   xm = (xa+xb)/2.;
   ym = (ya+yb)/2.;
   seglen(ccc[4][0], ccc[4][1], xm, ym, &Lc);

   Lmin = pow(10., 10.);
   for(i=1; i < nsgt; i++)/* i = 1 to skip the ref seg */
   {
   x1 = cordx[3*(seg[i][0]-1)+0];
   y1 = cordx[3*(seg[i][0]-1)+1];
   x2 = cordx[3*(seg[i][1]-1)+0];
   y2 = cordx[3*(seg[i][1]-1)+1];

   /* check that the line from mid of seg to (xp,yp) does not
   intersect any segment */
   res = interseg(xm,ym,xp,yp,x1,y1,x2,y2,&tch,&wth,&parr);
   if(res == 1)
   return(1);

   a11 = x2-x1;
   a12 = y2-y1;
   a21 = a12;
   a22 = -a11;
   c1 = xp-x1;
   c2 = yp-y1;
   D = a11*a22-a12*a21;
   t = (a22*c1-a12*c2)/D;
   s = (a11*c2-a21*c1)/D;
   x = x1+t*(x2-x1);
   y = y1+t*(y2-y1);

   if(fabs(t) <= pow(10.,-10.))
   t = 0.;
   if(fabs(t-1) <= pow(10.,-10.))
   t = 1.;
   if((t >= 0. && t <= 1.) && s < 0.)
   {
   seglen(xp,yp,x,y,&L1);
   Lmin = (L1 > Lmin) ? Lmin : L1; 
   }
   }

   if(Lmin < Lc)
   return(1);

return(-1);
}

/*Function to finally select the node from node[] and C (x,y) */
int selnode_drs( int iseg,int nod1, int nod2,
double ccc[][2], int node[])
{
int i;
int selnod;
double xa,ya,xb,yb, L1,L2,xn,yn;
int test;
int test2;
int selection = -1;


   xa=cordx[3*(nod1-1)+0];
   ya=cordx[3*(nod1-1)+1];
   xb=cordx[3*(nod2-1)+0];
   yb=cordx[3*(nod2-1)+1];

   cknode = 0;
   ckn = 0;
     
     for(i=0; node[i] != 0 ; i++)
     {
     cknode = 0; /* at bgin of each iter cknode = 0 */
     ckn=0;
     xn=cordx[3*(node[i]-1)+0];
     yn=cordx[3*(node[i]-1)+1];
     seglen(xa,ya,xn,yn,&L1);
     seglen(xb,yb,xn,yn,&L2);
     test=dirnode(nod1,nod2,xn,yn);
     test2=inter(iseg,xa,ya,xb,yb,xn,yn);

      if(test == -1 && test2 == -1)
      {
      selnod=node[i];
      cknode=1;
      selection = 1;
      goto out;
      }
      }

      for(i=0; i < 5; i++)
      {
      xn = ccc[i][0];
      yn = ccc[i][1];
      test2=inter(iseg,xa,ya,xb,yb,xn,yn);
	 if(test2 == -1) 
	 {
	 bdn1++ ;       
	 selnod=bdn1;
	 selection = 1;
	 ckn=1;
	 goto out;
	 }
      }  
   
 
   out:
      if(ckn==1 && selection == 1)
      {
      cordx[3*(bdn1-1)+0]=xn;  /* Load the cords of the selnod as a new
			      node has been gen  */
      cordx[3*(bdn1-1)+1]=yn;
      }

if(ckn != 1 && cknode != 1)
return(-1);

return(selnod);
}

/* function inter checks if the selected node does not intersects
   with any of the segments of the generated front. For intersection
the front returns 1 else -1 */
int inter(int iseg, double xa, double ya, 
double xb, double yb, double xn,
double yn)
{
int j;
int n1,n2;
int resa,resb;
double x1,y1,x2,y2;
int tcha,tchb, wtha, wthb, npara, nparb,tst;


     for(j=0; j < nsgt; j++)
	 {
	 n1= seg[j][0];
	 n2= seg[j][1];
	 x1=cordx[3*(n1-1)+0];
	 y1=cordx[3*(n1-1)+1];   
	 x2=cordx[3*(n2-1)+0];
	 y2=cordx[3*(n2-1)+1];

	 tst = -1;
	 if((seg[j][0] == seg[iseg][0] && seg[j][1] ==
	    seg[iseg][1]) ||
	    (seg[j][1] == seg[iseg][0] && seg[j][0] == seg[iseg][1]))
	 tst = 1;
	 
	 tcha = tchb = wtha = wthb = -1000; /* any value except -1, 0, 1 */
	 resa=interseg(xa,ya,xn,yn,x1,y1,x2,y2,&tcha,&wtha,&npara);   
	 resb=interseg(xb,yb,xn,yn,x1,y1,x2,y2,&tchb,&wthb,&nparb);
	
	if(iseg != j && tst != 1)
	{ 
	if(resa == 1 || resb == 1 || (tcha == 0 && tchb == 1)
	  || (tcha == 1 && tchb == 0) || wtha == 1 || wthb == 1 
	  || (resa == -1 && resb == -1))
	 return(1);
	}

 
	if((resa == -1 && resb == -1) || (npara == 1 && nparb == 1))
	return(1);
	else
	continue;

	}
     return(-1);
}


/* Function dirnode returns -1 if the nearest node[i] is on the
left of the seg else 1 */

int dirnode(int nod1, int nod2, double x3, double y3)
{
double x1,y1,x2,y2;
double a11,a12,a21,a22,D,c1,c2,t,s;

   x1=cordx[3*(nod1-1)+0];
   y1=cordx[3*(nod1-1)+1];
   x2=cordx[3*(nod2-1)+0];
   y2=cordx[3*(nod2-1)+1];


   a11=x2-x1;
   a12=y2-y1;
   a21=y2-y1;
   a22 = -1.*(x2-x1);

   D=a11*a22-a12*a21;

   if(D==0.0000001)
   return(1);

   c1=x3-x1;
   c2=y3-y1;
   
   t=(a22*c1-a12*c2)/D;
   s=(a11*c2-a21*c1)/D;

   if(s < 0.)
   return(-1);
   else
   return(1);
}

/* Function common_seg checks whether a generated segment is
   already common to a segment of the current front */

void common_seg(int iseg,int selnod, int *pseg1, int *pseg2)
{
int i,j,k;
int n11,n12;
   n11=seg[iseg][0];
   n12=selnod;
   
   n21=selnod;
   n22=seg[iseg][1];

   j=0;
   k=0;
   for(i=0; i < nsgt ; i++)
   {
    if(j == 0)
    {
       if((n11 == seg[i][0] && n12 == seg[i][1]) ||
       (n12 ==seg[i][0]  && n11 == seg[i][1]))
       {
       *pseg1 = i;
       j++;
       }
       else
       *pseg1 = -1;
    }
    if(k == 0)
    {
       if((n21 == seg[i][0] && n22 == seg[i][1]) ||
       (n22 ==seg[i][0]  && n21 == seg[i][1]))
       {
       *pseg2 = i;
       k++;
       }
       else
       *pseg2 = -1;
    }
   }
return;
}






/* Function to find a max of 5 nearest nodes to C except nod1 & 2*/
void nearnode_ult(int nod1,int nod2, double x, double y, int node[])
{
double L1;
double L2,xx2,yy2;
int i,j,l,k;
int nd;
node[0]=10000;
   
      
   k = 0;
   do
   {
   nd = 0;
   L1=pow(10.,33.);
      for(i=0; i < nsgt; i++)
      {

	 /* check if a node is not repeated */
	 for(j=0; j < 2; j++)
	 {
	 if(nd == seg[i][j])
	 goto next;
 
	 nd = seg[i][j];

	 /* check that the selected node is not already selected */
	 for(l=0; l < k; l++)
	 {
	 if(nd == node[l])
	 goto next;
	 }
	 
	 
	 /* check that the node is not the 1st or 2nd node of iseg */
	 if(nd != nod1 && nd != nod2)
	 {
	 xx2=cordx[3*(nd-1)+0];
	 yy2=cordx[3*(nd-1)+1];
	 seglen(x,y,xx2,yy2,&L2);
	    if(L1 > L2)
	    {

	    L1=L2;
	    node[k] = nd;
	    }
	 }
	 }
next:
      continue;
      }

   k++;

   }
   while( k < 10);
return;
}

/* Function to find a max of 5 nearest nodes to C except nod1 & 2
   and also observing the requirement that the nodes should be
   within a rad = del from the 1st pt in ccc[][]. The function
   returns [1] if no selection has been made else [-1] */
void nearnode(int nod1,int nod2, double x, double y, int node[], 
	 double del, int *ptest)
{
double L1;
double L2,xx2,yy2, rad;
int i,j,l,k,kount;
int nd;
node[0]=10000;
   
      
   k = 0;
   kount =0;
   do
   {
   nd = 0;
   L1=pow(10.,33.);
      for(i=0; i < nsgt; i++)
      {

	 /* check if a node is not repeated */
	 for(j=0; j < 2; j++)
	 {
	 if(nd == seg[i][j])
	 goto next;
 
	 nd = seg[i][j];

	 /* check that the selected node is not already selected */
	 for(l=0; l < k; l++)
	 {
	 if(nd == node[l])
	 goto next;
	 }
	 
	 /* check that the node is within rad = del */
	 xx2=cordx[3*(nd-1)+0];
	 yy2=cordx[3*(nd-1)+1];
	 rad = sqrt(pow(xx2-x,2.) + pow(yy2-y,2.));
	 
	 if(rad > del && fabs(rad - del) > pow(10.,-8.))
	 goto next;
	 
	 
	 /* check that the node is not the 1st or 2nd node of iseg */
	 if(nd != nod1 && nd != nod2)
	 {
	 xx2=cordx[3*(nd-1)+0];
	 yy2=cordx[3*(nd-1)+1];
	 seglen(x,y,xx2,yy2,&L2);
	    if(L1 > L2)
	    {

	    L1=L2;
	    node[k] = nd;
	    kount++;
	    }
	 }
	 }
next:
      continue;
      }

   k++;

   }
   while( k < 10);
   
   if(kount == 0)
   *ptest = 1;
   else
   *ptest = -1;
return;
}

/* Function to calculate the length of a segment 
void seglen(double x1,double y1, double x2, double y2, double *pL)
{
double L;
L=sqrt(pow(x2-x1,2.)+pow(y2-y1,2.));
*pL = L;
return;
}

*/


/* Function to calculate the a perpendicular cord at left 
of segment and at a distance delta1 from the segment */
void perpcord(double x1,double y1,double x2,double y2,double del1,
double ccc[][2])
{
int i;
double x3,y3,L,s[5];
double ang;

seglen(x1,y1,x2,y2,&L);

ang = acos(L/(2*del1));

s[0]=sin(ang)*del1/L;
s[1]=s[0]-.2*s[0];
s[2]=s[0]-.4*s[0];
s[3]=s[0]-.6*s[0];
s[4]=s[0]-.8*s[0];

x3=(x1+x2)/2.;
y3=(y1+y2)/2.;

   for(i=0; i < 5; i++)
   {
   ccc[i][0]=x3-s[i]*(y2-y1);
   ccc[i][1]=y3+s[i]*(x2-x1);
   }

return;
}

/*  function calculate the lengths of the front and the minimum L */
void renum(void)
{
int i,j;
int nod1,nod2,sgmin;
double x1,y1,x2,y2,L,Lmin=pow(10.,33.);

   for(i=0; i < nsgt; i++)
   {
   nod1=bseg[i][0];
   nod2=bseg[i][1];
   x1=cordx[3*(nod1-1)+0];
   y1=cordx[3*(nod1-1)+1];
   x2=cordx[3*(nod2-1)+0];
   y2=cordx[3*(nod2-1)+1];
   
   seglen(x1,y1,x2,y2,&L);
      if(fabs(Lmin-L) > pow(10.,-10.))
      {
	 if(Lmin > L)
	 {
	 Lmin=L;
	 sgmin=i;
	 }
      }
   }

   /* re-order seg[][] wrt sgmin */
    j=0;
   for(i=sgmin ; i < nsgt ; i++ )
   { 
    seg[j][0]=bseg[i][0];
    seg[j][1]=bseg[i][1];
    j++;                 /* j val not to be changed till end of
			    next loop */
   }

   for(i=0 ; i < sgmin; i++)
   {
    seg[j][0]=bseg[i][0];
    seg[j][1]=bseg[i][1];
    j++;
   }
return;
}

/* Function for calculating the intersection of two line segments 
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



*/


 /* Function parallel check that parallel line (x1,y1)(x2,y2)
    is within line (x3,y3)(x4,y4) and if such is the 
    case then the function returns 1 else -1 
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
  
*/

/* Function prune generates a perpcord of 1/20 (say) length  at the
middle side of each segment on the left of the same. the end cords
of this perp cord are checked in each element generated within that
front gen, and segments wich have these cords within an element are
deleted in the usual manner */
void prune(int selnod,int iseg,int en_i, int en_f)
{
int i,j,k;
int ns1,ns2;
int ne1,ne2,ne3;
double  xs1,ys1,xs2,ys2,xs,ys;
double xe1,ye1, xe2,ye2, xe3,ye3;   
double l1,l2,l3;

   addsegment(iseg,selnod);

   k=0;
   for(i=0; i < nsgt; i++)
   {
   ns1=seg[i][0];
   ns2=seg[i][1];
   xs1=cordx[3*(ns1-1)+0];
   ys1=cordx[3*(ns1-1)+1];

   xs2=cordx[3*(ns2-1)+0];
   ys2=cordx[3*(ns2-1)+1];

   incord(xs1,ys1,xs2,ys2,&xs,&ys);

      for(j=en_i; j < en_f; j++) /* scan the gen elements */
      {
      ne1=nodx[j][0];
      ne2=nodx[j][1];
      ne3=nodx[j][2];

      xe1=cordx[3*(ne1-1)+0];
      ye1=cordx[3*(ne1-1)+1];

      xe2=cordx[3*(ne2-1)+0];
      ye2=cordx[3*(ne2-1)+1];

      xe3=cordx[3*(ne3-1)+0];
      ye3=cordx[3*(ne3-1)+1];

      areacord(xe1,ye1,xe2,ye2,xe3,ye3, xs,ys,&l1,&l2,&l3);

      if(fabs(l1) <= pow(10.,-10.))
      l1 = 0.;
      if(fabs(l2) <= pow(10.,-10.))
      l2 = 0.;
      if(fabs(l2) <= pow(10.,-10.))
      l2 = 0.;

      if(fabs(l1 - 1.) <= pow(10.,-10.))
      l1 = 1.;
      if(fabs(l2 - 1.) <= pow(10.,-10.))
      l2 = 1.;
      if(fabs(l3 - 1.) <= pow(10.,-10.))
      l3 = 1.;


      if(l1 >= 0. && l2 >= 0. && l3 >= 0. && l1 <= 1. 
      && l2 <= 1. && l3 <= 1.)
      goto skip;
      }
   bseg[k][0]=seg[i][0];
   bseg[k][1]=seg[i][1];
   k++;
   skip:
   continue;
   }
   
   /* load  the pruned values */
   nsgt=k;
   for(i=0; i < nsgt; i++)
   {
   seg[i][0]=bseg[i][0];
   seg[i][1]=bseg[i][1];
   }

return;
}

/* Function addsegment add the true generated segment 
the function common is called to see if any of the gen segs
is common to the front if true 1 is returned otherwise -1 */
void addsegment(int iseg, int selnod)
{
int j;
int test1, test2;

      j=0;
      n11=seg[iseg][0];
      n12=selnod;
      test1=common1();
	       
      n21=selnod;
      n22=seg[iseg][1];
      test2=common2();
      
      if(test1 == 1 && test2 == -1)
      {
      j++;
      nsgt+=j;
      seg[(nsgt-1)][0]=n21;
      seg[(nsgt-1)][1]=n22;
      }

      if(test1 == -1 && test2 == 1)
      {
      j++;
      nsgt+=j;
      seg[(nsgt-1)][0]=n11;
      seg[(nsgt-1)][1]=n12;
      }

      if(test1 == -1 && test2 == -1)
      {
      j++;
      nsgt+=j;
      seg[(nsgt-1)][0]=n11;
      seg[(nsgt-1)][1]=n12;
      j++;
      nsgt++;
      seg[(nsgt-1)][0]=n21;
      seg[(nsgt-1)][1]=n22;
      }
      
      if(test1 == 1  && test2 == 1)
      j=0;
    
      /* nfrgen = j++; */
return;
}

/* Function to check if the gen side1  is common to front */

int common1(void)
{
int i;
   for(i=0; i < nsgt; i++)
   {
   if((n11==seg[i][0] && n12==seg[i][1]) ||
      (n12==seg[i][0] && n11==seg[i][1]))
   return(1);
   }
return(-1);
}

/* Function to check if the gen side2  is common to front */

int common2(void)
{
int i;
   for(i=0; i < nsgt; i++)
   {
   if((n21==seg[i][0] && n22==seg[i][1]) ||
      (n22==seg[i][0] && n21==seg[i][1]))
   return(1);
   }
return(-1);
}
   
/* Function areacord returns area cord of a trian elem */
void areacord(double x1, double y1, double x2, double y2,
double x3, double y3, double x, double y, double *pl1,
double *pl2,double *pl3)
{
double area,a1,a2,a3;

      area=(x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1);

      a1=(x2*y3-x3*y2)+x*(y2-y3)+y*(x3-x2);
      a2=(x3*y1-x1*y3)+x*(y3-y1)+y*(x1-x3);
      a3=(x1*y2-x2*y1)+x*(y1-y2)+y*(x2-x1);

      *pl1=a1/area;
      *pl2=a2/area;
      *pl3=a3/area;
return;
}

/* Function to calculate the a perpendicular cord at left 
of segment and at a distance 0.001*L from the segment */

void incord(double x1,double y1,double x2,double y2,
double *px, double *py)
{
double x3,y3,L,s;
seglen(x1,y1,x2,y2,&L);
s=.001;
x3=(x1+x2)/2.;
y3=(y1+y2)/2.;
*px=x3-s*(y2-y1);
*py=y3+s*(x2-x1);
return;
}
/* Function to calculate the a perpendicular cord at right 
of segment and at a distance 0.001*L from the segment */

void outcord(double x1,double y1,double x2,double y2,
double *px, double *py)
{
double x3,y3,L,s;
seglen(x1,y1,x2,y2,&L);
s=.001;
x3=(x1+x2)/2.;
y3=(y1+y2)/2.;
*px=x3+s*(y2-y1);
*py=y3-s*(x2-x1);
return;
}




void write_info(double cord_bgx[][3], int en)
{
/* this function write a file named filename.inf containing nnet training info*/
/* it writes to stream fp7 opened and finally closed in main() */


     fprintf(fp7,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %3d\n",
     cord_bgx[0][0], cord_bgx[0][1], 
     cord_bgx[1][0], cord_bgx[1][1],
     cord_bgx[2][0], cord_bgx[2][1],
     hhix[0], hhix[1], hhix[2], en);


return;
}
