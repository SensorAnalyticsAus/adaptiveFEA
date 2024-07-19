#define II 50000

#define JJ 70000

#define LL 50

#define MM 25

#define NN 30

#define LIL 10

#define BUFF 256

#define MAG 100.

#define LINEBUFF 80

typedef struct{int mats, id_no; double dels,cords[3][3],hhis[3];} IN;


typedef struct{int sign,mats,nbs,ilns,ne_i,id_no;
	       double cords[3][3], ratio_i;} OUT;

 int *ne_ii,*id, id_num;


/*  int signal; */

 time_t tima,timb,tim1,tim2;


FILE *fp0,*fp1,*fp5, *fp6, *fp7;

 char buf[BUFF];

 char nout[NN];

 char nout1[NN];

 char titel[32];

 int nn1, ne1, nb1,nm1;  /* assembled values */

 int nn,ne,nb,nm,nl,ndime,nstr,nbe,kmax,nmax,iln,nq;

 int stat, ready;

 int vec[JJ], knt_bg, knt_div;

 int name, inc_nn, inc_ne, ck_fail;

 int nod[JJ][4],nodl[LIL][3],nodq[LIL][5],
nco[MM][4],nsnod[LIL][50],sort[II][20],
mat[20],nc[30],ngs[LIL];

 int *nod_st;

 int nn_bg, ne_bg, nb_bg;


 int bdn1;

 int bdn2;

 int nsgt, en;

 int cknode = 0;

 int ckn =0;

 int n11, n12;

 int n22 ;

 int n21 ;

 int dbg,nodeb;

 int nod_bg[JJ][4], nco_bg[MM],nc_bg[30], n_iln;

 double cord_bg[II][3];
 
 double dt,damp;

 double cord[3*II];

 double *cord_st, *hhi;

 double amat[20][10],aload[30][3];

 double *hi, hmin, hmax, criter;

 double nita_bar,Fac;

 double *ratio_ii;
