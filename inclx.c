#define XX extern

#define II 15000

#define JJ 20000

#define  IIX 3000 

#define  JJX 5000 

#define MM 25

#define NN 30

#define LIL 10

#define BUFF 256

#define LINEBUFF 80

#define NSGT 4000

#define LL 1500 

#define KK 10    /* this is for bc[] which is presently not being used */ 

typedef struct{int mats, id_no; double dels,cords[3][3],hhis[3];} IN;


typedef struct{int sign,mats,nbs,ilns,ne_i,id_no;
	       double cords[3][3], ratio_i;} OUT;

XX FILE *fp7;

XX int *ne_ii,*id, id_num;


/* XX int signal; */

XX char buf[BUFF];

XX char nout[NN];

XX char nout1[NN];

XX char titel[32];

XX int nn1, ne1, nb1,nm1;  /* assembled values */

XX int nn,ne,nb,nm,nl,ndime,nstr,nbe,kmax,nmax,iln,nq;

XX int stat, ready;

/*XX int vec[JJ], knt_bg, knt_div;*/

XX int name, inc_nn, inc_ne, ck_fail;

XX int nod[JJ][4],nodl[LIL][3],nodq[LIL][5],
nco[MM][4],nsnod[LIL][50],sort[II][20],
mat[20],nc[30],ngs[LIL];

XX int *nod_st;

XX int nn_bg, ne_bg, nb_bg;


XX int bdn1;

XX int bdn2;

XX int nsgt, en;

XX int cknode;

XX int ckn;

XX int n11, n12;

XX int n22 ;

XX int n21 ;

XX int dbg,nodeb;

XX double dt,damp;

XX double cord[3*II];

XX double *cord_st/*, *hhi*/;

XX double amat[20][10],aload[30][3];

XX double hi[JJ], hmin, hmax, criter;

XX double nita_bar,Fac;

XX double *ratio_ii;

int bdn, bsgt;

int bseg[NSGT][2],seg[NSGT][2];

double dels, ratio;

double cordx[3*IIX],  dum[LL][2],  cordgen[LL][2],  Ss[LL];

int nodx[JJX][4]; 

double hhix[3];
