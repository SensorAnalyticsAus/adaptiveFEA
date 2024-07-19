/* this task requires minimum of 3 processors */
#define T800         0 
#define SUN          1
#define ITcnt        10   /* display iter counter */
#define DEBUGT800    0    /* with debug only 100 node/proc and 50 band */
#define Timer        1
#define Maxproc      20
#define RHSsub       100
#define MaxOk        RHSsub*50
#define TOL          1.e-01

/*
   This Worker task has 3 channels [0]
   channel_in [0] binded with proc no
   channels[1] for communicating with processor above
   channels[2] for communicating with processors below
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#if SUN
#include <memory.h>
#endif
#if T800
#include <alt.h>
#include <chan.h>
#if Timer
#include <timer.h>
#define  CLK    15625
#endif
#endif

#if T800
int iter;
#else
extern int iter;
#endif


extern unsigned nn2, msize, bsize;

int alloc_array(char *array_name, char **AR,
		 unsigned mat_size, size_t bytes);
#if  T800
void initialize(void);
#else
void initialize(unsigned msize, unsigned bsize, float *A, double *Minv,
		double *r, double *s, double *rs, double *p);
#endif
void calc_Minv(unsigned msize, unsigned bsize, float *A, double *Minv);
void calc_s(unsigned msize, double *Minv, double *r, double *s);
void dotp(unsigned msize, double *a, double *b, double *c);
int cg_iterations(unsigned msize, unsigned bsize, float *A,
		   double *p, double *r, double *s, double *Minv,
		   double *x, double rs);
void mat_vec_mult(unsigned msize, unsigned bsize, float *A,
		  double *p, double *Ap);
void mat_vec_mult_modified(unsigned msize, unsigned bsize, float *A,
			   double *x, double *y);
void vec_scale_add_x(unsigned msize, double *x, double a, double *p);
double vec_scale_add_r(unsigned msize, double *x, double a, double *p);
double converge(unsigned msize, double a, double *Ap);
void calc_p_from_s(unsigned msize, double *p, double *s, double c);

#if T800
void receive_PRE(void);
void receive_or_pass_ok(void);
unsigned global_to_local(unsigned rstart, unsigned nrow_g);
void rec_load_vec(void);
void reflect_back_the_data(void);
void clap_for_mult(void);
void collect_p(void);
void scaler_sum(double *scaler);
void send_results_x(void);
void send_iter(int i, double pAp, double conv_test, int end_of_iter_flag);
#if Timer
void send_proc_time(void);
#endif


struct {unsigned nproc, nn2, bsize, row_proc[Maxproc];int opt_proc;} PRE;
struct {unsigned proc_no, nrow, ncol, rstart, rend, bcflag, end_of_mess;
	float ok;} OK;
#if DEBUGT800
struct {unsigned rstart, rend; float A[MaxOk]; double r[RHSsub];} REF;
#endif
/*struct {unsigned destin, end_of_mess;} LOAD;*/

struct {unsigned iter, end_of_iter_flag; double pAp, conv_test;} ITER;

unsigned nn2, msize, bsize, nproc, proc_no, rstart,rend, bcflag;

unsigned rows_per_proc[Maxproc];

int opt_proc, clap_pt;

float *A;

double *r, *clap;

#if Timer
float proc_tim[2*Maxproc];
float t_cal, t_com, t_start, t1, t2;
#endif

CHAN **inp, **outp;

double
       *x,     /* A.x = b*/
       *Minv,  /* Inverse of diag. A */
       *p,       /* msize+bsize*/
       *Ap,    /* msize+bsize*/
       *s,
       rs;

main(argc, argv, envp, in_ports, ins, out_ports, outs)

int argc, ins, outs;
char *argv[], *envp[];
CHAN *in_ports[], *out_ports[];

#else
void  jcclib(unsigned msize, unsigned bsize, float *A, double *r)
#endif

{

register int i;
int allo, dum;
char filename[30], filename_sav[30], buf[90];

#if T800
#else
double
       *x,     /* A.x = b*/
       *Minv,  /* Inverse of diag. A */
       *p,     /* the projection vector */
       *Ap,
       *s,
       rs;
#endif

#if T800
inp = in_ports;
outp = out_ports;
proc_no = (unsigned) inp[0];
#endif

#if T800


       receive_PRE();

       allo = alloc_array("A", (char **) &A, msize*bsize, sizeof(float));
       allo = alloc_array("clap", (char **) &clap, nn2, sizeof(double));
       allo = alloc_array("r", (char **) &r, msize, sizeof(double));
       allo = alloc_array("x", (char **) &x, msize, sizeof(double));
       allo = alloc_array("Minv", (char **) &Minv, msize, sizeof(double));
       allo = alloc_array("p", (char **) &p, msize+bsize, sizeof(double));
       allo = alloc_array("Ap", (char **) &Ap, msize+bsize, sizeof(double));
       allo = alloc_array("s", (char **) &s, msize, sizeof(double));


       receive_or_pass_ok();

       rec_load_vec();

#if DEBUGT800
       reflect_back_the_data();
#endif


#else
       allo = alloc_array("x", (char **) &x, msize, sizeof(double));
       allo = alloc_array("Minv", (char **) &Minv, msize, sizeof(double));
       allo = alloc_array("p", (char **) &p, msize, sizeof(double));
       allo = alloc_array("Ap", (char **) &Ap, msize, sizeof(double));
       allo = alloc_array("s", (char **) &s, msize, sizeof(double));

       if(allo == -1) return;              /* out of memory error */

#endif

#if T800
#if Timer
t_com = 0.;
t_start = ((float)timer_now())/((float)CLK);

       /* determine mid-point of the pipeline */
       clap_pt = nproc/2 ;

       /* initialization procedure */
	for(i = 0; i < msize; x[i++] = 0.);
	for(i = 0; i < msize+bsize+10; p[i++] = 0.);

#endif
#endif

#if T800
	initialize();

#if Timer
	t1 = ((float)timer_now())/((float)CLK);

	collect_p();

	scaler_sum(&rs);

	t2 = ((float)timer_now())/((float)CLK);
	t_com += t2 - t1;
#endif
#else
	initialize(msize, bsize, A, Minv, r, s, &rs, p);
#endif


       /* begin the loop */
       iter = cg_iterations(msize, bsize, A, p, r, s, Minv, x, rs);


#if T800
#if Timer
      t1 = ((float)timer_now())/((float)CLK);

      send_results_x();

      t2 = ((float)timer_now())/((float)CLK);

      t_com += t2 - t1;

      send_proc_time();
#endif

#else
      memcpy(r,x,msize*sizeof(double));
#endif

#if T800
#else
free((char *)p);
free((char *)Ap);
free((char *)s);
free((char *)x);
#endif
return;
}

int alloc_array(char *array_name, char **AR, unsigned mat_size, size_t bytes)
{

if((*AR = calloc(mat_size, bytes)) == NULL)
return(-1);
else
return(1);
}

#if T800
void initialize (void)
#else
void initialize(unsigned msize, unsigned bsize, float *A,
		double *Minv, double *r, double *s, double *rs, double *p)
#endif
{
calc_Minv(msize, bsize, A, Minv);
calc_s(msize, Minv, r, s);

#if T800
dotp(msize, r, s, &rs);
#else
dotp(msize, r, s, rs);
#endif



#if T800
memcpy(&p[0], s, msize*sizeof(double));
#else
memcpy(p, s, msize*sizeof(double));
#endif
}

void calc_Minv(unsigned msize, unsigned bsize, float *A, double *Minv)
{
register int i,j;

for(i = 0; i < msize; i++)
Minv[i] = 1./A[bsize*i+0];

}

void calc_s(unsigned msize, double *Minv, double *r, double *s)
{
register int i;

for(i = 0; i < msize; i++)
s[i] = Minv[i]*r[i];
}

void dotp(unsigned msize, double *a, double *b, double *c)
{
register int i;

   *c = 0.;

for(i = 0; i < msize; i++)
   *c += a[i]*b[i];
}


int cg_iterations(unsigned msize, unsigned bsize, float *A,
		   double *p, double *r, double *s, double *Minv,
		   double *x, double rs)
{
register int i;
int allo;
int end_of_iter_flag, dum;

double
#if T800
#else
       *Ap,
#endif
       pAp,
       a,
       rs_new,
       conv_test,
       c;

#if T800
#else
       allo = alloc_array("Ap", (char **) &Ap, msize, sizeof(double));
#endif


       for(i = 0; i < 2*nn2; i++)
       {

#if T800
	 mat_vec_mult_modified(msize, bsize, A, &p[0], Ap);


#if Timer
	 t1 = ((float)timer_now())/((float)CLK);

	 clap_for_mult();

	 t2 = ((float)timer_now())/((float)CLK);
	 t_com += t2 - t1;
#endif
#else
	 mat_vec_mult_modified(msize, bsize, A, p, Ap);
#endif


#if T800
	 dotp(msize, &p[0], Ap, &pAp);

#if Timer
	 t1 = ((float)timer_now())/((float)CLK);

	 scaler_sum(&pAp);

	 t2 = ((float)timer_now())/((float)CLK);
	 t_com += t2 - t1;
#endif

#else
	 dotp(msize, p, Ap, &pAp);
#endif

	 a = rs/pAp;


#if T800
	 vec_scale_add_x(msize, x, 1.0*a, &p[0]);
#else
	 vec_scale_add_x(msize, x, +1.0*a, p);
#endif

#if T800
       conv_test = vec_scale_add_r(msize, r, -1.0*a, Ap);

#if Timer
       t1 = ((float)timer_now())/((float)CLK);

       scaler_sum(&conv_test);

       t2 = ((float)timer_now())/((float)CLK);
       t_com += t2 - t1;
#endif

       conv_test = sqrt(conv_test);

       if(conv_test < TOL) end_of_iter_flag = 1;
       else
       end_of_iter_flag = 0;

       if(proc_no == 0)
       send_iter(i+1, pAp, conv_test, end_of_iter_flag);


       if(end_of_iter_flag == 1) return;


#else
       if(vec_scale_add_r(msize, r, -1.0*a, Ap) < TOL) return(i+1);
#endif

       calc_s(msize, Minv, r, s);

       dotp(msize, r, s, &rs_new);

#if T800
#if Timer
       t1 = ((float)timer_now())/((float)CLK);

       scaler_sum(&rs_new);

       t2 = ((float)timer_now())/((float)CLK);
       t_com += t2 - t1;
#endif
#endif
       c = rs_new/rs;

       rs = rs_new;

#if T800
       calc_p_from_s(msize, &p[0], s, c);
#else
       calc_p_from_s(msize, p, s, c);
#endif

#if T800
#if Timer
       t1 = ((float)timer_now())/((float)CLK);

       collect_p();

       t2 = ((float)timer_now())/((float)CLK);
       t_com += t2 - t1;
#endif
#else
       if(fmod(i+1,10) == 0) printf("*** cg.iteration.no %d\n",i+1);
#endif
       }

return(i);
}

void mat_vec_mult(unsigned msize, unsigned bsize, float *A,
		  double *p, double *Ap)
{
register int i,j, jj, kk;
int bsize_lr, jstart, jsub, dum, di, dj;
double sum,Band;

bsize_lr = bsize - 1;

for(i = 0; i < msize; i++)
{
   sum = 0.;

   jj = i - bsize_lr; /*see if bw is not exceeded on left side of the diag.*/

   if(jj > 0) jstart = jj;
   else       jstart = 0;

   kk = msize - i;            /* look for bot. rt triangle of zeroes */

   if(kk < bsize) jsub = bsize - kk;
   else           jsub = 0;

   for(j = jstart; j < i + bsize - jsub; j++)
   {
     di = i; dj = j;
     if(di == dj) Band = A[bsize*di];
     else{
     if(di > dj) {dum = di; di = dj; dj = dum;}
     Band =  A[bsize*di+(dj-di)];
     }
   sum += Band*p[j];
   }
   Ap[i] = sum;
}
}


void mat_vec_mult_modified(unsigned msize, unsigned bsize, float *A,
			   double *x, double *y)
{
register int i, j, ii;
register double *ytmp, *xtmp;

  /* zero y[] */
  ytmp = (double *) &A[0];
  for(i = 0; i < msize+bsize; i++, *ytmp++ = 0.);


  /* Do main diagional */

  xtmp = x;
  ytmp = y;

  for( i=0; i<msize; i++, xtmp++, ytmp++)
    *ytmp = A[bsize*i+0] * (*xtmp);

  /* Loop over other diagionals. */
  for( j=1; j < bsize; j++ ) {

    /* Upper diagional. */

    xtmp = x+j;
    ytmp = y;

#if T800
    for( i=rstart, ii=rstart; i<nn2-j ; i++, ii++, ytmp++, xtmp++){

      if(ii == rend) break;
      *ytmp += A[bsize*(i-rstart)+j] * (*xtmp);
    }
#else
    for( i=0; i<msize-j; i++, ytmp++, xtmp++)
      *ytmp += A[bsize*i+j] * (*xtmp);
#endif



    /* Lower diagional. */

    xtmp = x;
    ytmp = y+j;

#if T800
    for( i=rstart, ii=rstart; i<nn2-j; i++, ii++, ytmp++, xtmp++){

      if(ii == rend) break;
      *ytmp += A[bsize*(i-rstart)+j] * (*xtmp);
    }
#else
    for( i=0; i<msize-j; i++, ytmp++, xtmp++)
      *ytmp += A[bsize*i+j] * (*xtmp);
#endif
  }                                     /* End of diagional loop. */

}


void vec_scale_add_x(unsigned msize, double *x, double a, double *p)
{
register int i;

for(i = 0; i < msize; i++)
 x[i] += a*p[i];
}


double vec_scale_add_r(unsigned msize, double *x, double a, double *p)
{
register int i;
double smod = 0.;

for(i = 0; i < msize; i++)
{
 x[i] += a*p[i];
 smod += pow(a*p[i],2.);
}

return(smod);
}


void calc_p_from_s(unsigned msize, double *p, double *s, double c)
{
register int i;
       for(i = 0; i < msize; i++)
       p[i] = s[i] + c*p[i];
}

#if T800

void receive_PRE(void)
{
unsigned msglen;

msglen = sizeof(PRE);

alt_wait(1,inp[1]);
chan_in_message(msglen, &PRE, inp[1]);

nn2 = PRE.nn2 ;
nproc = PRE.nproc;
bsize = PRE.bsize;
opt_proc = PRE.opt_proc;
msize = PRE.row_proc[proc_no];
memcpy(rows_per_proc, PRE.row_proc, sizeof(rows_per_proc));

if(nproc-1 > proc_no) chan_out_message(msglen, &PRE, outp[2]);
}

void receive_or_pass_ok(void)
{
unsigned msglen, i,j;

msglen = sizeof(OK);

for(;;)
{
alt_wait(1,inp[1]);
chan_in_message(msglen, &OK, inp[1]);    /* receive from above */

if(OK.end_of_mess == 1)
{
if(nproc - 1 > proc_no)
chan_out_message(msglen, &OK, outp[2]);  /* send to proc below */
break;
}

if(OK.proc_no == proc_no)
{
rstart = OK.rstart;
rend = OK.rend;
bcflag = OK.bcflag;
i = global_to_local(rstart, OK.nrow);
j = OK.ncol;

if(bcflag == 1)
A[i*bsize+j]  = OK.ok;
else
A[i*bsize+j] += OK.ok;

continue;                                /* wait for next message */
}
if(nproc - 1 > proc_no)
chan_out_message(msglen, &OK, outp[2]);  /* send to proc below */

}


}

unsigned global_to_local(unsigned rstart, unsigned nrow_g)
{
return(nrow_g - rstart);
}

#if DEBUGT800
void reflect_back_the_data(void)
{
register int knt;
unsigned msglen;


msglen = sizeof(REF);

REF.rstart = rstart;
REF.rend = rend;
memcpy(REF.A, A, sizeof(float)*msize*bsize);
memcpy(REF.r, r, sizeof(double)*msize);

chan_out_message(msglen, &REF, outp[1]);

for(knt = 0; knt < nproc - proc_no - 1; knt++)
{
alt_wait(1, inp[2]);
chan_in_message(msglen, &REF, inp[2]);
chan_out_message(msglen, &REF, outp[1]);
}

}

#endif

void rec_load_vec(void)
{
register int i;
int destin, msgin;


for(i = proc_no; i < nproc; i++)
{
alt_wait(1, inp[1]);
chan_in_word(&destin, inp[1]);
chan_in_word(&msgin, inp[1]);
chan_in_message(msgin, (char *) clap, inp[1]);

if(destin == proc_no)
{
memcpy(r, clap, msgin);
continue;
}
chan_out_word(destin, outp[2]);
chan_out_word(msgin, outp[2]);
chan_out_message(msgin, (char *) clap, outp[2]);
}

}



void clap_for_mult(void)  /* called after a mat_vec multiplication */
{
register int i;
register double *claptmp;
unsigned msglen, bsiz, idno;

bsiz = bsize - 1;

msglen = bsiz*sizeof(double);

if(proc_no == 0)
{
chan_out_message(msglen, (char *) &Ap[msize], outp[2]); /* send below */
return;
}


if(proc_no == nproc - 1)
{
alt_wait(1, inp[1]);
chan_in_message(msglen, (char *) clap, inp[1]);


claptmp = Ap;
for(i = 0; i < bsiz; i++, claptmp++) /* add to Ap[] */
*claptmp += clap[i];

return;
}

/* intermediate processors (the default) */
if(opt_proc == 1)
{
  if(proc_no % 2 == 0)
  {
  chan_out_message(msglen, (char *) &Ap[msize], outp[2]);

  alt_wait(1, inp[1]);
  chan_in_message(msglen, (char *) clap, inp[1]);

  claptmp = Ap;
  for(i = 0; i < bsiz; i++, claptmp++) /* add to Ap[] */
  *claptmp += clap[i];
  }
  else
  {
  alt_wait(1, inp[1]);
  chan_in_message(msglen, (char *) clap, inp[1]);

  claptmp = Ap;
  for(i = 0; i < bsiz; i++, claptmp++) /* add to Ap[] */
  *claptmp += clap[i];

  chan_out_message(msglen, (char *) &Ap[msize], outp[2]);
  }
}
else
{
alt_wait(1, inp[1]);
chan_in_message(msglen, (char *) clap, inp[1]);

claptmp = Ap;
for(i = 0; i < bsiz; i++, claptmp++) /* add to Ap[] */
*claptmp += clap[i];

chan_out_message(msglen, (char *) &Ap[msize], outp[2]);
}

return;
}


void collect_p(void)
{
register int i;
int msglen, msgin;


msglen = bsize*sizeof(double);


if(proc_no == nproc - 1)
{
msglen = (bsize < msize) ? bsize : msize;
msglen = msglen*sizeof(double);
chan_out_word(msglen, outp[1]);
chan_out_message(msglen, (char *) p, outp[1]);
return;
}

if(proc_no == 0)
{
alt_wait(1, inp[2]);
chan_in_word(&msgin, inp[2]);
chan_in_message(msgin, (char *) &p[msize], inp[2]);
return;
}

if(opt_proc == 0)
{
alt_wait(1, inp[2]);
chan_in_word(&msgin, inp[2]);
chan_in_message(msgin, (char *) &p[msize], inp[2]);

chan_out_word(msglen, outp[1]);
chan_out_message(msglen, (char *) p, outp[1]);
return;
}
else
{
   if((nproc - proc_no - 1) % 2 == 0)
   {
   chan_out_word(msglen, outp[1]);
   chan_out_message(msglen, (char *) p, outp[1]);

   alt_wait(1, inp[2]);
   chan_in_word(&msgin, inp[2]);
   chan_in_message(msgin, (char *) &p[msize], inp[2]);
   return;
   }
   else
   {
   alt_wait(1, inp[2]);
   chan_in_word(&msgin, inp[2]);
   chan_in_message(msgin, (char *) &p[msize], inp[2]);

   chan_out_word(msglen, outp[1]);
   chan_out_message(msglen, (char *) p, outp[1]);
   return;
   }

}

}



void scaler_sum(double *scaler)  /* called for collection of scaler*/
{
register int i;
int msglen, ch_no;
double rs_receive;

msglen = sizeof(double);

if(proc_no == 0)
{
chan_out_message(msglen, (char *) scaler, outp[2]); /* send below */

alt_wait(1, inp[2]);
chan_in_message(msglen, (char *) &rs_receive, inp[2]);/*receive from below*/
*scaler = rs_receive;
return;
}

if(proc_no == nproc - 1)
{
chan_out_message(msglen, (char *) scaler, outp[1]);

alt_wait(1, inp[1]);
chan_in_message(msglen, (char *) &rs_receive, inp[1]);
*scaler = rs_receive;
return;
}

if((int) proc_no - clap_pt < 0)
{

alt_wait(1, inp[1]);
chan_in_message(msglen, (char *) &rs_receive, inp[1]);

rs_receive += *scaler;

chan_out_message(msglen, (char *) &rs_receive, outp[2]);


alt_wait(1, inp[2]);
chan_in_message(msglen, (char *) &rs_receive, inp[2]);

*scaler = rs_receive;

chan_out_message(msglen, (char *) &rs_receive, outp[1]);

return;
}


if((int) proc_no - clap_pt > 0)
{

alt_wait(1, inp[2]);
chan_in_message(msglen, (char *) &rs_receive, inp[2]);

rs_receive += *scaler;

chan_out_message(msglen, (char *) &rs_receive, outp[1]);


alt_wait(1, inp[1]);
chan_in_message(msglen, (char *) &rs_receive, inp[1]);

*scaler = rs_receive;

chan_out_message(msglen, (char *) &rs_receive, outp[2]);

return;
}



 /* this is the clap point */
 for(i = 0; i < 2; i++)
 {
  ch_no  = alt_wait(2, inp[1],inp[2]) + 1;
  chan_in_message(msglen, (char *) &rs_receive, inp[ch_no]);

  *scaler += rs_receive;
 }

  chan_out_message(msglen, (char *) scaler, outp[2]);
  chan_out_message(msglen, (char *) scaler, outp[1]);

}

void send_results_x(void)
{
register int knt;
int rin, rout, msgin, msgout, msglen;

msglen = sizeof(double)*msize;


if(proc_no == nproc - 1)       /* last processor */
{
rout = rstart;
msgout = msglen;
chan_out_word(rout,   outp[1]);
chan_out_word(msgout, outp[1]);
chan_out_message(msgout, (char *) x, outp[1]);
return;
}


alt_wait(1, inp[2]);
chan_in_word(&rin,  inp[2]);
chan_in_word(&msgin,inp[2]);
chan_in_message(msgin, (char *) &clap[rin], inp[2]);

memcpy(&clap[rstart], x, msglen);
rout   = rstart;
msgout = msgin+msglen;
chan_out_word(rout,   outp[1]);
chan_out_word(msgout, outp[1]);
chan_out_message(msgout, (char *) &clap[rstart], outp[1]);
return;
}




/* only for 1st procesor*/
void send_iter(int i, double pAp, double conv_test, int end_of_iter_flag)
{
register int knt;
unsigned msglen;

if(i % ITcnt == 0 || end_of_iter_flag == 1)
{
ITER.iter = i;
ITER.pAp = pAp;
ITER.conv_test = conv_test;
ITER.end_of_iter_flag = end_of_iter_flag;

msglen = sizeof(ITER);

chan_out_message(msglen, &ITER, outp[1]);
}

}

#if Timer

void send_proc_time(void)
{
register int i;
int msglen;


msglen = 2*nproc*sizeof(float);



if(proc_no != nproc - 1)
{
alt_wait(1, inp[2]);
chan_in_message(msglen, (char *) proc_tim, inp[2]);
t_cal = t2 - t_start - t_com;
/*if(t_cal == 0) t_cal = 1;
if(t_com == 0) t_com = 1;*/
proc_tim[2*proc_no + 0] = t_cal;
proc_tim[2*proc_no + 1] = t_com;
chan_out_message(msglen, (char *) proc_tim, outp[1]);
}
else
{
t_cal = t2 - t_start - t_com;
/*if(t_cal == 0) t_cal = 1;
if(t_com == 0) t_com = 1;*/
proc_tim[2*proc_no + 0] = t_cal;
proc_tim[2*proc_no + 1] = t_com;
chan_out_message(msglen, (char *) proc_tim, outp[1]);
}

}
#endif

#endif
