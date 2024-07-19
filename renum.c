#define T800  0 
#define SUN   1

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MM 25
#define LM 5

extern int nb, iln, nco[MM][4];
extern struct {int nc; double aload[3];} NCLOAD[LM];

void alloc_array(char *array_name, char **AR, unsigned mat_size, size_t bytes);
void identify_bc_load(int nnn, int j, int *ncoo, int *ncc);
void load_bc_load(int *ncoo, int *ncc, int nn);

#if T800
int  *jmem,
     *jnt,
     *joint,
     *memjt,
     *newjt,
     *ncoo,
     *ncc;


static double *cordd, x1, y1;
#endif

void /*FUNCTION*/ setup(int nn, int ne, double *cord, int *nod)
{
int  i, i_, idiff, idime, ielem, ielem_, ii, ii_, iii, 
	 iii_, ik, ik_, inode, ipoin, ipoin_, j, j_, jj, jj_, jjt,
	 jnti, jsub, k, k4, k5, max_, 
	 mem1, minmax, minsav, ndiff, nnn;

#if T800
#else
int  *jmem,
     *jnt,
     *joint,
     *memjt,
     *newjt,
     *ncoo,
     *ncc;


static double *cordd, x1, y1;
#endif

    alloc_array("jem",(char **) &jmem, nn, sizeof(int));
    alloc_array("jnt",(char **) &jnt, nn, sizeof(int));
    alloc_array("joint",(char **) &joint, nn, sizeof(int));
    alloc_array("memjt",(char **) &memjt, 14*ne,sizeof(int));
    alloc_array("newjt",(char **) &newjt, nn, sizeof(int));
    alloc_array("cordd",(char **) &cordd, nn*2, sizeof(double));
    alloc_array("ncoo", (char **) &ncoo, MM, sizeof(int));
    alloc_array("ncc", (char **) &ncc, LM, sizeof(int));  


	idiff = 0;
	for( j = 1; j < nn; jmem[j++] = 0)

	for( j = 1; j <= ne; j++ ){
		j_ = j - 1;
		for( i = 1; i <= 3; i++ ){
			i_ = i - 1;
			jnti = nod[i_+ 4*j_];
			if( jnti == 0)
				goto L_60;
			jsub = (jnti - 1)*8;
			for( ii = 1; ii <= 3; ii++ ){
				ii_ = ii - 1;
				if( ii == i )
					goto L_40;
				jjt = nod[ii_ + 4*j_];
				if( jjt == 0)
					goto L_50;
				mem1 = jmem[jnti - 1];
				if( mem1 == 0)
					goto L_30;
				for( iii = 1; iii <= mem1; iii++ ){
					iii_ = iii - 1;
					if( memjt[jsub + iii_] == jjt )
						goto L_40;
					}
L_30:
				jmem[jnti - 1] = jmem[jnti - 1] + 1;
				memjt[jsub + jmem[jnti - 1] - 1] = jjt;
				if( abs( jnti - jjt ) > idiff )
					idiff = abs( jnti - jjt );
L_40:
				;
				}
L_50:
			;
			}
L_60:
		;
		}
	minmax = idiff;
	minsav = idiff;
    printf("MINMAX INITIAL = %d\n", minmax);

	for( ik = 1; ik <= nn; ik++ ){
		ik_ = ik - 1;
		for( j = 1; j <= nn; j++ ){
			j_ = j - 1;
			joint[j_] = 0;
			newjt[j_] = 0;
			}
		max_ = 0;
		i = 1;
		newjt[0] = ik;
		joint[ik_] = 1;
		k = 1;
L_130:
		k4 = jmem[newjt[i - 1] - 1];
		if( k4 == 0 )
			goto L_145;
		jsub = (newjt[i - 1] - 1)*8;
		for( jj = 1; jj <= k4; jj++ ){
			jj_ = jj - 1;
			k5 = memjt[jsub + jj_];
			if( joint[k5 - 1] > 0 )
				goto L_140;
			k = k + 1;
			newjt[k - 1] = k5;
			joint[k5 - 1] = k;
			ndiff = abs( i - k );
			if( ndiff >= minmax )
				goto L_160;
			if( ndiff > max_ )
				max_ = ndiff;
L_140:
			;
			}
		if( k == nn )
			goto L_150;
L_145:
		i = i + 1;
		goto L_130;
L_150:
		minmax = max_;
		for( j = 1; j <= nn; j++ ){
			j_ = j - 1;
			jnt[j_] = joint[j_];
			}
L_160:
		;
		}
    printf(" MINMAX = %d\n", minmax);

	if( minmax == minsav )
		goto L_5666;
	/* *** NEW NUMBERING SYSTEM                                  */
	for( j = 1; j <= nn; j++ ){
		j_ = j - 1;
		nnn = jnt[j_];
		x1 = cord[0 + 2*j_];
		y1 = cord[1 + 2*j_];
		cordd[0 + 2*(nnn - 1)] = x1;
		cordd[1 + 2*(nnn - 1)] = y1;
        identify_bc_load(nnn,j,ncoo,ncc);
		for( i = 1; i <= ne; i++ ){
			i_ = i - 1;
			for( ii = 1; ii <= 3; ii++ ){
				ii_ = ii - 1;
				if( nod[ii_ + 4*i_] == j )
					nod[ii_ +4*i_] = -nnn;
				}
			}
		}
	for( i = 1; i <= ne; i++ ){
		i_ = i - 1;
		for( j = 1; j <= 3; j++ ){
			j_ = j - 1;
			nod[j_ + 4*i_] = abs( nod[j_ + 4*i_] );
			}
		}
	for( i = 1; i <= nn; i++ ){
		i_ = i - 1;
		for( j = 1; j <= 2; j++ ){
			j_ = j - 1;
			cord[j_ + 2*i_] = cordd[j_ + 2*i_];
			}
		}
L_5666:

load_bc_load(ncoo, ncc, nn);


free((char *)jmem);
free((char *)jnt);
free((char *)joint);
free((char *)memjt);
free((char *)newjt);
free((char *)cordd);
	return;
} /* end of function */


void identify_bc_load(int nnn, int j, int *ncoo, int *ncc)
{
register int i;

for(i = 0; i < nb; i++)
if(j == nco[i][0])
ncoo[i] = nnn;

for(i = 0; i < iln; i++)
if(j == NCLOAD[i].nc)
ncc[i] = nnn;
}

void load_bc_load(int *ncoo, int *ncc, int nn)
{
register int i,j;
int kiln;
double sum;

for(i = 0; i < nb; i++)
nco[i][0] = ncoo[i];

kiln = 0;
for(i = 0; i < iln; i++)
{
  if(ncc[i] == nn) {NCLOAD[i].nc = ncc[i]; kiln++; continue;}
  
  sum = 0.;
  for(j = 0; j < 3; j++) sum += NCLOAD[i].aload[j];

  if(fabs(sum) <= pow(1., -10.)) continue;

  NCLOAD[i].nc = ncc[i];

  kiln++;
}

iln = kiln;

NCLOAD[iln].nc = 15000;

for(i = 0; i < 3; i++)
NCLOAD[iln].aload[i] = 0.;

iln++;

}



#if T800

void alloc_array(char *array_name, char **AR, unsigned mat_size, size_t bytes)
{

if((*AR = calloc(mat_size, bytes)) == NULL)
{printf("out of memory for %s\n", array_name); exit(1);}
}

#endif
