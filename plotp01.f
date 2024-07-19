c*********** plotscreen program
c      implicit double precision (a-h,o-z)                               
c      logical*1 triang,cable,numyes,autoscale,nuele,
c     .string, err, stres
c     real sfs
c     character  *2 calcomp,laser,screen
c     parameter    (screen='s',laser='l',calcomp='c')
c      parameter    (nin=7, nout=8)
c
      include 'plot.cmn'
      logical*1 views
c
      open(nin,status='old',file='input.dat')
      open(nout,status='unknown',file='plotout.dat')
      kbug=0

      write(6,*)'output device is postscript file pltf*.plt'

      write(6,*)'Brightness controll factor [0.]...[1.]'
      read(5,*) brightness

      write(6,*)'plot membrane elements ? (t/f)' 
      read(5,*)triang

      write(6,*)'plot element errors or Domain dec? (t/f)'
      read(5,*)err
      if(err) then
      write(6,*)'is this for domain decomposition plot ? (t/f)'
      read(5,*) dec
      if(dec) goto 10 
      write(6,*)'if color error plot w.r.t. max error is required [1]'
      read(5,*)norm
      end if

      if(err.eqv..false.) then   
      write(6,*)'plot stress ? (t/f)'
      read(5,*)stres
      endif

10    write(6,*)'plot cable elements ?    (t/f)'
      read(5,*)cable

      write(6,*)'plot g strings ?   (t/f)'
      read(5,*)string

      write(6,*)'number nodes ?           (t/f)'
      read(5,*)numyes

      write(6,*)'number elements ? (t/f)'
      read(5,*)nuele

      write(6,*)'element reduction factor ? (0.0 - 1.0)'
      read(5,*)redfact

      write(6,*)'auto scale ?             (t/f)'
      read(5,*)autoscale

      if (.not.autoscale) then
          write(6,*)'scale factor ?'
          read(5,*)sfs
      endif

      call data                                                         

      write(6,*)'plot single view ? t/f'
      read(5,*)views

      if(views)then 
      write (6,*)'view no   1/2/3/4 ?'
      read(5,*)nomb
      goto(11,12,13,14)nomb
      endif

  11  nd1=1                                                             
      nd2=2                                                             
      call plot(7,nd1,nd2,.false.) 
      if(views.eqv..true.)goto 122

  12  nd1=1                                                             
      nd2=3                                                             
      call plot(8,nd1,nd2,.false.)
      if(views.eqv..true.)goto 122

  13  nd1=2                                                             
      nd2=3                                                             
      call plot(9,nd1,nd2,.false.)
      if(views.eqv..true.)goto 122

  14  call plot(10,1,1,.true.)
      
122   write(6,*) aminerr,amaxerr,aminst,amaxst

      end

         
      subroutine plot(ichan,nd1,nd2,isoyes)
c****  this is to form a graphics subroutine                            

      include 'plot.cmn'


      double precision  linthick
      parameter        (linthick = .2,  pen = 2)

      kbug=0

c******* open channel to plotfile

      call opench (ichan,twdowx,twdowy)

c******* set up graphics coordinate array

      call setga  (isoyes,nd1,nd2)

c******* find max and min coordinate values

      call maxmin (xmin,xmax,ymin,ymax,xrange,yrange,sf,
     .             twdowx,twdowy,cmxoff2,cmyoff2)

c******  calculate the shift in origin and the scalefactor for ps

      call setscale(cmxoff2, cmyoff2)

c***** open graphics window and clear screen

      call openg (linthick,twdowx,twdowy,cmxoff2,sf,
     .            cmyoff2,xmin,ymin,xmax,ymax)

        if(triang) then

        call psplottri(sf,cmxoff2,cmyoff2)
        call psplotquad(sf)

        end if

        if(cable) call pltcbl(sf)

        if(string) call pltgst(sf)

        call  nodnum (kbug,sf,cmxoff2,cmyoff2)
   
        call  elemno(sf,cmxoff2,cmyoff2)

        if(dec.eqv..false.) call indexerrstr

        call closing

      return                                                            
      end


c  ************************************************************************
c  *                   orthogonal viewing transformation                  *
c  ************************************************************************

       subroutine orthog(x,y,z,u1,v1)

       include 'plot.cmn'

       u1=cang(2)*cang(3)*x-cang(2)*sang(3)*y+sang(2)*z
       v1=(sang(1)*sang(3)-cang(1)*sang(2)*cang(3))*x+(sang(1)*cang(3)+
     +   cang(1)*sang(2)*sang(3))*y+cang(1)*cang(2)*z

       return
       end




c*******************************************************************************

                 subroutine opench (ichan,twdowx,twdowy)

c******   open channel to plotfile 

      include      'plot.cmn'

      integer       linscr, linlas
      parameter    (linscr = 5, linlas = 1)

          plotline=linlas
          twdowx=250.00d+00
          twdowy=200.00d+00
        ich=ichan-6

        goto(102,104,106,108) ich
 102      plotfile='pltf1.plt'
          goto 110
 104      plotfile='pltf2.plt'
          goto 110
 106      plotfile='pltf3.plt'
          goto 110
 108      plotfile='pltf4.plt'
 110    continue
     

      return
      end


c*******************************************************************************

                 subroutine setga (isoyes,nd1,nd2)

c******   set up graphics coordinate array 

      include    'plot.cmn'

      write(6,*)' isoyes =  ',isoyes
      if (isoyes) then
          write(6,*)'input vertical angle'
          read(5,*)a1
          write(6,*)'input horizontal angle'
          read(5,*)a3
          a2=0.0d0
          pi=4.0d+00*atan(1.0d+00)
          a1=a1*pi/180.0d0
          a2=a2*pi/180.0d0
          a3=a3*pi/180.0d0
          cang(1)=cos(a1)
          sang(1)=sin(a1)
          cang(2)=cos(a2)
          sang(2)=sin(a2)
          cang(3)=cos(a3)
          sang(3)=sin(a3)
          do 10 i=1,nn
            call orthog(cord(i,1),cord(i,2),cord(i,3),
     +         tcor(i,1),tcor(i,2))
c            if (i.lt.5) then
c            endif
 10       continue
        else
          do 20 i=1,nn
          tcor(i,1)=cord(i,nd1)
          tcor(i,2)=cord(i,nd2)
 20       continue
      endif

      return
      end


c******************************************************************************

          subroutine maxmin (xmin,xmax,ymin,ymax,xrange,yrange,sf,
     .                       twdowx,twdowy,cmxoff2,cmyoff2)

c******     find max and min coordinate values

      include   'plot.cmn'

      xmax=tcor(1,1)                                                  
      xmin=tcor(1,1)                                                  
      ymin=tcor(1,2)                                                  
      ymax=tcor(1,2)                                                  
      do 678 j=1,nn                                                     
      xd=tcor(j,1)                                                    
      yd=tcor(j,2)                                                    
      if(xd.lt.xmin)xmin=xd                                             
      if(xd.gt.xmax)xmax=xd                                             
      if(yd.lt.ymin)ymin=yd                                             
      if(yd.gt.ymax)ymax=yd                                             
 678  continue                                                          
      xrange=xmax-xmin
      yrange=(ymax-ymin) * 1.1
      if (autoscale) then
          if (xrange.eq.0.0d+00) then
              sf1=1.0d+10
          else
              sf1=(twdowx-20.0d+00)/xrange
          endif
          if (yrange.eq.0.0d+00) then
              sf2=1.0d+10
          else
              sf2=(twdowy-20.0d+00)/yrange
          endif
          if (sf1.lt.sf2) then
              sf=sf1
          else
              sf=sf2
          endif
      else
          sf=sfs
      endif
      cmxrng=sf*xrange
      cmyrng=sf*yrange
      cmxmin=sf*xmin
      cmymin=sf*ymin
      cmxoff1=(twdowx-cmxrng)/2.0d+00
      cmyoff1=(twdowy-cmyrng)/3.0d+00
      cmxoff2=cmxoff1-cmxmin
      cmyoff2=cmyoff1-cmymin
      write(nout,67)xmin,xmax,ymin,ymax                                    
 67   format(//,' max min values:',4(3x,f10.4))                         

      return
      end

c*****************************************************************************

                          subroutine setscale(cmxoff2, cmyoff2)
      include 'plot.cmn'

             scalefac = 72./25.4/.24

             shiftx = cmxoff2*scalefac 
             shifty = cmyoff2*scalefac
      return
      
      end


c******************************************************************************

      subroutine  openg (linthick,twdowx,twdowy,cmxoff2,sf,
     .                   cmyoff2,xmin,ymin,xmax,ymax)
                 
c******  open graphics window and clear screen

      include  'plot.cmn'

      double precision   linthick

      if (kbug.gt.2) write(6,5) plotfile
 5    format(1x,a15,'  opened ')
c***** draw border
c----- call sr pboder to set up ps instructions to draw border 
      call pborder(cmxoff2,cmyoff2,sf)
      return
      end
c*****************************************************************************

                 subroutine pborder(cmxoff2,cmyoff2,sf)

       include 'plot.cmn'
       character*36 til
       open(10, status = 'unknown', file = plotfile)
       write(10,1000)
1000   format('%!ps-adobe-2.0 epsf-2.0',/,' %%boundingbox: 0 0 580 
     +  829',/,' /l {lineto} def /m {moveto} 
     +  def 577 14 translate 90 rotate ')
       write(6,*)'enter the title'
	   read(5,1202) til
 1202  format(a36)
	   write(10,1203) til
1203   format(' ',' /helvetica findfont 20 scalefont setfont 370 500
     + m',/,' (',a36,') show 2 setlinewidth stroke')
       write(10,1001)
1001   format(' ','.24 .24 scale newpath')
       return
       end

c******************************************************************************

                 subroutine  nodnum (kbug,sf,cmxoff2,cmyoff2)

c******   number nodes

      include   'plot.cmn'

      if (numyes) then

       ixo = int(shiftx)
       iyo = int(shifty)

       write(10,1002)
1002   format(' ',' /helvetica findfont
     + 10 .24 div scalefont setfont ')

          do 56 jt=1,nn 
c          if(jt.eq.1) then                                                    
c          write(10,900)
c900       format(' ','  (1  ) show')
c          else
          x = sf*tcor(jt,1)*scalefac + shiftx
          y = sf*tcor(jt,2)*scalefac + shifty
          ix = int(x)
          iy = int(y)
          write(10,901)ix,iy,jt
901       format(' ',i5,i5,' m  (',i3,') show ')
c          end if
            if (kbug.gt.2) write(nout,55) jt, (tcor(jt,j), j=1,2)
 55         format('  tcor(',i3,',*) :',2(4x,f12.5))
 56       continue                                                          
      endif

      return
      end

c*******************************************************************

             subroutine elemno(sf,cmxoff2,cmyoff2)

      include 'plot.cmn'
      kbug = 0
c******** mark elements if required
      if(nuele) then

        do 285 j=1,ne
          xdif=0.0d0
          ydif=0.0d0
          do 367 jj=1,3
            xdif=xdif+tcor(nod(j,jj),1)
            ydif=ydif+tcor(nod(j,jj),2)
  367     continue
          xdif= sf*xdif/3.0d0*scalefac + shiftx
          ydif= sf*ydif/3.0d0*scalefac + shifty
          ixdif = int(xdif)
          iydif = int(ydif)

          write(10,1004)ixdif,iydif,j
1004     format(' ',' /courier findfont 100 scalefont setfont ',i5,i5,
     +   ' m  (',i5,')',' show')
          if (kbug.gt.2) write(nout,28) xdif, ydif
  28      format('  xdif, ydif : ',2(4x,f12.5))
  285   continue
      endif

      return

      end

             subroutine  psplottri(sf,cmxoff2,cmyoff2)

c******   plot triangular elements in ps file

      include      'plot.cmn'

      kbug = 0

c---------------------------------- plot elements
      do 164 j=1,ne
        n3=nod(j,3)
        x3=tcor(n3,1)
        y3=tcor(n3,2)
        x = sf*x3*scalefac + shiftx
        y = sf*y3*scalefac + shifty
        ix3 = int(x)
        iy3 = int(y)
        if (redfact.gt.0.0) then
c ---------------------------------- find centre of element
          xdif=0.0d0
          ydif=0.0d0
          do 369 jj=1,3
            xdif=xdif+tcor(nod(j,jj),1)
            ydif=ydif+tcor(nod(j,jj),2)
 369      continue
          xc=xdif/3.0d0
          yc=ydif/3.0d0
c ---------------------------------- modify x3,y3
          x3=x3 + (xc-x3)*redfact
          y3=y3 + (yc-y3)*redfact
          x3 = sf*x3*scalefac + shiftx
          y3 = sf*y3*scalefac + shifty
          ix3 = int(x3)
          iy3 = int(y3)
        endif
        write(10,1005) ix3,iy3
1005    format(' ',i5,i5,' m ')
        do 165 i=1,3
          np=nod(j,i)
          xq=tcor(np,1)
          yq=tcor(np,2)
          xq = sf*xq*scalefac + shiftx
          yq = sf*yq*scalefac + shifty
          ixq = int(xq)
          iyq = int(yq)
          if (redfact.gt.0.0) then
            np = nod(j,i)
            xq = tcor(np,1)
            yq = tcor(np,2)
            xq=xq + (xc-xq)*redfact
            yq=yq + (yc-yq)*redfact
            xq = sf*xq*scalefac + shiftx
            yq = sf*yq*scalefac + shifty
            ixq = int(xq)
            iyq = int(yq)
          endif
          write(10,1006)ixq,iyq
1006      format(' ',i5,i5,' l ')          
165       continue

          if(err.or.stres) then

          if(dec) then
          callshadedom(j)
          goto 164
          endif


          if(norm.eq.1) then
          call shaderr1(j)        
          else
          call shaderr(j)
          endif  

          call shadste(j)
     
          else
          
          write(10,1205)
1205      format(' ','stroke')

          end if
          
164   continue

      return
      end


             subroutine  psplotquad(sf)

c******   plot rect elements in ps file

      include      'plot.cmn'

      kbug = 0

c---------------------------------- plot elements
      do 164 j=1,nq
        n4=nod(j,4)
        x4=tcor(n4,1)
        y4=tcor(n4,2)
        x = sf*x4*scalefac + shiftx
        y = sf*y4*scalefac + shifty
        ix = int(x)
        iy = int(y)
        if (redfact.gt.0.0) then
c ---------------------------------- find centre of element
          xdif=0.0d0
          ydif=0.0d0
          do 369 jj=1,4
            xdif=xdif+tcor(nod(j,jj),1)
            ydif=ydif+tcor(nod(j,jj),2)
 369      continue
          xc=xdif/4.0d0
          yc=ydif/4.0d0
c ---------------------------------- modify x3,y3
          x4=x4 + (xc-x4)*redfact
          y4=y4 + (yc-y4)*redfact
          x4 = sf*x4*scalefac + shiftx
          y4 = sf*y4*scalefac + shifty
          ix = int(x4)
          iy = int(y4)
        endif
        write(10,1005) ix,iy
1005    format(' ',i5,i5,' m ')
        do 165 i=1,4
          np=nod(j,i)
          xq=tcor(np,1)
          yq=tcor(np,2)
          xq = sf*xq*scalefac + shiftx
          yq = sf*yq*scalefac + shifty
          ixq = int(xq)
          iyq = int(yq)
          if (redfact.gt.0.0) then
            np = nod(j,i)
            xq = tcor(np,1)
            yq = tcor(np,2)
            xq=xq + (xc-xq)*redfact
            yq=yq + (yc-yq)*redfact
            xq = sf*xq*scalefac + shiftx
            yq = sf*yq*scalefac + shifty
            ixq = int(xq)
            iyq = int(yq)
          endif
          write(10,1006)ixq,iyq
1006      format(' ',i5,i5,' l ')          
165       continue
          

          
          if(err.or.stres) then
      
          if(norm.eq.1) then
          call shaderr1(j)
          else
          call shaderr(j)
          end if

          call shadste(j)
         
          else
          write(10,1204)
 1204     format(' ','stroke')
          end if          
        

164   continue

      
      if(norm.eq.1) then
c      
c      write(10,1850) amaxer
c1850  format(' ','/times-roman findfont 50 scalefont setfont',/,
c     .' 1200 0 m (relative error plot max error=',f5.1,'%) show')
      
      else
      if(err) then
      if(dec) goto 200
      write(10,1802)
1802  format(' ','/times-roman findfont 50 scalefont setfont',/,
     .' 1200 0 m (error plot) show')
      end if
      end if

c      if(stres) then
c         if(ipv.eq.1) chstr = 'xx'
c         if(ipv.eq.2) chstr = 'yy'
c         if(ipv.eq.3) chstr = 'xy'
c         if(ipv.eq.4) chstr = 'max-p'
c         if(ipv.eq.5) chstr = 'min-p'
c	  write(10,1800) chstr,amaxst
c1800  format(' ','0 10 m /times-roman findfont 50 
c     .scalefont setfont',/,' (100% stress ',a5,' = ',e12.5,') show')
c      write(10,1801)
c1801  format(' ','1200 10 m (stress plot) show')
c      end if


200   return
      end

c******    shade the elements with respect to specified parameter

      subroutine shaderr(j)

      include 'plot.cmn'

      if(err) then
       if(dec) errorr(j) = 10.*errorr(j)
         if(abs(errorr(j)).ge. 1.) shade = 0. 
         if(abs(errorr(j)).lt. 1.) shade = .1
         if(abs(errorr(j)).lt. 0.9) shade = .2
         if(abs(errorr(j)).lt. 0.8) shade = .3
         if(abs(errorr(j)).lt. 0.7) shade = .4
         if(abs(errorr(j)).lt. 0.6) shade = .5
         if(abs(errorr(j)).lt. 0.5) shade = .6
         if(abs(errorr(j)).lt. 0.4) shade = .7
         if(abs(errorr(j)).lt. 0.3) shade = .8
         if(abs(errorr(j)).lt. 0.2) shade = .9
         if(abs(errorr(j)).lt. 0.1) shade = 1.
          
      write(10,1008) shade
1008  format(' ','closepath gsave ',f4.1,' setgray fill
     + grestore stroke ')
      write(10,1009)
1009  format(' ','0.0 setgray')
      end if

      return

      end
c******    shade the elements with respect to specified parameter

      subroutine shaderr1(j)

      include 'plot.cmn'

      
c         if(abs(errorr(j)/amaxer).ge. 1.) shade = 0. 
c         if(abs(errorr(j)/amaxer).lt. 1.) shade = .1
c         if(abs(errorr(j)/amaxer).lt. 0.9) shade = .2
c         if(abs(errorr(j)/amaxer).lt. 0.8) shade = .3
c         if(abs(errorr(j)/amaxer).lt. 0.7) shade = .4
c         if(abs(errorr(j)/amaxer).lt. 0.6) shade = .5
c         if(abs(errorr(j)/amaxer).lt. 0.5) shade = .6
c         if(abs(errorr(j)/amaxer).lt. 0.4) shade = .7
c         if(abs(errorr(j)/amaxer).lt. 0.3) shade = .8
c         if(abs(errorr(j)/amaxer).lt. 0.2) shade = .9
c         if(abs(errorr(j)/amaxer).lt. 0.1) shade = 1.



      call getrgberr(errorr(j))


c      write(10,1008) shade
c1008  format(' ','closepath gsave ',f4.1,' setgray fill
c     + grestore stroke ')

      write(10,1008) red, green, blue
1008  format(' ','closepath gsave ',3(1X,f4.1),' setrgbcolor fill
     + grestore stroke ')

      write(10,1009)
1009  format(' ','0.0 setgray')
      

      return

      end



      subroutine getrgberr(evalue)

      include 'plot.cmn'

      amaxstr = amaxer

      btf = 1.- brightness

      red    = 1.
      green  = 0.0
      blue   = 1.0
      colfac1 = 1. - red
      colfac2 = 1. - green
      colfac3 = 1. - blue
      red =   1. - (evalue/amaxstr*colfac1)**btf
      green = 1. - (evalue/amaxstr*colfac2)**btf
      blue = 1. - (evalue/amaxstr*colfac3)**btf
      return
      end


c******    shade the elements with respect to specified stress

      subroutine shadste(i)

      include 'plot.cmn'

      if(stres) then

c         if(abs(abs(ste(i,ipv))/amaxst).ge. 1.) shade = 0.
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 1.) shade = .1
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.9) shade = .2
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.8) shade = .3
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.7) shade = .4
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.6) shade = .5
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.5) shade = .6
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.4) shade = .7
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.3) shade = .8
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.2) shade = .9
c         if(abs(abs(ste(i,ipv))/amaxst).lt. 0.1) shade = 1.


      call getrgbstr(ste(i,ipv))
       
c       write(6,*) red, green
c      write(10,1008) shade
c1008  format(' ','closepath gsave ',f4.1,' setgray fill
c     + grestore stroke ')

      write(10,1008) red, green, blue
1008  format(' ','closepath gsave ',3(1X,f4.1),' setrgbcolor fill
     + grestore stroke ')


      write(10,1009)
1009  format(' ','0.0 setgray')
      end if

      return

      end


c****  SET RGB COlOURS FOR STRESSES

      subroutine getrgbstr(esvalue)

      include 'plot.cmn'

      amaxstr = abs(aminst)
      if(abs(amaxst).gt.abs(aminst)) amaxstr = abs(amaxst)


      btf = 1.- brightness

      if(esvalue.lt.0.) then 
      red    = 1.
      green  = 0.
      blue   = 0.
      colfac1 = 1. - red
      colfac2 = 1. - green
      colfac3 = 1. - blue
      amaxstr = -1.*amaxstr
      red =   1. - (esvalue/amaxstr*colfac1)**btf
      green = 1. - (esvalue/amaxstr*colfac2)**btf
      blue = 1. - (esvalue/amaxstr*colfac3)**btf
      else
      red    = 0.
      green  = 1.
      blue   = 0.
      colfac1 = 1. - red
      colfac2 = 1. - green
      colfac3 = 1. - blue
      red =   1. - (esvalue/amaxstr*colfac1)**btf
      green = 1. - (esvalue/amaxstr*colfac2)**btf
      blue = 1. - (esvalue/amaxstr*colfac3)**btf
      endif
      return
      end




c****   plot cable elements
      
       subroutine pltcbl(sf)
   
       include 'plot.cmn'

       
       if(nl .gt. 0) then

       write(10,1011)
1011   format(' ','stroke newpath ')

       do  461 j = 1,nl
          n3 = nodl(j,2)
          x3 = sf*tcor(n3,1)*scalefac + shiftx
          y3 = sf*tcor(n3,2)*scalefac + shifty
          ix = int(x3)
          iy = int(y3)

          nq = nodl(j,1)    
          xq = sf*tcor(nq,1)*scalefac + shiftx
          yq = sf*tcor(nq,2)*scalefac + shifty
          ixq = int(xq)
          iyq = int(yq)
          write(10,1013)ix,iy, ixq,iyq
1013      format(' ',i5,i5,' m ',i5,i5,' l ')

461   continue
      write(10,*)' 1   setlinecap'
      write(10,*)' 20  setlinewidth stroke'
      end if
      return

      end


c***** plot geodesic strings

      subroutine pltgst(sf)

      include 'plot.cmn'

      if(nstr.gt. 0) then

      write(10,1014)
1014  format(' ','stroke newpath')       

      do 58 j = 1,nstr

         nsu = ngs(j)
         nr = nsnod(j,1)
         x3 = sf*tcor(nr,1)*scalefac + shiftx
         y3 = sf*tcor(nr,2)*scalefac + shifty
         ix = int(x3)
         iy = int(y3)
         write(10,1015) ix,iy
1015     format(' ',i5,i5,' m')

         do 456 jx = 2, nsu + 1
            nr = nsnod(j,jx)
            x3 = sf*tcor(nr,1)*scalefac + shiftx
            y3 = sf*tcor(nr,2)*scalefac + shifty
            ix = int(x3)
            iy = int(y3)
            write(10,1016) ix,iy
1016        format(' ',i5,i5,' m ')
456     continue
58    continue

      write(10,*)' closepath 3 setlinewidth'
      end if

      return
   
      end         



c****************************************************************************

                            subroutine closing

       include 'plot.cmn'

      if(err.or.stres) then
      write(10,1010)
      else
      write(10,1007)
      end if
1007  format(' ','  stroke 4.2 4.2 scale',/,' showpage ')
1010  format(' ',' showpage')

      return
   
      end




c******************************************************************************

                            subroutine data 

c******   read input data from input.dat
      include 'plot.cmn'
c****  data input routine                                               
      read (nin,10)titel                                                 
   10 format(a32)                                                       
      write(nout,20)titel                                                   
   20 format(      /   /,' data and results  ',/,1x,a32)                
      read(nin,30)dt,damp
   30 format(bz,2f10.5)
      write(nout,40)dt,damp
   40 format(/,' timestep  dt=',f10.5,5x,'damp=',f10.5)
      read(nin,50) nn,ne,nb,nm,nl,ndime,nstr,nbe,kmax,nmax,iln
      read(nin,50) nq
   50 format(bz,2i6,9i5)                                                   
      write(nout,60) nn,ne,nb,nm,nl,ndime,nstr,nbe,kmax,nmax                 
   60 format(/,' no of nodes=',i5,/,' no of elements=',i5,/,            
     1' no of restrained boundary nodes=',i5,/,' no of material types=',
     2i5,/,' no of links=',i5,/,' no of dimensions=',i5,/,' no of strs='
     3,i5,/,' no beams=',i5,/,' no of itgrps=',i5,/,' no of itbtwp=',i5)
      write(nout,65) nq
   65 format(/,' no of quadrilaterals=',i5)
      write(nout,70)iln                                                   
   70 format(' o if no loads=',i5)                                      
      do 80 l=1,nn                                                      
      read(nin,90) n,(cord(n,i),i=1,ndime)
   80 end do                                 
   90 format(bz,i5,3f15.5) 
      do 100 n=1,nn                                                     
      write(nout,110) n,(cord(n,i),i=1,ndime )                           
  100 end do
  110 format(' node=',i5,10x,'cords=',3(2x ,f10.4))                     
      if(ne.gt.0) then
        do 120 l=1,ne                                                     
        read(nin,130)nxx,(nod(l,i),i=1,4)                                     
  120   end do
  130   format(bz,5i5)                                                    
        write(nout,140)(n,(nod(n,i),i=1,4),n=1,ne)                             
  140   format(' element=',i3,2x,'node1=',i3,2x,'node2=',i3,2x,'node3='   
     1          ,i3,2x,'mat no=',i5)
      endif                                              
      if(nl.gt.0) then
        do 150 l=1,nl                                                     
        read(nin,160) n,(nodl(l,i),i=1,3)                                    
  150   end do
  160   format(bz,4i5)                                                    
        write(nout,170)(n,(nodl(n,i),i=1,3),n=1,nl)                           
  170   format(' link=',i3,2x,'node1=',i3,2x,'node2=',i3,2x,              
     1         'mat no=',i5)
      endif
      if(nq.gt.0) then
        write(nout,171) 
  171   format(/,' quadrilateral elements ',
     .   /,' element          nodal topology             material no')
        do 172   iq = 1,nq
        read(nin,174) n,(nodq(iq,i),i=1,5)                                    
  172   end do
  174   format(bz,6i5)                                                    
        write(nout,176)(n,(nodq(n,i),i=1,5),n=1,nq)
  176   format(3x,i5,5x,4i5,5x,i5)
      endif
      write(nout,180)
  180 format( / ,' material data')                                      
      do 190 l=1,nm                                                     
      read(nin,200) n,mat(n),(amat(n,i),i=1,5)                             
      write(nout,210) n,mat(n),(amat(n,i),i=1,5)
  190 continue
  200 format(bz,2i5,5e12.5)                                             
  210 format('  material number',i2,3x,i5,5x,5(3x,e12.5))               
      do 220 n=1,nb                                                      
      read(nin,230)(nco(n,i),i=1,(ndime+1))                                 
  220 end do
  230 format(bz,4i5)                                                    
      write(nout,240)                                                        
  240 format(/,' boundary conditions')                                  
      do 250 n=1,nb                                                     
      write(nout,260)(nco(n,i),i=1,(ndime+1))                               
  250 end do
  260 format(  ' node=',i5,5x,3(5x,i2))                                 
c*** read in the applied loads put in f array                           
c*** remember must have last line for last node even if no applied load 
      write(nout,270)                                                      
  270 format(///, ' loading conditions'/)                               
      if(iln.gt.0)  then
        iln=0                                                             
  280   iln=iln+1                                                         
        i=iln                                                             
        read(nin,290)nc(i),(aload(i,ii),ii=1,ndime)                          
  290   format(i5,3f15.5)                                                 
        write(nout,300)nc(i),(aload(i,ii),ii=1,ndime)                         
  300   format(1h ,'node=',i5,5x,3(4x,f15.5))                             
        if (nc(i).lt.nn) goto 280
      endif
      if(nstr.gt.0) call gsdata

c----- open and read error file

      if(err) then
      open(95,status = 'old',file = 'error.dat')
      do 2000 k = 1,ne
      read(95,*)errorr(k)
         if(norm.eq.1) then
         call maxerr
         else
         errorr(k) = errorr(k)/100.
         end if
2000  continue
      end if

c----- open and read stress file

      if(stres) then
      open(90,status = 'old', file = 'stress.dat')
      write(6,*)'enter 1 2 3 4 5 for sx sy sxy smax smin'
      read(5,*) ipv

      do 2001 k = 1,ne
      read(90,1500)(ste(k,l),l = 1,5)
2001  end do
1500  format(e12.5,4(1x,e12.5))
      end if      

      call maxstr

      return                                                            
      end        


c****************************************************************************

                         subroutine maxstr
      include 'plot.cmn'
      
      amaxst = -1.0e+30    

      do 2000 i = 1, ne
          if(amaxst.lt.ste(i,ipv)) amaxst = ste(i,ipv)
2000  end do

      aminst = amaxst
      
      do 3002 i = 1, ne
          if(aminst.gt.ste(i,ipv)) aminst = ste(i,ipv) 
3002  end do

      return
      end

c****************************************************************************

                         subroutine maxerr
      include 'plot.cmn'
      
      amaxer = errorr(1)    

      do 2000 i = 1, ne
          if(amaxer.lt.abs(errorr(i))) amaxer = abs(errorr(i))
2000  end do

      aminer = amaxer

      do 3002 i = 1, ne
          if(aminer.gt.abs(errorr(i))) aminer = abs(errorr(i))
3002  end do

      return
      end
                                                       
c******************************************************************************
 
                          subroutine gsdata

c******read data for geodesic strings

      include 'plot.cmn'

c**** reads in string data                                              dre15150
      do 58 j=1,nstr                                                    dre15160
      read(nin,781)ngs(j)                                               dre15170
  781 format(i5)                                                        dre15180
      write(nout,782)j,ngs(j)                                           dre15190
  782 format(' string no=',i5,'   no of links=',i5,/,                   dre15200
     .' string node:')                                                  dre15210
      read(nin,784)(nsnod(j,jd),jd=1,(ngs(j)+1))                        dre15220
  784 format (8i5)                                                   
      write(nout,785) (nsnod(j,jd),jd=1,(ngs(j)+1))                        dre15240
  785 format(1x,8i5)
   58 continue                                                          dre15250
      return                                                            dre15540
      end                                                               dre15550


c***** paint the legend
 
      subroutine index

      include 'plot.cmn'


      if(err.or.stres) then

      sc = 1./.24

      write(10,1018)sc,sc
1018  format(' ',f4.1,f4.1,' scale /times-roman findfont 12 scalefont
     + setfont 711 331 translate')     

      do 777 i = 0, 10
      ai = 1. - i/10.
      jj = i*10
      kk = (i+1)*10
      if(kk.ge.100) kk = 100
      write(10,1017)ai, jj, kk
1017  format(' ',' 1. setlinewidth ',f4.1,' setgray   0 0 m 
     +12 0 l 12 12 l 0 12 l closepath gsave',/,' fill grestore 0 setgray 
     + stroke 17 0 m (',i3,'-',i3,'%) show 0 15 translate')
777   continue

      end if

      return
      end 

c******    shade the elements with respect to specified parameter

      subroutine shadedom(j)

      include 'plot.cmn'

         icerr = 100*errorr(j) 
         
        
         if(icerr.eq.1) then
         red = 0.
         green = 0.
         blue = 1.
         goto 3000
         endif

         if(icerr.eq.2) then
         red = 0.
         green = 1.
         blue = 0.
         goto 3000
         endif

         if(icerr.eq.3) then
         red = 0.9
         green = 0.
         blue = 0.9
         goto 3000
         endif

         if(icerr.eq.4) then
         red = 0.9
         green = 0.9
         blue = 0.
         goto 3000
         endif


         if(icerr.eq.5) then
         red = 0.
         green = 0.9
         blue = 0.9
         goto 3000
         endif

         if(icerr.eq.6) then
         red = 0.9
         green = 0.9
         blue = 0.9
         goto 3000
         endif


         if(icerr.eq.7) then
         red = 0.9
         green = 0.7
         blue = 0.7
         goto 3000
         endif

         if(icerr.eq.8) then
         red = 0.9
         green = 0.6
         blue = 0.6
         goto 3000
         endif


         if(icerr.eq.9) then
         red = 0.9
         green = 0.4
         blue = 0.4
         goto 3000
         endif


         if(icerr.eq.10) then
         red = 0.8
         green = 0.8
         blue = 0.9
         goto 3000
         endif


         if(icerr.eq.11) then
         red = 0.2
         green = 0.6
         blue = 0.9
         goto 3000
         endif

         if(icerr.eq.12) then
         red = 0.7
         green = 0.4
         blue = 0.9
         goto 3000
         endif


         if(icerr.eq.13) then
         red = 0.7
         green = 0.9
         blue = 0.7
         goto 3000
         endif


         if(icerr.eq.14) then
         red = 0.7
         green = 0.9
         blue = 0.4
         goto 3000
         endif


         if(icerr.eq.15) then
         red = 0.5
         green = 0.7
         blue = 0.5
         goto 3000
         endif


         if(icerr.eq.16) then
         red = 0.5
         green = 0.5
         blue = 0.5
         goto 3000
         endif

   
3000  write(10,1008) red, green, blue
1008  format(' ','closepath gsave ',3(1x,f4.1),' setrgbcolor fill
     + grestore stroke ')
      write(10,1009)
1009  format(' ','0.0 setgray')
      

      return

      end



c***** paint the color legend for error and stress
 
      subroutine indexerrstr

      include 'plot.cmn'


      if(err.or.stres) then

      sc = 1./.24

      write(10,1018)sc,sc
1018  format(' ',f4.1,f4.1,' scale /times-roman findfont 12 scalefont
     + setfont 690 331 translate')     

      do 777 i = 0, 10

      if(stres) then 

      amaxstr = abs(aminst)
      if(abs(amaxst).gt.abs(aminst)) amaxstr = abs(amaxst)
      if (i.lt.5) ai = (5-i)/5.0*(-amaxstr)
      if(i.gt.5)  ai = (i-5.)/5.0*amaxstr  
      if(i.eq.5)  ai = 0     
      call getrgbstr(ai)

      else

      ai = i/10.*100.
      call getrgberr(ai)

      endif

      jj = i*10

      if(err) then
      write(10,1017)red, green, blue, jj
1017  format(' ',' 1. setlinewidth ',f4.1, f4.1, f4.1,
     +' setrgbcolor 0 0 m 
     +12 0 l 12 12 l 0 12 l closepath gsave',/,' fill grestore 0 setgray 
     + stroke 17 0 m (',i3,'%) show 0 15 translate')
      
      else

      write(10,2017)red, green, blue, ai
2017  format(' ',' 1. setlinewidth ',f4.1, f4.1, f4.1,
     +' setrgbcolor 0 0 m 
     +12 0 l 12 12 l 0 12 l closepath gsave',/,' fill grestore 0 setgray 
     + stroke 17 0 m (',e8.2,') show 0 15 translate')      

      endif

777   continue

      end if

      return
      end 
