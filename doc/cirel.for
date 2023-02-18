c compilar con /assume:byterecl en MS Developer Studio 97 para que asuma
c recl=2 en bytes y no en palabras de 4 bytes (subrutina readsad)
      use dflib
	use dfport
	parameter(deg=57.2957795,ispline=100,maxnp=4)
	parameter(maxnsta=25,maxndat=30000,maxmm=200)
	parameter(maxnslo=601,maxlwin=200)
	character*60 temp,farray,fdat1,fdat2
	logical resfile
	real ss(maxnp),ds(maxnp)
	real xsta(maxnsta),ysta(maxnsta),dsta(maxnsta)
	real dx(maxnsta,maxnsta),dy(maxnsta,maxnsta),dt(maxnsta,maxnsta)
	real dat(maxndat),data1(maxnsta,maxndat),data2(maxnsta,maxndat)
	real w(maxmm),t(maxmm),cc(maxmm),ccpen(maxmm)
	real sx(maxnslo),sy(maxnslo),ff(maxnslo,maxnslo)
	real tspline(maxmm*ispline),ccspline(maxmm*ispline)
	integer ii(2),ktime(7),msta(maxnsta),status*2
	data msam1/30/,msam2/30/
      data cut/1.16/
c NAME OF THE DATA FILE
      call getarg(1,temp,status)
        if(status.lt.0)stop 'Error: need input file'
c READ DATA FILE
      open(21,file=temp)
        read(21,'(a)')fdat1
        read(21,'(a)')fdat2
        read(21,'(a)')farray
	  read(21,*)slo0,azi0
        read(21,*)ni,nl,fl,fh
        read(21,'(a)')temp
	  np=0
	  do i=2,len_trim(temp)
	    if(temp(i-1:i-1).ne.' '.and.temp(i:i).eq.' ')np=np+1
        enddo
        np=(np+1)/2
	  if(np.gt.maxnp)stop 'ERROR: too many slowness domains'
        read(temp,*)(ss(i),ds(i),i=1,np)
        read(21,*)nini1,nini2,nwin,lwin,advance
c          if(nwin.gt.maxnwin)stop 'ERROR: too many windows'
          if(lwin.gt.maxlwin)stop 'ERROR: window too long'
          if(advance.gt..75)stop 'ERROR: small overlapping'
      close(21) 
c NAME OF THE OUTPUT FILE
      call getarg(2,temp,status)
        if(status.lt.0)stop 'ERROR: need ouput file'
      open(27,file=temp)
c NAME OF THE RESIDUAL FILE
      call getarg(3,temp,status)
        if(status.lt.0)then
          resfile=.false.
        else
          resfile=.true.
        endif  
      if(resfile)open(28,file=temp)
c READ ARRAY FILE (EAST,NORTH) [KM]
      open(21,file=farray)
	i=0
      do while(.not.eof(21))
	  i=i+1
        read(21,*,end=100)xsta(i),ysta(i)
      end do
100	nsta=i-1
	close(21)
c READ DATA
      call readsad(fdat1,nsta1,ndat1,ktime,fsam1,data1)
      call readsad(fdat2,nsta2,ndat2,ktime,fsam2,data2)
	write(*,*)nsta,nsta1,nsta2
	if(nsta.ne.nsta1)stop'ERROR: number of stations in array file'
      if(nsta1.ne.nsta2)stop'ERROR: different number of traces'
	if(fsam1.ne.fsam2)stop'ERROR: different sampling rate'
c BASELINE AND FILTER
      do 3 i=1,nsta
	  do 2 j=1,ndat1
2	    dat(j)=data1(i,j)
        call base(dat,ndat1,maxndat)
        call filter(dat,ndat1,maxndat,fl,fh,fsam1)
	  do 3 j=1,ndat1
3	    data1(i,j)=dat(j)
      do 5 i=1,nsta
	  do 4 j=1,ndat2
4	    dat(j)=data2(i,j)
        call base(dat,ndat2,maxndat)
        call filter(dat,ndat2,maxndat,fl,fh,fsam2)
	  do 5 j=1,ndat2
5	    data2(i,j)=dat(j)
c DELAYS OF MASTER EVENT ARRIVALS TO ARRAY STATIONS
      slo0x=slo0*sin(azi0/deg)
      slo0y=slo0*cos(azi0/deg)
      do 1 i=1,nsta
	  msta(i)=int(fsam1*(slo0x*(xsta(i)-xsta(1))
     +                    +slo0y*(ysta(i)-ysta(1))))
1     continue
c DEFINE WINDOW
      do 6 l=1,lwin
c BOX
	  w(l)=1.0
c        w=hanning(lwin)
c        w=ones(1,lwin);iw=round(0.1*lwin);
c        w(1:iw)=cos((1-([1:iw]-1)/iw)*pi/2);
c        w(lwin-[1:iw]+1)=cos(([iw:-1:1]/iw)*pi/2);
6     continue
c TIME VECTORS FOR DATA AND INTERPOLATED DATA
      mm=msam1+msam2+1
	mmspline=ispline*(mm-1)+1
      do 7 i=1,mm
        t(i)=float(i-1)/fsam1
7     continue
      do 8 i=1,mmspline
        tspline(i)=float(i-1)/ispline/fsam1
8     continue
c LOOP IN WINDOWS
      do 22 kk=1,nwin
        write(*,105)kk,' OF ',nwin
105     format(i4,a4,i4)
        kwin=(kk-1)*nint(lwin*advance)
	  ni1=nini1+kwin
	  ni2=nini2+kwin
c CROSS-CORRELATION COEFICIENT AS A FUNCTION OF DELAY
        do 13 ista=1,nsta
          c11=0.
          do 91 l=1,lwin
	      if(msta(ista)+ni1.ge.0)then 
              c11=c11+data1(ista,ni1+msta(ista)+l)
     +  	           *data1(ista,ni1+msta(ista)+l)
     +	           *w(l)
	      else
	        stop 'ERROR: no data in window'
	      endif
91        continue
c         do 10 m1=1,mm
          m1=msam1
          do 10 m2=1,mm
	  	  c12=0.
		  c22=0.
            do 9 l=1,lwin
              c12=c12+data1(ista,ni1+msta(ista)+l)
     +               *data2(ista,ni2+msta(ista)+m2-m1-1+l)
     +    	       *w(l)
              c22=c22+data2(ista,ni2+msta(ista)+m2-m1-1+l)
     +               *data2(ista,ni2+msta(ista)+m2-m1-1+l)
     +               *w(l)
9           continue
            cc(m2)=c12/sqrt(c11*c22)
10        continue
c INTERPOLATE CROSS-CORRELATION FUNCTION (CUBIC SPLINE)
          pen1=(cc(2)-cc(1))/(t(2)-t(1))
          penmm=(cc(mm)-cc(mm-1))/(t(mm)-t(mm-1))
          call spline(t,cc,mm,pen1,penmm,ccpen)
	    do 11 i=1,mmspline
	      call splint(t,cc,ccpen,mm,tspline(i),ccspline(i))
11        continue
c FIND MAXIMUM OF THE INTERPOLATED CROSS-CORRELATION
          ccsplinemax=0.
	    tsplinemax=0.
          do 12 i=1,mmspline
	      if(ccspline(i).gt.ccsplinemax)then
	        ccsplinemax=ccspline(i)
	        tsplinemax=tspline(i)
	      endif
12        continue
c DELAY BETWEEN EQS AT STATION ISTA
          dsta(ista)=tsplinemax-float(msam1)/fsam1
13      continue
c DIFFERENCES IN SPACE AND TIME AMONG STATIONS
        do 14 i=1,nsta-1
	  do 14 j=i+1,nsta
          dx(i,j)=xsta(i)-xsta(j)
	    dy(i,j)=ysta(i)-ysta(j)
          dt(i,j)=dsta(i)-dsta(j)
14      continue
c INITIALIZE
	  cte=1e3*sqrt(2.0/nsta/(nsta-1))
        s0x=0.
	  s0y=0.
c INITIALIZE RESIDUAL FF
        do 15 ix=1,maxnslo
        do 15 iy=1,maxnslo
          ff(ix,iy)=666.e6
15      continue
c LOOP IN SIZE AND STEP OF THE SLOWNESS DOMAIN
        do 19 k=1,np
c NUMBER OF NODES
	    ns=2*nint(ss(k)/ds(k))+1
	    if(ns.gt.maxnslo)stop 'ERROR: ns > maxnslo'
c CALCULATE RANGES OF EAST AND NORTH SLOWNESS
          do 16 is=1,ns 
            sx(is)=s0x-ss(k)+(is-1)*ds(k)
	      sx(is)=ds(np)*nint(sx(is)/ds(np))
            sy(is)=s0y-ss(k)+(is-1)*ds(k)
	      sy(is)=ds(np)*nint(sy(is)/ds(np))
16        continue
c EVALUATE RESIDUAL FOR EACH SLOWNESS IN THE DOMAIN
          do 18 ix=1,ns
          do 18 iy=1,ns
	      res=0.
            do 17 i=1,nsta-1
	      do 17 j=i+1,nsta
              res=res+(sx(ix)*dx(i,j)+sy(iy)*dy(i,j)-dt(i,j))**2
17          continue
c RESIDUAL [ms]
            ff(ix,iy)=cte*sqrt(res)
c WRITE TO RESIDUAL FILE          
            if(resfile.and.k.eq.np)
     +	   	write(28,110)slo0x+sx(ix),slo0y+sy(iy),ff(ix,iy)
110         format(2f9.5,f12.5)   
18        continue
c FIND MINIMUM RESIDUAL FOR CURRENT SLOWNESS DOMAIN
          ii=minloc(ff)
          s0x=sx(ii(1))
	    s0y=sy(ii(2))
19      continue
c CALCULATE SLOWNESS, AZIMUTH, AND MINIMUM RESIDUAL
        call r2p(slo0x+s0x,slo0y+s0y,slo,azi)
        ffmin=ff(ii(1),ii(2))
        tt=(nini1+.5*lwin+kwin)/fsam1
c WRITE RESULTS
        write(27,111)tt,slo,azi,ffmin
111     format(f10.3,1x,f7.4,1x,f6.1,1x,f10.4)
c        write(27,111)tt,s0x,s0y,slo,azi,ffmin
c111     format(f10.3,1x,f9.5,1x,f9.5,1x,f9.5,1x,f7.2,1x,f10.4)
c ERROR LIMITS
c        l=0
c        c=cut*ffmin
c        write(temp,'(a12,i2.2)')'cirel_error.',kk
c        open(26,file=temp)
c        do 21 ix=1,ns
c        do 21 iy=1,ns-1
c	     if(ff(ix,iy).gt.c.and.ff(ix,iy+1).le.c)then
c 	      l=l+1
c	      rx=slo0x+sx(ix)
c	      ry=slo0y+sy(iy)
c            write(26,112)rx,ry
c          elseif(ff(ix,iy).le.c.and.ff(ix,iy+1).gt.c)then
c	      l=l+1
c	      rx=slo0x+sx(ix)
c	      ry=slo0y+sy(iy+1)
c            write(26,112)rx,ry
c          endif
c21      continue
c        close(26)
c112     format(2f10.5)
c CLOSE LOOP IN WINDOWS
22    continue      
c CLOSE FILES
      if(resfile)close(28)
      close(27)
c END
	end