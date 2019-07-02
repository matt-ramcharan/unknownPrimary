cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SVM using OPTIMA                               c
c     Dataset: Sonar Classification                  c
c     Program version: V1.0                          c
c     Written by Colin Campbell (using a QP          c
c     routine by Andre Tits )                        c
c     Now with added kernels (Hillinger, etc)        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      program optima
c
      integer iwsize,nwsize,nparam,np,nf,nineq,neq,nn,nat
      character*50 fil,filtr,filts,filtar,filconf,filalp
c
c Use the program optpar.f to find the mimimum size
c of the workspaces below.
c
      parameter(np=1000)
      parameter(iwsize=1000)
      parameter(nwsize=100000)
      parameter(nat=1000)
c
c number of objective functions
c
      parameter(nf=1)
c
c number of inequality contraints
c
      parameter(nineq=0)
c
c number of equality constraints
c
      parameter(neq=1)
c
      integer iw(iwsize)
      double  precision alpha(np),bl(np),bu(np),
     *        f(nf+1),g(nineq+neq+1),w(nwsize),
     *        ri(nat,np)
c
      double precision ker(1000,1000),y(1000)     
c
      external obj32,cntr32,grob32,grcn32
c
      integer mode,iprint,miter,neqn,nineqn,inform
      integer n,isample,ikernel,ndegree
      integer ic,icount,infile,nl
      double precision bigbnd,eps,epsneq,udelta,dd
      double precision rtarget,target(np),arg,inp(nat,np)
      double precision z,cc,thsum,ss,out,errt,generr,gen
      double precision sigma,akappa,thres,con,cp,cn
      double precision objeps,objrep,gLgeps,lambda
      double precision sum(np),sump,sumn,bias,margin
c
      common /data/ker,y
      common /counter/ic,icount
c
      open(1,file='svm.in',status='old')
c
c (1) nparam=training sample size
c
       read(1,*) nparam
c
c (2) Test sample size
c
       read(1,*) isample
c
c (3) n=number of attributes
c
       read(1,*) n
c
c (4) Length of input file name
c
       read(1,*) nl
c
c (5) The name of the input file
c
       read(1,*) fil
       filtr=fil(:nl)//'.tr'
       filts=fil(:nl)//'.ts'
       filtar=fil(:nl)//'.out'
       filconf=fil(:nl)//'.details'
       filalp=fil(:nl)//'.alp'
c
c input/output files
c training data
       open(2,file=filtr,status='unknown')
c test data
       open(3,file=filts,status='unknown')
c target, actual output
       open(4,file=filtar,status='unknown')
c training error, test error, etc
       open(5,file=filconf,status='unknown')
c alpha distribution
       open(6,file=filalp,status='unknown')
c
       write(4,123) 
123    format('Pattern/target value/actual output')
       write(6,124)
124    format('Pattern/alpha value')
c
c Parameters used in the SVM
c
c ikernel = choice of kernel
c ikernel = 0 : Linear
c ikernel = 1 : Gaussians
c ikernel = 2 : polynomial
c ikernel = 3 : tanh
c ikernel = 4 : chi^2
c ikernel = 5 : L1 norm
c ikernel = 6 : Hillinger
c
c (6) C soft margin
c
       read(1,*) bigbnd
c
c (7) lambda soft margin
c 
       read(1,*) lambda
c
c (8) Choice of kernel
c
       read(1,*) ikernel
c
c (9) Parameter for Gaussian kernels
c
       read(1,*) sigma
c
c Parameter for polynomial kernels (degree)
c
       ndegree=3
c
c Parameters for tanh kernels
c
       akappa=1.0d0/dble(n)
       thres=-1.0d0
c
c iprint  = use 1 ... 5 for full diagnostics
c miter   = max. number of iterations
c bigbnd  = upper bound
c icount  = Output L ever icount calls
c eps     = stops if increase in L less than this amount
c
       mode=100
       iprint=0
       miter=1000
c       bigbnd=1.d+10
       epsneq=0.d0
       udelta=0.d0
       eps=1.0d-5
c
       ic=0
       icount=1
c
c End of specification of the problem
c*************************************************
c Array check
c
       if(nparam.gt.1000)then
         print*,'Dataset too large for program.'
         print*,'Use SMO algorithm.'
       endif
c
       cp=0.0D0
       cn=0.0D0
       do 10 i=1,nparam
         read(2,*) (ri(j,i),j=1,n),rtarget
         if(rtarget.gt.0.5d0)then
           y(i)=1.0d0
           cp=cp+1.0d0
         else
           y(i)=-1.0d0
           cn=cn+1.0d0
         endif
10     continue
c
       do 20 is=1,isample
         read(3,*) (inp(j,is),j=1,n),rtarget
         if(rtarget.gt.0.5d0)then
           target(is)=1.0d0
         else
           target(is)=-1.0d0
         endif
20     continue
c
c Find kernel matrix 
c
       do 30 im1=1,nparam
         do 40 im2=1,nparam
c
c ikernel=1 : Gaussian kernels
c          
           if(ikernel.eq.1)then
             arg=0.0d0
             do 50 j=1,n
               arg=arg+(ri(j,im1)-ri(j,im2))**2
50           continue
             arg=-arg/(2.0d0*sigma*sigma)
             if(arg.gt.-20.0d0)then
               ker(im1,im2)=dexp(arg)
             else
               ker(im1,im2)=0.0d0 
             endif
c
c ikernel=2 : Polynomial kernels
c
            elseif(ikernel.eq.2)then
              arg=0.0d0
              do 55 j=1,n
                arg=arg+ri(j,im1)*ri(j,im2)
55            continue
c              arg=arg/dble(n)
              ker(im1,im2)=arg**ndegree
c
c ikernel=3 : tanh updating function
c
            elseif(ikernel.eq.3)then
              arg=0.0d0
              do 65 j=1,n
                arg=arg+ri(j,im1)*ri(j,im2)
65            continue
              arg=akappa*arg+thres
              ker(im1,im2)=dtanh(arg)
c
c ikernel=4 : chi^2 kernel
c
            elseif(ikernel.eq.4)then
              arg=0.0d0
              do 75 j=1,n
                dd=(ri(j,im1)-ri(j,im2))**2
                dd=dd/(ri(j,im1)+ri(j,im2))
                arg=arg+dd
75            continue
              arg=-arg/(2.0d0*sigma*sigma)
              if(arg.gt.-20.0d0)then
                ker(im1,im2)=dexp(arg)
              else
                ker(im1,im2)=0.0d0
              endif
c
c ikernel=5 : L1 norm 
c
            elseif(ikernel.eq.5)then
              arg=0.0d0
              do 85 j=1,n  
                dd=ri(j,im1)-ri(j,im2)
                if(dd.gt.0.0d0)then
                  arg=arg+dd 
                else
                  arg=arg-dd
                endif
85            continue  
              arg=-arg/(2.0d0*sigma*sigma)
              if(arg.gt.-20.0d0)then
                ker(im1,im2)=dexp(arg)
              else
                ker(im1,im2)=0.0d0
              endif
c
c Hillinger   
c
            elseif(ikernel.eq.6)then
              arg=0.0d0  
              do 95 j=1,n
                dd=dsqrt(ri(j,im1))-dsqrt(ri(j,im2))
                arg=arg+dd*dd
95            continue
              arg=-arg/(2.0d0*sigma*sigma)
              if(arg.gt.-20.0d0)then
                ker(im1,im2)=dexp(arg)
              else
                ker(im1,im2)=0.0d0
              endif   
c
           endif  
40       continue
30     continue
c
c Add lambda soft margin
c
       do 35 im=1,nparam
         ker(im,im)=ker(im,im)+lambda
35     continue
c
c***********************************************
c Initialise alpha's so constraint satisfied.
c
      ss=0.0d0
      do 60 i=1,(nparam-1)
        if(y(i).gt.0.0D0)then
          alpha(i)=1.0d0/cp
          ss=ss+1.0d0/cp
        else
          alpha(i)=1.0d0/cn
          ss=ss+1.0d0/cn
        endif
        bl(i)=0.0d0
        bu(i)=bigbnd
60    continue
      if(y(nparam).gt.0.5d0)then
        alpha(nparam)=-ss
      else
        alpha(nparam)=ss
      endif
      bl(nparam)=0.0d0
      bu(nparam)=bigbnd
c
c**********************************
c Optimisation routine
c
      call FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,
     *           miter,inform,bigbnd,eps,epsneq,udelta,bl,bu,alpha,
     *           f,g,iw,iwsize,w,nwsize,obj32,cntr32,grob32,grcn32)
c
c******************************************
c (a) Training error
c
       z=0.0d0
       sump=1.0d10
       sumn=-1.0d10
c
       do 400 i=1,nparam
         ss=0.0d0
         do 410 j=1,nparam
             ss=ss+alpha(j)*y(j)*ker(j,i)
410      continue
         sum(i)=ss
c
         if(y(i).gt.0.5d0)then
           if(ss.lt.sump)then
             sump=ss
           endif
         else
           if(ss.gt.sumn)then
             sumn=ss
           endif
         endif
c
400    continue
c
       bias=0.5d0*(sump+sumn)
       margin=0.5d0*(sump-sumn)
c
       do 450 i=1,nparam
c
         thsum=sum(i)-bias
c         print*,y(i),thsum
         write(3,*) i,y(i),thsum
c
         if(thsum.gt.0.0d0)then
           out=1.0d0
         else
           out=-1.0d0
         endif       
c
         cc=y(i)*out
         if(cc.gt.0.0d0) then
           z=z+1.0d0
         endif
c
450    continue
c
       errt=1.0d0-z/dble(nparam)
c
c******************************************
c (b) Generalisation error
c
       gen=0.0d0
c
       do 500 is=1,isample
c
         ss=0.0d0
         do 510 i=1,nparam           
           if(ikernel.eq.1)then
             arg=0.0d0
             do 520 j=1,n
               arg=arg+(inp(j,is)-ri(j,i))**2
520          continue
             arg=-arg/(2.0d0*sigma*sigma)
             if(arg.gt.-20.0d0)then
               ss=ss+alpha(i)*y(i)*dexp(arg)
             endif
           elseif(ikernel.eq.2)then
             arg=0.0d0
             do 525 j=1,n
               arg=arg+inp(j,is)*ri(j,i)
525          continue
c             arg=arg/dble(n)             
             ss=ss+alpha(i)*y(i)*(arg**ndegree)
           elseif(ikernel.eq.3)then
             arg=0.0d0
             do 535 j=1,n
               arg=arg+inp(j,is)*ri(j,i)
535          continue
             arg=akappa*arg+thres
             ss=ss+alpha(i)*y(i)*dtanh(arg)
           elseif(ikernel.eq.4)then
             arg=0.0d0
             do 545 j=1,n
               dd=(inp(j,is)-ri(j,i))**2
               dd=dd/(inp(j,is)+ri(j,i))
               arg=arg+dd
545          continue
             arg=-arg/(2.0d0*sigma*sigma)
             if(arg.gt.-20.0d0)then  
               ss=ss+alpha(i)*y(i)*dexp(arg)
             endif
           elseif(ikernel.eq.5)then
             arg=0.0d0
             do 555 j=1,n
               dd=inp(j,is)-ri(j,i)
               if(dd.gt.0.0d0)then
                 arg=arg+dd
               else
                 arg=arg-dd
               endif
555          continue
             arg=-arg/(2.0d0*sigma*sigma)
             if(arg.gt.-20.0d0)then
               ss=ss+alpha(i)*y(i)*dexp(arg)
             endif
           elseif(ikernel.eq.6)then
             arg=0.0d0
             do 565 j=1,n
               dd=dsqrt(inp(j,is))-dsqrt(ri(j,i))
               arg=arg+dd*dd
565          continue
             arg=-arg/(2.0d0*sigma*sigma)
             if(arg.gt.-20.0d0)then
               ss=ss+alpha(i)*y(i)*dexp(arg)
             endif
c
           endif
510     continue
c     
         thsum=ss-bias         
         if(thsum.gt.0.0d0)then
           out=1.0d0
         else
           out=-1.0d0
         endif
c         print*,target(is),out,thsum
         write(3,*) target(is),thsum
         cc=target(is)*out
         if(cc.gt.0.0d0) then
           gen=gen+1.0d0
         endif
500    continue
       gen=gen/dble(isample)
       generr=1.0d0-gen
c
c ******** prints out on screen *************
c
       print*,'Training error=',errt
       print*,'Generalisation error=',generr
       print*,'Margin=',margin
       print*,'Bias=',bias
       write(5,221) nparam
221    format('Number of training patterns=',i8)
       write(5,222) n
222    format('Number of attributes=',i8)
       write(5,223) sigma
223    Format('Sigma value if Gaussian kernel used:',d16.8)
       write(5,224) ikernel
224    format('Choice of kernel=',i8)
       write(5,226) errt
226    format('Training error=',d16.8)
       write(5,227) generr
227    format('Test error=',d16.8)
c
c ********************************************
c
c Check alpha() all positive and equality condition
c satisfied.
c
       con=0.0D0
       do 800 i=1,nparam
         if(alpha(i).lt.0.0D0)then
           if(alpha(i).lt.-1.0D-5)then
             print*,'al negative',alpha(i),i
             pause
           else
             alpha(i)=0.0D0
           endif
         endif
         con=con+y(i)*alpha(i)
c         print*,alpha(i)
         write(6,*) i,alpha(i)
800    continue
       print*,'sum_i alpha_i y_i=',con
       write(5,*)
       write(5,228) con
228    format('sum_i alpha_i y_i=',d16.8)
c
       stop
       end
c
c************************************************
c Objective Function
c
      subroutine obj32(nparam,j,alpha,fj)
      integer nparam,j,k,l
      double precision alpha(nparam),fj
c
      double precision ker(1000,1000),y(1000)
c
      double precision f1,f2
      integer ic,icount
c
      common /data/ker,y
      common /counter/ic,icount
c
      f1=0.0D0
      f2=0.0D0
c
      do 10 k=1,nparam
        f1=f1+alpha(k)
        do 20 l=1,nparam
          f2=f2+alpha(k)*alpha(l)*y(k)*y(l)*ker(k,l)
20      continue
10    continue
c
c Progaram minimises objective function hence flip 
c sign to maximise
c
      fj=-(f1-0.5D0*f2)
c    
      ic=ic+1
      if(ic.eq.icount)then
        print*,'W= ',-fj
        ic=0
      endif
c
      return
      end
c
c**********************************************
c Gradients
c     
      subroutine grob32(nparam,j,alpha,gradfj,dummy)
      integer nparam,j,k
      double precision alpha(nparam),gradfj(nparam),
     *       f1,f2
c
      double precision ker(1000,1000),y(1000)
c
      external dummy
c
      common /data/ker,y
c
      do 10 i=1,nparam
        f1=0.d0
        do 20 k=1,nparam
          f1=f1+alpha(k)*y(k)*ker(k,i)
20      continue
        f1=y(i)*f1
c
c Program minimises hence flip sign to maximise 
c (see objective function above)
c
        gradfj(i)=-(1.d0-f1)
c
10    continue
      return
      end
c   
c**********************************
c Constraint conditions
c  
       subroutine cntr32(nparam,j,alpha,gj)
       integer nparam,j,k
       double precision alpha(nparam),gj
       double precision ker(1000,1000),y(1000)
c
       common /data/ker,y
c
       gj=0.0d0
       do 10 k=1,nparam
          gj=gj+alpha(k)*y(k)
10     continue
       return
       end
c
c*********************************
c Gradients of the constraint conditions 
c
      subroutine grcn32(nparam,j,alpha,gradgj,dummy)
      integer nparam,j
      double precision alpha(nparam),gradgj(nparam)
c
      double precision ker(1000,1000),y(1000)
c
      external dummy
c
      common /data/ker,y
c
      do 10 i=1,nparam
        gradgj(i)=y(i)
10    continue
c
      return
      end
c
      subroutine FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,
     *                  miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x,
     *                  f,g,iw,iwsize,w,nwsize,obj,constr,gradob,gradcn)
c                                                                      
c     implicit real*8(a-h,o-z)
      integer nparam,nf,neqn,nineqn,nineq,neq,mode,iprint,miter,inform,
     *        iwsize,nwsize
      integer iw(iwsize)
      double  precision bl(1),bu(1),x(1),
     *        f(1),g(1),w(1)
c     double  precision bl(nparam),bu(nparam),x(nparam),
c    *        f(nf),g(nineq+neq),w(nwsize)
c
c   When nf=0 and/or nineq+neq=0, the dimension for f and/or g
c   should be at least 1 due to F77 standard. See more details 
c   below.
c
      double  precision bigbnd,eps,epseqn,udelta
      external obj,constr,gradob,gradcn
c
c**********************************************************************c
c                                                                      c
c brief specification of various arrays and parameters in the calling  c
c sequence. See manual for more detailed description.                  c
c                                                                      c
c nparam : number of variables                                         c
c nf     : number of objective functions                               c
c nineqn : number of nonlinear inequality constraints                  c
c nineq  : number of inequality constraints                            c
c neqn   : number of nonlinear equality constraints                    c
c neq    : number of equality constraints                              c
c mode   : mode=CBA specifies job options as described below:          c
c          A = 0 : ordinary minimax problems                           c
c            = 1 : ordinary minimax problems with each individual      c
c                  function replaced by its absolute value, ie,        c
c                  an L_infty problem                                  c
c          B = 0 : monotone decrease of objective function             c
c                  after each iteration                                c
c            = 1 : monotone decrease of objective function after       c
c                  at most four iterations                             c
c          C = 1 : during line search, the function that rejected      c
c                  the previous step size is checked first;            c
c                  all functions of the same type ("objective" or      c
c                  "constraints") as the latter will then be checked   c
c                  first                                               c
c          C = 2 : all contraints will be checked first at every trial c
c                  point during the line search                        c
c iprint : print level indicator with the following options            c
c          iprint=0: no normal output except error information         c
c                    (this option is imposed during phase 1)           c
c          iprint=1:  a final printout at a local solution             c
c          iprint=2:  a brief printout at the end of each iteration    c
c          iprint=3:  detailed infomation is printed out at the end    c
c                     of each iteration for debugging purpose          c
c          iprint=10*N+M: N any positive integer, M=2 or 3.            c
c                     Information corresponding to iprint=M will be    c
c                     displayed at every 10*Nth iterations at the last c
c                     iteration                                        c
c miter  : maximum number of iterations allowed by the user to solve   c
c          the problem                                                 c
c inform : status report at the end of execution                       c
c          inform= 0:normal termination                                c
c          inform= 1:no feasible point found for linear constraints    c
c          inform= 2:no feasible point found for nonlinear constraints c
c          inform= 3:no solution has been found within miter iterations
c          inform= 4:stepsize is smaller than machine precision before c
c                    a successful new iterate is found                 c
c          inform= 5:failure of the QP solver in attempting to         c
c                    construct d0. A more robust QP solver may succeed.c
c          inform= 6:failure of the QP solver in attempting to         c
c                    construct d1. A more robust QP solver may succeed.c
c          inform= 7:inconsistent input data                           c
c          inform= 8:new iterate essentially identical to previous     c
c                    iterate, though stopping criteria not satisfied   c
c          inform= 9:penalty parameter too large, unable to satisfy    c
c                    nonlinear equality constraint                     c
c bigbnd : plus infinity                                               c
c eps    : stopping criterion that ensures at a solution, the norm of  c
c          the Newton direction vector is smaller than eps             c
c epseqn : tolerance of the violation of nonlinear equality constraintsc
c          allowed by the user at an optimal solution                  c
c udelta : perturbation size in computing gradients by finite          c
c          difference and the true perturbation is determined by       c
c          sign(x_i) X max{udelta, rteps X max{1, |x_i|}} for each     c
c          component of x, where rteps is the square root of machine   c
c          precision                                                   c
c bl     : array of dimension nparam,containing lower bound of x       c
c bu     : array of dimension nparam,containing upper bound of x       c
c x      : array of dimension nparam,containing initial guess in input c
c          and final iterate at the end of execution                   c
c f      : array of dimension max{1,nf}, containing objective values   c
c          at x in output                                              c
c g      : array of dimension max{1,nineq+neq}, containing constraint  c
c          values at x in output                                       c
c iw     : integer working space of dimension iwsize                   c
c iwsize : length of integer array iw                                  c
c w      : double precision working space of dimension nwsize.         c
c          at output, it contains lagrange multipliers                 c
c nwsize : length of double precision array w                          c
c obj    : subroutine that returns the value of objective functions    c
c          one upon each call                                          c
c constr : subroutine that returns the value of constraints            c
c          one upon each call                                          c
c gradob : subroutine that computes gradients of f, alternatively      c
c          it can be replaced by grobfd that computes finite           c
c          difference approximations                                   c
c gradcn : subroutine that computes gradients of g, alternatively      c
c          it can be replaced by grcnfd that computes finite           c
c          difference approximations                                   c
c                                                                      c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c                                                                      c
c  The purpose of FFSQP is to solve general nonlinear constrained      c
c  minimax optimization problems of the form                           c
c                                                                      c
c   (A=0 in mode)     minimize    max_i f_i(x)   for i=1,...,n_f       c
c                        or                                            c
c   (A=1 in mode)     minimize    max_j |f_i(x)|   for i=1,...,n_f     c
c                       s.t.      bl   <= x <=  bu                     c
c                                 g_j(x) <= 0,   for j=1,...,nineqn    c
c                                 A_1 x - B_1 <= 0                     c
c                                                                      c
c                                 h_i(x)  = 0,   for i=1,...,neqn      c
c                                 A_2 x - B_2  = 0                     c
c                                                                      c
c                                                                      c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c
      integer i,io,ipp,iter,j,ncallg,ncallf,ncnstr,nclin,nctotl,leniw,
     *        lenw,nwx,nwbl,nwbu,nwgrg,nwgpf,nwpenp,nwa,nwcvec,nwhess,
     *        nwcla,nww,nrowa,modd,nppram,iwnob,iwncn,iwia,iwisp,iwist,
     *        iwiw,nwdi,nwd,nwff,nwgrf,nwclla,nwhes1,nwsp,nwbak,nwsg,M,
     *       maxit,nob,nobL,nnineq,info,idummy,nshift,max0,modem,lstype,
     *       nstop,initvl,nn,nnn,nwgm,ipsp,ipspan,ipyes,iprnto,mod
      double  precision epsmac,QLeps,small,xi,gi,gmax,dummy,big,tolfea,
     *        rteps,epskt,upert,valnom,dsqrt,dmax1
      logical feasbl,feasb,prnt,nspace,Linfty,nAD,rolis1,d0is0
      common  /fsqpp1/nnineq,M,ncallg,ncallf,modd,lstype,nstop,
     *        /fsqpp2/io,ipp,ipsp,ipyes,info,idummy,iter,initvl,
     *        /fsqpp3/epsmac,rteps,upert,valnom,
     *        /fsqpp4/rolis1,d0is0,
     *        /fsqpq1/big,tolfea,/fsqpq2/maxit
      common  /CMACHE/QLeps
c
c     compute the machine precision
c
      io=6
      rolis1=.false.
c     iwiw=6*nparam+8*max0(1,nineq+neq)+7*max0(nf,1)+30
c     i=nineq+neq+1
c     nww=4*nparam**2+5*i*nparam+3*(nf+1)*nparam+26*(nparam+nf+1)
c    *   +45*i+100
c     if(iwsize.ge.iwiw.and.nwsize.ge.nww) goto 10
c       if(iwsize.lt.iwiw) write(io,9906) iwiw
c       if(nwsize.lt.nww)  write(io,9907) nww
c       info=7
c       goto 9000
c
 10   iter=0
      nstop=1
      nn=nineqn+neqn
      epsmac=small()
      QLeps=epsmac
      tolfea=epsmac*1.d+02
      big=bigbnd
      rteps=dsqrt(epsmac)
      upert=udelta
c
      i=mod(iprint,10)
      ipspan=max0(iprint-i,1)
      iprnto=iprint
      if(iprint.ge.10) iprint=i
      if(iprint.lt.2) ipspan=1
      if(ipspan.lt.10) ipyes=0
      nob=0
      gmax=-bigbnd
      info=0
      ipsp=ipspan
      ipp=iprint
      ncnstr=nineq+neq
      nnineq=nineq
c
c     check input data 
c
      if(iprint.gt.0) write(io,9900)
      call check(nparam,nf,Linfty,nAD,nineq,nineqn,neq,neqn,
     *           mode,modem,nwa,eps,bigbnd,bl,bu)
      if(info.eq.7) goto 9000
      lstype=nwa
c
      maxit=max0(max0(miter,10*max0(nparam,ncnstr)),1000)
      feasbl=.true.
      feasb=.true.
      prnt=.false.
      nspace=.false.
      nppram=nparam+1
      nshift=nparam**2+nppram**2
c
c  check whether x is within bounds
c
      do 100 i=1,nparam
        xi=x(i)
        if(bl(i).le.xi.and.bu(i).ge.xi) goto 100
        feasbl=.false.
        goto 110
 100  continue
 110  nclin=ncnstr-nn
c
c  check whether linear constraints are feasible
c
      if(nclin.eq.0) goto 210
      do 200 i=1,nclin
        j=i+nineqn
        if(j.le.nineq) then
          call constr(nparam,j,x,gi)
          if(gi.le.epsmac) goto120
          feasbl=.false.
        else if(j.gt.nineq) then
          call constr(nparam,j+neqn,x,gi)
          if(dabs(gi).le.epsmac) goto 120
          feasbl=.false.
        endif
 120    g(j)=gi
 200  continue
 210  if(feasbl) goto 240
      if(iprint.le.0) goto 230
        write(io,9901)
        call sbout1(io,nparam,'x                   ',dummy,x,2,1)
        prnt=.true.
 230  nctotl=nparam+nclin
      leniw=max0(2*nparam+2*nctotl+3,2*nclin+2*nparam+6)
      if(leniw.le.iwsize)then
        leniw=iwsize
      else
        write(io,9906) leniw
        info=7
        nspace=.true.
      endif
      nwx=1
      nwbl=nwx+nparam
      nwbu=nwbl+nctotl+4
      nwgrg=nwbu+nctotl+2
      nwa=nwgrg+nclin*nparam+1
      nwcvec=nwa+nparam*nclin+1
      nwhess=nwcvec+nparam
      nwcla=nwhess+nparam*nparam
      nww=nwcla+nctotl+nparam
      lenw=2*nparam**2+10*nparam+2*nctotl+1
      if((nww+lenw).le.nwsize) then
        lenw=nwsize-nww
        if(.not.nspace) goto 235
        write(io,9909)
        goto 9000
      else
        write (io,9907) nww+lenw
        write(io,9909)
        info=7
        goto 9000
      endif
c
c     attempt to generate a point satisfying all linear constraints
c 
 235  nrowa=max0(nclin,1)
      call initpt(nparam,nineqn,neq,neqn,nclin,nctotl,nrowa,x,bl,bu,
     *            iw,leniw,w(nwx),w(nwbl),w(nwbu),g(nineqn+1),w(nwgrg),
     *            w(nwa),w(nwcvec),w(nwhess),w(nwcla),w(nwbl+nparam+3),
     *            w(nww),lenw,constr,gradcn)
      if(info.ne.0) goto 9000
 240  do 245 i=1, neq-neqn
 245    g(nineq+neqn+i)=g(nineq+i) 
      if(nn.ne.0) goto 510
      goto 605
c
 290    do 300 i=1,nob
 300      w(i+nineqn+nshift)=w(i+nshift)
        nob=0
c
 510  continue
      if(info.eq.-1) goto 540
        do 520 i=1,nineqn
          call constr(nparam,i,x,w(i+nineqn+nshift))
          if(w(i+nineqn+nshift).gt.0.d0) feasb=.false.
 520    continue
        ncallg=nineqn
        if(feasb) goto 540
c
c     set indxob(i) in phase 1
c
        do 530 i=1,nineqn
          nob=nob+1
          iw(nob)=i
          w(nob+nshift)=w(i+nineqn+nshift)
          gmax=dmax1(gmax,w(nob+nshift))
 530    continue
      goto 580
 540    do 550 i=1,nineqn
          g(i)=w(i+nineqn+nshift)
          iw(nineqn+i+1)=i
 550    continue
        do 560 i=1,neq-neqn
          g(i+nineq+neqn)=g(i+nineq)
 560    continue
        do 570 i=1,neqn
          j=i+nineq
          call constr(nparam,j,x,g(j))
          iw(nineqn+nineqn+i+1)=j
 570    continue
        ncallg=ncallg+neqn
 580  continue
c
 605  if(iprint.le.0.or..not.feasb.or.prnt) goto 610
        write(io,9902)
        call sbout1(io,nparam,'x                   ',dummy,x,2,1)
        prnt=.true.
 610  if(nob.ne.0) goto 620
      if(iprint.le.0) goto 615
        if(info.eq.0) goto 613
        write(io,9904) ncallg
        if(ipp.eq.0) write(io,9910) iter
        if(ipp.gt.0) write(io,9910) iter-1
        if(ipp.eq.0) iter=iter+1
 613    if(.not.feasb.or.feasbl) goto 614
          write(io,9903)
          call sbout1(io,nparam,'x                   ',dummy,x,2,1)
 614    if(info.eq.0.and.prnt.and.feasb) goto 615
          write(io,9903)
          call sbout1(io,nparam,'x                   ',dummy,x,2,1)
 615  feasb=.true.
      feasbl=.true.
 620  nspace=.false.
      if(ipp.le.0.or.feasb.or.prnt) goto 630
        write(io,9901)
        call sbout1(io,nparam,'x                   ',dummy,x,2,1)
        prnt=.true.
 630  if(nob.eq.0) nob=1
c
c     set indxcn(1)--indxcn(ncnstr)
c
      if(feasb) nnn=nn
      if(.not.feasb) nnn=0
      do 700 i=1,nnn
 700    iw(nob+i)=iw(nineqn+i+1)
 710  do 800 i=1,nineq-nineqn
 800    iw(nob+nnn+i)=nineqn+i
      do 805 i=1,neq-neqn
        if(feasb) iw(nob+nineq+neqn+i)=nineq+neqn+i
        if(.not.feasb) iw(nineq+i)=nineq+neqn+i
 805  continue
      if(.not.feasb) goto 810
        nob=nf
        info=0
        ipp=iprint
        ipsp=ipspan
        modd=modem
        epskt=eps
        if(Linfty) nobL=2*nob
        if(.not.Linfty) nobL=nob
        if(nob.ne.0.or.neqn.ne.0) goto 910
        write(io,9908)
      goto 9000
 810    ipp=0
        ipsp=1
        modd=0
        nobL=nob
        info=-1
        epskt=1.d-10
 910  nctotl=nppram+ncnstr+max0(nobL,1)
      iwnob=1
      if(feasb) iwncn=iwnob+1
      if(.not.feasb) iwncn=iwnob+nob
      iwia=iwncn+ncnstr
      iwisp=iwia+nn+max0(nob,1)
      iwist=iwisp+nnineq-nineqn+1
      iwiw=iwist+nn+max0(nob,1)
      leniw=2*(ncnstr+max0(nobL,1))+2*nppram+6
c
      if((iwiw+leniw).le.iwsize) then
        leniw=iwsize-iwiw
      else
        write (io,9906) iwiw+leniw
        info=7
        nspace=.true.
      endif
      M=4
      if(modem.eq.1.and.nn.eq.0) M=3
      nwhess=1
      nwhes1=nwhess+nparam**2
      nwff=nwhes1+nppram**2
      nwx=nwff+max0(nob,1)+1
      nwdi=nwx+nppram
      nwd=nwdi+nppram
      nwgm=nwd+nppram
      nwgrg=nwgm+max0(1,4*neqn)
      nwgrf=nwgrg+ncnstr*nparam+1
      nwgpf=nwgrf+nparam*max0(nob,1)+1
      nwpenp=nwgpf+nparam
      nwa=nwpenp+neqn+1
      nwbl=nwa+(ncnstr+max0(nobL,1))*(nppram+1)
      nwbu=nwbl+nctotl+4
      nwcla=nwbu+nctotl+2
      nwclla=nwcla+nctotl+nppram
      nwcvec=nwclla+nctotl
      nwsp=nwcvec+nppram
      nwbak=nwsp+M+1
      nwsg=nwbak+max0(nob,1)+ncnstr+1
      nww=nwsg+neqn+1
      lenw=2*nppram*nppram+10*nppram+6*(ncnstr+max0(nobL,1)+1)
c
      if((nww+lenw).le.nwsize) then
        lenw=nwsize-nww
        if(.not.nspace) goto 920
        write(io,9909)
        goto 9000
      else
        write (io,9907) nww+lenw
        write(io,9909)
        info=7
        goto 9000
      endif
c
 920  do 1000 i=nwx,nwx+nparam-1
 1000   w(i)=x(i-nwx+1)
      w(nwx+nparam)=gmax
      if(.not.feasb) goto 1150
        do 1100 i=1,neqn
          if(g(i+nineq).gt.0d0) w(nwsg+i-1)=-1.d0
          if(g(i+nineq).le.0d0) w(nwsg+i-1)=1.d0
 1100   continue
c
c     either attempt to generate a point satisfying all constraints 
c     or try to solve the original problem
c
 1150 nrowa=max0(ncnstr+max0(nobL,1),1)
      call FFSQP1(miter,nparam,nob,nobL,nineqn,neq,neqn,ncnstr,nctotl,
     *            nrowa,feasb,epskt,epseqn,bl,bu,iw(iwnob),iw(iwncn),
     *            iw(iwia),iw(iwisp),iw(iwist),iw(iwiw),leniw,w(nwx),
     *            w(nwdi),w(nwd),g,w(nwgm),w(nwgrg),w(nwff),w(nwgrf),
     *            w(nwgpf),w(nwpenp),w(nwa),w(nwbl),w(nwbu),w(nwcla),
     *            w(nwclla),w(nwcvec),w(nwbl+nparam+3),w(nwhess),
     *            w(nwhes1),w(nwsp),w(nwbak),w(nwsg),w(nww),lenw,
     *            obj,constr,gradob,gradcn)
      do 1200 i=1,nparam
 1200   x(i)=w(nwx+i-1)
      if(info.eq.-1) goto 290
      if(info.eq.0.or.feasb) goto 1220
        info=2
        write(io,9905)
      goto 9000
 1220 do 1300 i=1,nf
 1300   f(i)=w(nwff+i-1)
      if(nobL.le.1) idummy=0
      if(nobL.gt.1) idummy=1
      if(nf.eq.1) nob=0
      do 1400 i=1,nparam+ncnstr+nob
        j=i
        if(i.gt.nparam.and.i.le.(nparam+ncnstr)) 
     *    j=nparam+iw(iwncn+i-nparam)
        if(i.le.nparam) then
          w(i)=w(nwclla+j-1)
        else if(i.gt.nparam) then
          if(i.le.(nparam+ncnstr)) j=nparam+iw(iwncn+i-1-nparam)
          w(i)=w(nwclla+j-1+idummy)
        endif
 1400 continue
      if(nf.eq.1) w(nparam+ncnstr+1)=1.d0
c
 9000 inform=info
      iprint=iprnto
      return
 9900 format(1x,// 1x,'   FFSQP Version 3.7b  (Released January 1998)'
     *   /   1x,'            Copyright (c) 1989 --- 1998         '
     *   /   1x,'               J.L. Zhou, A.L. Tits,            '
     *   /   1x,'                 and C.T. Lawrence              '
     *   /   1x,'                All Rights Reserved             ',//)
 9901 format(1x,'The given initial point is infeasible for inequality',
     *       /10x,'constraints and linear equality constraints:')

 9902 format(1x,'The given initial point is feasible for inequality',
     *   /8x,'constraints and linear equality constraints:')
 9903 format(1x,'Starting from the generated point feasible for',
     *   ' inequality',
     *   /10x,'constraints and linear equality constraints:')
 9904 format(1x,'To generate a point feasible for nonlinear inequality',
     *   /1x,'constraints and linear equality constraints,',
     *   ' ncallg = ',i10)
 9905 format(1x,'Error: No feasible point is found for nonlinear',
     *   ' inequality',
     *    /8x,'constraints and linear equality constraints'/)
 9906 format(1x,'iwsize should be bigger than', i20)
 9907 format(1x,'nwsize should be bigger than', i20)
 9908 format(1x,'current feasible iterate with no objective specified'/)
 9909 format(1x,/)
 9910 format(43x,'iteration = ',i10)
      end
      subroutine FFSQP1(miter,nparam,nob,nobL,nineqn,neq,neqn,ncnstr,
     *                  nctotl,nrowa,feasb,epskt,epseqn,xl,xu,indxob,
     *                  indxcn,iact,iskip,istore,iw,leniw,x,di,d,g,gm,
     *                  gradg,f,gradf,grdpsf,penp,a,bl,bu,clamda,
     *                  cllamd,cvec,bj,hess,hess1,span,backup,signeq,
     *                  w,lenw,obj,constr,gradob,gradcn)
c
c     FFSQP : main routine for the optimization
c
c     implicit real*8(a-h,o-z)
      integer miter,nparam,nob,nobL,nineqn,neq,neqn,ncnstr,nctotl,nrowa,
     *        leniw,lenw
      integer indxob(1),indxcn(1),iact(1),iskip(1),
     *        istore(1),iw(leniw)
c     integer indxob(nob),indxcn(ncnstr),iact(nob+nineqn+neqn),iskip(1),
c    *        istore(nineqn+nob+neqn),iw(leniw)
      double  precision epskt,epseqn
      double  precision xl(nparam),xu(nparam),x(nparam+1),di(nparam+1),
     *        d(nparam+1),g(1),gm(1),gradg(nparam,1),
     *        f(1),gradf(nparam,1),grdpsf(nparam),penp(1),
     *        a(nrowa,1),bl(1),bu(1),clamda(1),
     *        cllamd(1),cvec(nparam+1),bj(1),
     *        hess(nparam,nparam),hess1(1),span(1),
     *        backup(1),signeq(1),w(lenw)
c     double  precision xl(nparam),xu(nparam),x(nparam+1),di(nparam+1),
c    *        d(nparam+1),g(ncnstr),gm(4*neqn),gradg(nparam,ncnstr),
c    *        f(nob),gradf(nparam,nob),grdpsf(nparam),penp(neqn),
c    *        a(nrowa,1),bl(nctotl),bu(nctotl),clamda(nctotl+nparam+1),
c    *        cllamd(nctotl),cvec(nparam+1),bj(nrowa),
c    *        hess(nparam,nparam),hess1(nparam+1,nparam+1),span(1),
c    *        backup(nob+ncnstr),signeq(neqn),w(lenw)
      external obj,constr,gradob,gradcn
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,info,ipd,iter,nstop,
     *        initvl,ipspan,ipyes,lstype
      double  precision bigbnd,tolfea,epsmac,rteps,udelta,valnom
      logical dlfeas,local,update,first,rolis1,d0is0
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,iter,initvl,
     *        /fsqpp3/epsmac,rteps,udelta,valnom,
     *        /fsqpp4/rolis1,d0is0,
     *        /fsqpq1/bigbnd,tolfea,
c    *        /fsqp1/rentry,
     *        /fsqplo/dlfeas,local,update,first
c
c     bj(1+) is equivalent to bl(nparam+3+)
c
      integer i,ng,iskp,nfs,ncf,ncg,nn,non,nstart,nrst,ncnst1,nctot1
      double  precision Cbar,Ck,dbar,fmax,fM,fMp,steps,d0nm,dummy,
     *        sktnom,scvneq,grdftd,dmax1,psf
c
      initvl=1
      first=.true.
      nrst=0
      ipd=0
      if(nobL.gt.1) ng=1
      if(nobL.le.1) ng=0
      if(iter.eq.0) call diagnl(nparam,1.d0,hess)
      if(.not.feasb) goto 5
        first=.true.
        if(iter.gt.0) iter=iter-1
        if(iter.ne.0) call diagnl(nparam,1.d0,hess)
 5    Cbar=1.d-02
      Ck=Cbar
      dbar=5.0d0
      nstart=1
      ncallf=0
      nstop=1
      nfs=0
      non=miter
      if(mode.eq.0) goto 10
        nfs=M
        non=0
 10   if(feasb) then
        nn=nineqn+neqn
        ncnst1=ncnstr
        nctot1=nctotl
      else 
        nn=0
        ncnst1=ncnstr-nineqn-neqn
        nctot1=nnineq-nineqn+neq-neqn+nparam
        if(nob.gt.1) nctot1=nctot1+1
      endif
      scvneq=0.d0
      do 100 i=1,ncnst1
        valnom=g(indxcn(i))
        backup(i)=valnom
        if(feasb.and.i.gt.nineqn.and.i.le.nn) then
          gm(i-nineqn)=valnom*signeq(i-nineqn)
          scvneq=scvneq+dabs(valnom)
        endif
        if(.not.feasb.or.i.gt.nn) goto 20
          iact(i)=indxcn(i)
          istore(i)=0
          if(i.gt.nineqn) penp(i-nineqn)=2.d0
 20     call gradcn(nparam,indxcn(i),x,gradg(1,indxcn(i)),constr)
 100  continue
      call nullvc(nparam,grdpsf)
      psf=0.d0
      if(.not.feasb.or.neqn.eq.0) goto 110
        call resign(nparam,neqn,psf,grdpsf,penp,g(nnineq+1),
     *              gradg(1,nnineq+1),signeq,12,12)
 110  fmax=-bigbnd
      if(feasb.and.nob.eq.0) then
        fmax=0.d0
        fMp=0.d0
      endif
      do 140 i=1,nob
        if(.not.feasb) goto 120
          iact(nn+i)=i
          istore(nn+i)=0
          call obj(nparam,i,x,f(i))
          valnom=f(i)
          backup(i+ncnst1)=valnom
          call gradob(nparam,i,x,gradf(1,i),obj)
          ncallf=ncallf+1
          if(nobL.ne.nob) fmax=dmax1(fmax,-f(i))
        goto 130
 120      valnom=f(i)
          iact(i)=i
          istore(i)=0
          call gradcn(nparam,indxob(i),x,gradf(1,i),constr)
 130    fmax=dmax1(fmax,f(i))
 140  continue
      fM=fmax
      fMp=fM-psf
      span(1)=fM
c
      if(iprint.lt.3.or..not.first.or.ipyes.gt.0) goto 600
        do 300 i=1,nob
          if(.not.feasb) goto 250
            if(nob.gt.1) 
     *        call sbout2(io,nparam,i,'gradf(j,',')',gradf(1,i))
            if(nob.eq.1) 
     *        call sbout1(io,nparam,'gradf(j)            ',
     *                    dummy,gradf(1,1),2,2)
          goto 300
 250        call sbout2(io,nparam,indxob(i),'gradg(j,',')',gradf(1,i))
 300    continue
 310    if(ncnstr.eq.0) goto 410
        do 400 i=1,ncnst1
 400      call sbout2(io,nparam,indxcn(i),'gradg(j,',')',
     *                gradg(1,indxcn(i)))
        if(neqn.eq.0) goto 410
        call sbout1(io,nparam,'grdpsf(j)           ',dummy,grdpsf,2,2)
        call sbout1(io,neqn,'P                   ',dummy,penp,2,2)
 410    do 500 i=1,nparam
 500      call sbout2(io,nparam,i,'hess (j,',')',hess(1,i))
c
c     main loop of the algorithm
c
 600  nstop=1
 601  continue
        call out(miter,nparam,nob,nineqn,nn,neqn,ncnst1,x,g,
     *           f,fmax,fM,psf,steps,sktnom,d0nm,feasb)
        if(nstop.ne.0) goto 810
        if(.not.feasb) goto 809
          do 700 i=1,ncnst1
 700        g(i)=backup(i)
          do 800 i=1,nob
 800        f(i)=backup(i+ncnst1)
          if(neqn.eq.0) goto 809
            do 805 i=1,neqn
              cllamd(nparam+ng+nnineq+i)=(cllamd(nparam+ng+nnineq+i)-
     *          penp(i))*signeq(i)
 805        continue
 809      return
 810    continue
        if(ipspan.ge.10.and.iprint.ge.2.and.ipyes.eq.0) 
     *    write(io,9900) iter
        call dir(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnst1,nctot1,nrowa,
     *           feasb,steps,epskt,epseqn,sktnom,scvneq,Ck,d0nm,grdftd,
     *           xl,xu,indxob,indxcn,iact,iskp,iskip,istore,iw,leniw,
     *           x,di,d,g,gradg,f,fmax,fM,fMp,psf,gradf,grdpsf,penp,a,
     *           bl,bu,clamda,cllamd,cvec,bj,hess,hess1,w,lenw,
     *           backup,signeq,obj,constr)
        if(nstop.eq.0) goto 601
        first=.false.
        if(update.or.d0is0) goto 820
        call step(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnst1,ncg,ncf,
     *            indxob,indxcn,iact,iskp,iskip,istore,feasb,grdftd,
     *            f,fmax,fM,fMp,psf,penp,steps,scvneq,bu,x,di,d,g,w,
     *            backup,signeq,obj,constr)
        if(nstop.eq.0) goto 601
 820    call hesian(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnst1,nctot1,
     *              nfs,nstart,feasb,bigbnd,bu,x,f,fmax,fM,fMp,psf,
     *              gradf,grdpsf,penp,g,gm,gradg,indxob,indxcn,cllamd,
     *              bl,clamda,di,hess,d,steps,nrst,signeq,span,
     *              obj,constr,gradob,gradcn,
     *              hess1,cvec,bj,w,lenw,iw,leniw) 
        if(nstop.eq.0 .or. mode.eq.0) goto 601
        if(d0nm.gt.dbar) Ck=dmax1(dble(0.5*Ck),Cbar)
        if(d0nm.le.dbar.and.dlfeas) Ck=Ck
        if(d0nm.le.dbar.and..not.dlfeas.and..not.rolis1) Ck=10.0*Ck
      goto 601
 9900 format(1x,9hiteration,t22,i22)
      end
      subroutine check(nparam,nf,Linfty,nAD,nineq,nnl,neq,neqn,
     *                 mode,modem,lstype,eps,bigbnd,bl,bu)
c
c     FFSQP : check input data
c
c     implicit real*8(a-h,o-z)
      integer nparam,nf,nineq,nnl,neq,neqn,mode,modem,lstype
      double  precision bigbnd,eps
      double  precision bl(nparam),bu(nparam)
      logical Linfty,nAD
c
      integer io,iprint,ipspan,ipyes,info,idum1,idum2,idum3
      double  precision epsmac,dummy1,dummy2,dummy3
      common  /fsqpp2/io,iprint,ipspan,ipyes,info,idum1,idum2,idum3,
     *        /fsqpp3/epsmac,dummy1,dummy2,dummy3
c
      integer i
      double  precision bli,bui
c
      if (nparam.le.0)  
     *  call error('nparam should be positive!              ',info,io)
      if (nf.lt.0)     
     *  call error('nf     should not be negative!          ',info,io)
      if (nnl.lt.0)     
     *  call error('nineqn should not be negative!          ',info,io)
      if (nineq.lt.nnl) 
     *  call error('nineq  should be no smaller than nineqn!',info,io)
      if (neqn.lt.0)  
     *  call error('neqn   should not be negative!          ',info,io)
      if (neq.lt.neqn)  
     *  call error('neq    should not be smaller than neqn  ',info,io)
      if (nparam.le.(neq-neqn))
     *  call error('Need nparam >number of linear equalities',info,io)
      if (nparam.lt.neq) then
        call error('WARNING: nparam < neq                   ',info,io)
        info=0
      endif
      if (iprint.lt.0.or.iprint.gt.3)     
     *  call error('iprint is not a valid number            ',info,io)
      if (eps.gt.epsmac) goto 10
      call error('eps    should be bigger than epsmac!    ',info,io)
      write(io,9902) epsmac
 10   if(mode.ge.200) then
        lstype=2
        mode=mode-100
      else
        lstype=1
      endif
      if (.not.(mode.eq.100.or.mode.eq.101.or.
     *          mode.eq.110.or.mode.eq.111))
     *  call error('mode   is not properly specified!       ',info,io)
      if (info.eq.0) goto 20
      write(io,9903)
      goto 9000
c
 20   do 30 i=1,nparam
        bli=bl(i)
        bui=bu(i)
        if (bli.le.bui) goto 25
        write(io,9901)
        info=7
 25     if (info.ne.0) goto 9000
        if (bli.lt.(-bigbnd)) bl(i)=-bigbnd
        if (bui.gt.bigbnd)    bu(i)=bigbnd
 30   continue
c
      i=mode-100
      if(i.lt.10) then
        modem=0
      else 
        modem=1
        i=i-10
      endif
      if(i.eq.0) Linfty=.false.
      if(i.eq.1) Linfty=.true.
c
 9000 return
 9901 format(1x,'lower bounds should be smaller than upper bounds',/)
 9902 format(1x,'epsmac = ',e22.14,' which is machine dependent',/)
 9903 format(1x,'Error: Input parameters are not consistent',/)
      end
      subroutine initpt(nparam,nnl,neq,neqn,nclin,nctotl,nrowa,x0,
     *                  bndl,bndu,iw,leniw,x,bl,bu,g,gradg,a,cvec,hess,
     *                  clamda,bj,w,lenw,constr,gradcn)
c
c     FFSQP : generation of a feasible point satisfying
c             simple bounds and linear constraints
c
c     implicit real*8(a-h,o-z)
      integer nparam,nnl,neq,neqn,nclin,nctotl,nrowa,leniw,lenw
      integer iw(leniw)
      double  precision x0(nparam),bndl(nparam),bndu(nparam),x(nparam),
     *        bl(1),bu(1),g(1),gradg(nparam,1),
     *        a(nrowa,1),cvec(nparam),hess(nparam,nparam),
     *        clamda(1),bj(1),w(lenw)
c     double  precision x0(nparam),bndl(nparam),bndu(nparam),x(nparam),
c    *        bl(nctotl),bu(nctotl),g(nclin),gradg(nparam,nclin),
c    *        a(nrowa,nparam),cvec(nparam),hess(nparam,nparam),
c    *        clamda(nctotl+nparam),bj(nclin),w(lenw)
      external constr,gradcn
c
c     bj(1) is equivalent to bl(nparam+3)
c
      integer io,iprint,ipspan,ipyes,info,ipd,idum,idum2,maxit,
     *        nnineq,M,id2,id3,id4,id5,id6
      double  precision epsmac,rteps,udelta,valnom,big,tolfea,objeps,
     *        objrep,gLgeps
      logical xisnew
      common  /fsqpp1/nnineq,M,id2,id3,id4,id5,id6
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,idum,idum2,
     *        /fsqpp3/epsmac,rteps,udelta,valnom,
     *        /fsqpq1/big,tolfea,/fsqpq2/maxit
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
      integer i,j,infoql,mnn
      double  precision x0i
c
      info=1
      do 10 i=1,nclin
        valnom=g(i)
        j=i+nnl
        if(j.le.nnineq) call gradcn(nparam,j,x0,gradg(1,i),constr)
        if(j.gt.nnineq) 
     *    call gradcn(nparam,j+neqn,x0,gradg(1,i),constr)
 10   continue
      do 20 i=1,nparam
        x0i=x0(i)
        bl(i)=bndl(i)-x0i
        bu(i)=bndu(i)-x0i
        cvec(i)=0.d0
 20   continue
      do 30 i=nclin,1,-1
 30     bj(nclin-i+1)=-g(i)
      do 60 i=nclin,1,-1
        do 50 j=1,nparam
 50       a(nclin-i+1,j)=-gradg(j,i)
 60   continue
      call diagnl(nparam,1.d0,hess)
      call nullvc(nparam,x)
C
      mnn=nrowa+2*nparam
      iw(1)=1
      call QL0001(nclin,neq-neqn,nrowa,nparam,nparam,mnn,hess,cvec,A,
     *            bj,bL,bU,X,clamda,io,infoql,0,w,lenw,iw,leniw)
      if(infoql.ne.0) goto 90
      do 70 i=1,nparam
 70     x0(i)=x0(i)+x(i)
      xisnew=.true.
      do 80 i=1,nclin
        j=i+nnl
        if(j.le.nnineq) call constr(nparam,j,x0,g(i))
        if(j.gt.nnineq) call constr(nparam,j+neqn,x0,g(i))
 80   continue
      info=0
 90   if(info.eq.1.and.iprint.ne.0) write(io,1000)
 1000 format(1x,'Error: No feasible point is found for the',
     *                 ' linear constraints',/)
      return
      end
      subroutine dir(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nctotl,
     *               nrowa,feasb,steps,epskt,epseqn,sktnom,scvneq,Ck,
     *               d0nm,grdftd,xl,xu,indxob,indxcn,iact,iskp,iskip,
     *               istore,iw,leniw,x,di,d,g,gradg,f,fmax,fM,fMp,psf,
     *               gradf,grdpsf,penp,a,bl,bu,clamda,cllamd,cvec,bj,
     *               hess,hess1,w,lenw,backup,signeq,obj,constr)
c
c     FFSQP : computation of a search direction
c
c     implicit real*8(a-h,o-z)
      integer nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nctotl,nrowa,
     *        iskp,leniw,lenw
      integer indxob(1),indxcn(1),iact(1),iskip(1),
     *        istore(1),iw(leniw)
c     integer indxob(nob),indxcn(ncnstr),iact(nob+nineqn+neqn),iskip(1),
c    *        istore(nineqn+nob+neqn),iw(leniw)
      double  precision steps,epskt,epseqn,sktnom,Ck,d0nm,grdftd,
     *        fmax,fM,fMp,psf,scvneq
      double  precision xl(nparam),xu(nparam),x(nparam+1),di(nparam+1),
     *        d(nparam+1),g(1),gradg(nparam,1),f(1),
     *        gradf(nparam,1),grdpsf(nparam),penp(1),
     *        a(nrowa,nparam+1),bl(1),bu(1),clamda(1),cllamd(1),
     *        cvec(nparam+1),bj(nrowa),hess(nparam,nparam),
     *        hess1(nparam+1,nparam+1),w(lenw),
     *        backup(1),signeq(1)
c     double  precision xl(nparam),xu(nparam),x(nparam+1),di(nparam+1),
c    *        d(nparam+1),g(ncnstr),gradg(nparam,ncnstr),f(nob),
c    *        gradf(nparam,nob), 
c    *        grdpsf(nparam),penp(neqn),a(nrowa,nparam+1),bl(nctotl),
c    *        bu(nctotl),clamda(nctotl+nparam+1),cllamd(nctotl),
c    *        cvec(nparam+1),bj(nrowa),hess(nparam,nparam),
c    *        hess1(nparam+1,nparam+1),w(lenw),
c    *        backup(nob+ncnstr),signeq(neqn)
      external obj,constr
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,ipspan,ipyes,info,
     *        ipd,iter,nstop,initvl,lstype
      double  precision epsmac,rteps,udelta,valnom,bigbnd,tolfea,
     *        objeps,objrep,gLgeps
      logical dlfeas,local,update,first,lqpsl,ld0,rolis1,d0is0,
     *        xisnew
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,iter,initvl,
     *        /fsqpp3/epsmac,rteps,udelta,valnom
     *        /fsqpp4/rolis1,d0is0,
     *        /fsqpq1/bigbnd,tolfea,
     *        /fsqplo/dlfeas,local,update,first,  
     *        /fsqpqp/lqpsl,ld0
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
c     bj(1) is equivalent to bl(nparam+3)
c
      integer i,j,k,kk,ncg,ncf,nqprm0,nclin0,nctot0,infoqp,nqprm1,ncl,
     *        nclin1,nctot1,ncc,nff,nrowa0,nrowa1,ninq,nobb,nobbL,nncn,
     *        max0
      double  precision fmxl,vv,dx,dmx,dnm1,dnm,v0,v1,vk,temp1,temp2,
     *        theta,rhol,rhog,rho,grdfd0,grdfd1,dummy,grdgd0,grdgd1,
     *        thrshd,sign,scaprd,slope,lfuscp,dsqrt,dmin1,dmax1,dabs,
     *        adummy(1),dnmtil
      logical ltem1,ltem2,needd1
c
      ncg=0
      ncf=0
      iskp=0
      ncl=nnineq-nineqn
      local=.false.
      update=.false.
      lqpsl=.false.
      thrshd=tolfea
      needd1=.true.
      rolis1=.false.
c
      if(nobL.le.1) goto 10
        nqprm0=nparam+1
        nclin0=ncnstr+nobL
      goto 20
 10     nqprm0=nparam
        nclin0=ncnstr
 20   nctot0=nqprm0+nclin0
      vv=0.d0
      nrowa0=max0(nclin0,1)
      do 25 i=1,ncnstr
        if(feasb) then
          if(i.gt.nineqn.and.i.le.nnineq) iskip(nnineq+2-i)=i
          iw(i)=i
        else if(.not.feasb) then
          if(i.le.ncl) iskip(ncl+2-i)=nineqn+i
          if(i.le.ncl) iw(i)=nineqn+i
          if(i.gt.ncl) iw(i)=nineqn+neqn+i
        endif
 25   continue
      do 27 i=1,nob
 27     iw(ncnstr+i)=i
      ld0=.true.
      call nullvc(nparam,cvec)
      d0is0=.false.
      call dqp(nparam,nqprm0,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nclin0,
     *         nctot0,nrowa0,infoqp,iw,leniw,x,di,xl,xu,feasb,f,fmax,
     *         gradf,grdpsf,g,gradg,a,cvec,bl,bu,clamda,cllamd,bj,
     *         hess,hess1,di,w,lenw,vv,0)
      ld0=.false.
      if(infoqp.eq.0) goto 30
        info=5
        if(.not.feasb) info=2
        nstop=0
        goto 9000
c
c    reorder indexes of constraints and objectives
c
 30   if(nn.le.1) goto 45  
      j=1
      k=nn
      do 40 i=nn,1,-1
        if(lfuscp(cllamd(nqprm0+indxcn(i)),thrshd).ne.0) then
          iact(j)=indxcn(i)
          j=j+1
        else
          iact(k)=indxcn(i)
          k=k-1
        endif
 40   continue
 45   if(nobL.le.1) goto 60
      j=nn+1
      k=nn+nob
      do 50 i=nob,1,-1
        kk=nqprm0+ncnstr
        ltem1=lfuscp(cllamd(kk+i),thrshd).ne.0
        ltem2=nobL.ne.nob.and.(lfuscp(cllamd(kk+i+nob),thrshd).ne.0)
        if(ltem1.or.ltem2) then
          iact(j)=i
          j=j+1
        else
          iact(k)=i
          k=k-1
        endif
 50   continue
c
 60   vv=f(iact(nn+1))
      d0nm=dsqrt(scaprd(nparam,di,di))
      if(.not.first.or.nclin0.ne.0) goto 110
        dx=dsqrt(scaprd(nparam,x,x))
        dmx=dmax1(dx,1.d0)
        if(d0nm.le.dmx) goto 110
        do 100 i=1,nparam
 100      di(i)=di(i)*dmx/d0nm
        d0nm=dmx
 110  call matrvc(nparam,nparam,nparam,nparam,hess,di,w)
      if(nn.eq.0) grdftd=-scaprd(nparam,w,di)
      sktnom=dsqrt(scaprd(nparam,w,w))
      if(gLgeps.gt.0.d0.and.sktnom.le.gLgeps) goto 115
      if(d0nm.gt.epskt) goto 120
 115  if(neqn.ne.0.and.scvneq.gt.epseqn) goto 120
        nstop=0
        if(.not.feasb) info=2
        if(iprint.lt.3.or.ipyes.gt.0) goto 9000
        if(nobL.le.1) nff=1
        if(nobL.gt.1) nff=2
        call sbout1(io,nparam,'multipliers  for  x ',dummy,cllamd,2,2)
        if(ncnstr.ne.0) call sbout1(io,ncnstr,'             for  g ',
     *                              dummy,cllamd(nparam+nff),2,2)
        if(nobL.gt.1) call sbout1(io,nob,'             for  f ',
     *                            dummy,cllamd(nparam+nff+ncnstr),2,2)
        goto 9000
 120  if(iprint.lt.3.or.ipyes.gt.0) goto 125
        call sbout1(io,nparam,'d0                  ',dummy,di,2,2)
        call sbout1(io,0,'d0norm              ',d0nm,adummy,1,2)
        call sbout1(io,0,'ktnorm              ',sktnom,adummy,1,2)
 125  temp1=dmin1(0.5d0*epskt,0.1d-1*rteps)
      if(neqn.eq.0.or.scvneq.le.epseqn.or.d0nm.gt.temp1) goto 127
        d0is0=.true.
        goto 9000
c
c     single objective without nonlinear constraints requires
c     no d1 and dtilde; multi-objectives without nonlinear 
c     constraints requires no d1
c
 127  call nullvc(nparam,w)
      if(nn.ne.0) grdftd=slope(nob,nobL,neqn,nparam,feasb,f,gradf,
     *                         grdpsf,di,w,fmax,dummy,0)
      if(nn.eq.0.and.nobL.le.1) goto 1130 
      if(nn.ne.0) goto 130
        dnm=d0nm
        rho=0.d0
        rhog=0.d0
        goto 310
c
c     compute modified first order direction d1
c
 130  if (mode.eq.1) then
        needd1=.false.
        vk=dmin1(Ck*d0nm*d0nm,d0nm)
        do 131 i=1,nn
          grdgd0=scaprd(nparam,gradg(1,indxcn(i)),di)
          temp1=vk+g(indxcn(i))+grdgd0
          if(temp1.le.0.d0) goto 131
          needd1=.true.
          goto 132
 131    continue
      endif
 132  if (needd1) then
        nqprm1=nparam+1
        if(mode.eq.0) nclin1=ncnstr+max0(nobL,1)
        if(mode.eq.1) nclin1=ncnstr
        nctot1=nqprm1+nclin1
        nrowa1=max0(nclin1,1)
        ninq=nnineq
        call di1(nparam,nqprm1,nob,nobL,nineqn,neq,neqn,ncnstr,nclin1,
     *           nctot1,nrowa1,infoqp,mode,iw,leniw,x,di,xl,xu,f,fmax,
     *           gradf,grdpsf,g,gradg,cvec,a,bl,bu,clamda,bj,hess1,d,
     *           w,lenw)
        if(infoqp.eq.0) goto 140
          info=6
          if(.not.feasb) info=2
          nstop=0
          goto 9000
 140    dnm1=dsqrt(scaprd(nparam,d,d))
        if(iprint.lt.3.or.ipyes.gt.0) goto 145
          call sbout1(io,nparam,'d1                  ',dummy,d,2,2)
          call sbout1(io,0,'d1norm              ',dnm1,adummy,1,2)
      else
        dnm1=0.d0
        do 141 i=1,nparam
 141      d(i)=0.d0
      endif
 145  if(mode.eq.1) goto 150
        v0=d0nm**2.1
        v1=dmax1(dble(0.5),dble(dnm1**2.5))
        rho=v0/(v0+v1)
        rhog=rho
      goto 250
 150    rhol=0.d0
        if(.not.needd1) goto 210
        do 200 i=1,nn
          grdgd0=scaprd(nparam,gradg(1,indxcn(i)),di)
          grdgd1=scaprd(nparam,gradg(1,indxcn(i)),d)
          temp1=vk+g(indxcn(i))+grdgd0
          temp2=grdgd1-grdgd0
          if(temp1.le.0.d0) goto 200
          if(temp2.ge.0.d0) goto 190
          rhol=dmax1(rhol,-temp1/temp2)
          if(rhol.lt.1.d0) goto 200
 190        rhol=1.0d0
            rolis1=.true.
            goto 210
 200    continue
 210    theta=0.2d0
        if(rhol.ne.0.d0) goto 220
c
c       to check if rhol is reset
c
          rhog=0.d0
          rho=0.d0
          dnm=d0nm
        goto 310
 220    if(nobL.gt.1) goto 230
          grdfd0=grdftd
          if(nob.eq.1) grdfd1=scaprd(nparam,gradf(1,1),d)
          grdfd1=grdfd1-scaprd(nparam,grdpsf,d)
          temp1=grdfd1-grdfd0
          if(temp1.le.0.d0) then
            rhog=rhol
          else
            rhog=dmin1(rhol,(theta-1.d0)*grdfd0/temp1)
          endif
        goto 240
 230      rhog=slope(nob,nobL,neqn,nparam,feasb,f,gradf(1,1),grdpsf,
     *               di,d,fmax,theta,mode)
          rhog=dmin1(rhol,rhog)
 240    rho=rhog
        if (steps.eq.1.d0.and.rhol.lt.0.5d0) rho=rhol
 250  continue
      do 300 i=1,nparam
        if (rho.ne.rhog) cvec(i)=di(i)
        di(i)=(1.d0-rho)*di(i)+rho*d(i)
 300  continue
      dnm=dsqrt(scaprd(nparam,di,di))
      if(iprint.lt.3.or.mode.eq.1.or.nn.eq.0.or.ipyes.gt.0) goto 310
        call sbout1(io,0,'rho                 ',rho,adummy,1,2)
        call sbout1(io,nparam,'d                   ',dummy,di,2,2)
        call sbout1(io,0,'dnorm               ',dnm,adummy,1,2)
 310  continue
 320  do 400 i=1,nob
 400    bl(i)=f(i)
      if (rho.eq.1.d0) goto 510
      if(nn.eq.0.or.iprint.ne.3.or.mode.eq.0.or.ipyes.gt.0) goto 410
        call sbout1(io,0,'Ck                  ',Ck,adummy,1,2)
        call sbout1(io,0,'rhol                ',rho,adummy,1,2)
        call sbout1(io,nparam,'dl                  ',dummy,di,2,2)
        call sbout1(io,0,'dlnorm              ',dnm,adummy,1,2)
 410  if(mode.eq.0) goto 510
        local=.true.
        call step(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,ncg,ncf,
     *            indxob,indxcn,iact,iskp,iskip,istore,feasb,grdftd,
     *            f,fmax,fM,fMp,psf,penp,steps,scvneq,bu,x,di,d,g,w,
     *            backup,signeq,obj,constr)
        if(update) goto 9000
        local=.false.
        if(rho.eq.rhog.or.nn.eq.0) goto 510
        do 500 i=1,nparam
 500      di(i)=(1-rhog)*cvec(i)+rhog*d(i)
        dnm=dsqrt(scaprd(nparam,di,di))
 510  if (nn.eq.0.or.iprint.lt.3.or.mode.eq.0.or.ipyes.gt.0) goto 520
        call sbout1(io,0,'rhog                ',rhog,adummy,1,2)
        call sbout1(io,nparam,'dg                  ',dummy,di,2,2)
        call sbout1(io,0,'dgnorm              ',dnm,adummy,1,2)
 520  if(rho.ne.0.d0) grdftd=slope(nob,nobL,neqn,nparam,feasb,bl,
     *                             gradf,grdpsf,di,d,fmax,theta,0)
      if(mode.eq.1.and.rho.eq.rhog) goto 610
      do 600 i=1,nparam
 600    bu(i)=x(i)+di(i)
      xisnew=.true.
 610  if(rho.ne.rhog) ncg=0
      ncc=ncg+1
      fmxl=-bigbnd
      ninq=ncg
      nncn=ncg
      j=0
c
c     iskip(1) --- iskip(iskp) store the indexes of linear inequality
c     constraints that are not to be used to compute d~
c     iskip(nnineq-nineqn+1) --- iskip(nnineq-ncn+1-iskp) store those 
c     that are to be used to compute d~
c
      do 700 i=ncc,ncnstr
        if(i.le.nn) then
          kk=iact(i)
        else
          kk=indxcn(i)
        endif
        if(kk.le.nineqn.or.kk.gt.nnineq) goto 615
          iskip(ncl+1-j)=kk
          j=j+1
 615    if(kk.gt.nnineq) goto 617
        temp1=-0.2d0*(dnm*dsqrt(scaprd(nparam,gradg(1,kk),gradg(1,kk))))
        temp2=cllamd(nqprm0+kk)
        if(temp2.eq.0.d0.and.g(kk).lt.temp1) goto 620
 617      ninq=ninq+1
          iw(ninq)=kk
          if(feasb.and.kk.le.nineqn) istore(kk)=1
          call constr(nparam,kk,bu,g(kk))
          if(.not.feasb.or.feasb.and.kk.gt.(nnineq+neqn)) goto 700
          if(kk.le.nineqn) nncn=ninq
          fmxl=dmax1(fmxl,g(kk))
          if(.not.feasb) goto 618
          if(kk.le.nineqn.or.kk.gt.nnineq.and.kk.le.(nnineq+neqn))
     *       ncallg=ncallg+1
 618      if(dabs(fmxl).gt.bigbnd) goto 1130
        goto 700
 620      if(kk.le.nineqn) goto 700
          iskp=iskp+1
          iskip(iskp)=kk
          j=j-1
 700  continue
      if(neqn.ne.0) call resign(nparam,neqn,psf,grdpsf,penp,g(nnineq+1),
     *                          gradg(1,nnineq+1),signeq,10,20)
      ninq=ninq-neq
      if(.not.feasb) ninq=ninq+neqn
      if(ncg.eq.0) goto 810
      do 800 i=1,ncg
        iw(i)=iact(i)
        if(iact(i).le.nineqn) istore(iact(i))=1
        fmxl=dmax1(fmxl,g(iact(i)))
        if(dabs(fmxl).gt.bigbnd) goto 1130
 800  continue
 810  if(nobL.gt.1) goto 820
        iw(1+ninq+neq)=1
        nobb=nob
        goto 1110
 820  if(rho.ne.rhog) ncf=0
      nff=ncf+1
      nobb=ncf
      sign=1.d0
      fmxl=-bigbnd
      if(cllamd(nqprm0+ncnstr+iact(nn+1)).lt.0.d0) sign=-1.d0
      do 1000 i=nff,nob
        kk=iact(nn+i)
        if(.not.feasb) kk=iact(i)
        if(feasb) k=nn+1
        if(.not.feasb) k=1
        do 900 j=1,nparam
 900      w(j)=sign*gradf(j,iact(k))-gradf(j,kk)
        temp1=dabs(f(kk)-sign*vv)
        temp2=dnm*dsqrt(scaprd(nparam,w,w))
        if(temp1.eq.0.d0.or.temp2.eq.0.d0) goto 910
        temp1=temp1/temp2
        temp2=cllamd(nqprm0+ncnstr+kk)
        if(temp2.eq.0.d0.and.temp1.gt.0.2d0) goto 1000
 910    nobb=nobb+1
        iw(nobb+ninq+neq)=kk
        if(feasb)      istore(nineqn+kk)=1
        if(.not.feasb) istore(kk)=1
        if(.not.feasb) goto 920
          call obj(nparam,kk,bu,f(kk))
          ncallf=ncallf+1
          if(nobL.ne.nob) fmxl=dmax1(fmxl,-f(kk))
        goto 930
 920      call constr(nparam,indxob(kk),bu,f(kk))
          ncallg=ncallg+1
 930    fmxl=dmax1(fmxl,f(kk))
        if(dabs(fmxl).gt.bigbnd) goto 1130
 1000 continue
      if(ncf.eq.0) goto 1110
      do 1100 i=1,ncf
        iw(ninq+neq+i)=iact(i+nn)
        istore(nineqn+iact(i+nn))=1
        fmxl=dmax1(fmxl,f(iact(i+nn)))
        if(nobL.ne.nob) fmxl=dmax1(fmxl,-f(iact(i+nn)))
        if(dabs(fmxl).gt.bigbnd) goto 1130
 1100 continue
 1110 call matrvc(nparam,nparam,nparam,nparam,hess,di,cvec)
      vv=-dmin1(0.01d0*dnm,dnm**2.5)
c
c     compute a correction dtilde to d=(1-rho)d0+rho*d1
c
      if(nobL.ne.nob) nobbL=2*nobb
      if(nobL.eq.nob) nobbL=nobb
      if(nobbL.le.1) goto 1115
        nqprm0=nparam+1
        nclin0=ninq+neq+nobbL
      goto 1117
 1115   nqprm0=nparam
        nclin0=ninq+neq
 1117 nctot0=nqprm0+nclin0
      nrowa0=max0(nclin0,1)
      i=ninq+neq
      call dqp(nparam,nqprm0,nobb,nobbL,nncn,neq,neqn,nn,i,nclin0,
     *         nctot0,nrowa0,infoqp,iw,leniw,x,di,xl,xu,feasb,f,fmxl,
     *         gradf,grdpsf,g,gradg,a,cvec,bl,bu,clamda,cllamd,bj,
     *         hess,hess1,d,w,lenw,vv,1)
      if(infoqp.ne.0) goto 1130
      dnmtil=dsqrt(scaprd(nparam,d,d))
      if(dnmtil.gt.dnm) goto 1130
      if(dnmtil.eq.0.d0) goto 1119
        do 1118 i=1,nineqn+nob
 1118     istore(i)=0
 1119 if(iprint.lt.3.or.ipyes.gt.0) goto 9000
        call sbout1(io,nparam,'dtilde              ',dummy,d,2,2)
        call sbout1(io,0,'dtnorm              ',dnmtil,adummy,1,2)
        goto 9000
c
 1130 do 1200 i=1,nparam
 1200   d(i)=0.d0
      dnmtil=0.d0
 9000 return
      end
c
      subroutine dqp(nparam,nqpram,nob,nobL,nineqn,neq,neqn,nn,ncnstr,
     *               nclin,nctotl,nrowa,infoqp,iw,leniw,x0,di,xl,xu,
     *               feasb,f,fmax,gradf,grdpsf,g,gradg,a,cvec,bl,bu,
     *               clamda,cllamd,bj,hess,hess1,x,w,lenw,vv,job)
c     implicit double precision(a-h,o-z)
      integer nparam,nqpram,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nclin,
     *        nctotl,nrowa,infoqp,leniw,lenw,job
      integer iw(leniw)
      double  precision fmax,vv
      double  precision x0(nparam),di(1),xl(nparam),xu(nparam),
     *        f(1),gradf(nparam,1),grdpsf(nparam),g(1),
     *        gradg(nparam,1),
     *        a(nrowa,1),cvec(1),bl(1),bu(1),clamda(1),
     *        cllamd(1),bj(1),hess(nparam,nparam),
     *        hess1(nparam+1,nparam+1),x(1),w(lenw)
c     double  precision x0(nparam),di(nqpram),xl(nparam),xu(nparam),
c    *        f(nob),gradf(nparam,nob),grdpsf(nparam),g(ncnstr),
c    *        gradg(nparam,ncnstr),
c    *        a(nrowa,nqpram),cvec(nqpram),bl(nctotl),bu(nctotl),
c    *        clamda(nctotl+nqpram),cllamd(nctotl),bj(nrowa),
c    *        hess(nparam,nparam),hess1(nparam+1,nparam+1),
c    *        x(nqpram),w(lenw)
      logical feasb
c
      integer io,idum1,idum2,idum3,idum4,idum5,idum6,idum7
      double  precision bigbnd,dummy,epsmac,rteps,dummy1,dummy2
      common  /fsqpp2/io,idum1,idum2,idum3,idum4,idum5,idum6,idum7,
     *        /fsqpp3/epsmac,rteps,dummy1,dummy2,
     *        /fsqpq1/bigbnd,dummy
c
c     bj(1) is equivalent to bl(nparam+3)
c
c     job=0 : compute d0; job=1 : compute  d~
c
      integer i,ii,j,iout,mnn,nqnp
      double  precision x0i,xdi
c
      iout=io
      do 100 i=1,nparam
        x0i=x0(i)
        if(job.eq.1) xdi=di(i)
        if(job.eq.0) xdi=0.d0
        bl(i)=xl(i)-x0i-xdi
        bu(i)=xu(i)-x0i-xdi
        cvec(i)=cvec(i)-grdpsf(i)
 100  continue
      if(nobL.le.1) goto 110
        bl(nqpram)=-bigbnd
        bu(nqpram)=bigbnd
 110  ii=ncnstr-nn
c
c     constraints are assigned to a in reverse order
c
      do 300 i=1,ncnstr
        x0i=vv
        if(i.le.(neq-neqn).or.(i.gt.neq.and.i.le.(ncnstr-nineqn)))
     *    x0i=0.d0
        if(.not.feasb) x0i=0.d0
        bj(i)=x0i-g(iw(ncnstr+1-i))
        do 200 j=1,nparam
 200      a(i,j)=-gradg(j,iw(ncnstr+1-i))
        if(nobL.gt.1) a(i,nqpram)=0.d0
 300  continue
      if(nobL.le.1) goto 510
      do 500 i=1,nob
        ii=ncnstr+i
        bj(ii)=fmax-f(iw(ii))
        if(nobL.gt.nob) bj(ii+nob)=fmax+f(iw(ii))
        do 400 j=1,nparam
          a(ii,j)=-gradf(j,iw(ii))
          if(nobL.gt.nob) a(ii+nob,j)=gradf(j,iw(ii))
 400    continue
        a(ii,nqpram)=1.d0
        if(nobL.gt.nob) a(ii+nob,nqpram)=1.d0
 500  continue
      cvec(nqpram)=1.d0
      goto 610
 510  if(nob.eq.0) goto 610
      do 600 i=1,nparam
 600    cvec(i)=cvec(i)+gradf(i,1)
 610  call matrcp(nparam,hess,nparam+1,hess1)
      call nullvc(nqpram,x)
c
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
c  The following modification is done inside QP0001
c  for the ease of interfacing with QPSOL
c
c     if(hess1(nqpram,nqpram).lt.qleps) hess1(nqpram,nqpram)=qleps
C
      iw(1)=1
      mnn=nclin+2*nqpram
      call QL0001(nclin,neq-neqn,nrowa,nqpram,nparam+1,mnn,hess1,cvec,A,
     *            bj,bL,bU,X,clamda,iout,infoqp,0,w,lenw,iw,leniw)
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      if(infoqp.ne.0.or.job.eq.1) goto 9000
      do 700 i=1,nqpram
        ii=nclin+i
        if(clamda(ii).eq.0.d0.and.clamda(ii+nqpram).eq.0.d0) then
          goto 700
        else if(clamda(ii).ne.0.d0) then
          clamda(ii)=-clamda(ii)
        else 
          clamda(ii)=clamda(ii+nqpram)
        endif
 700  continue
      nqnp=nqpram+ncnstr
      do 800 i=1,nctotl
        if(i.le.nqpram) then
          ii=nclin+i
        else if(i.gt.nqpram.and.i.le.nqnp) then
          ii=nqnp+1-i
        else if(i.gt.nqnp) then
          ii=i-nqpram
        endif
        cllamd(i)=clamda(ii)
 800  continue
      if(nobL.eq.nob) goto 9000
      do 900 i=1,nob
        ii=i+nqpram+ncnstr
        cllamd(ii)=cllamd(ii)-cllamd(ii+nob)
 900  continue
 9000 return
      end
c
      subroutine di1(nparam,nqpram,nob,nobL,nineqn,neq,neqn,ncnstr,
     *               nclin,nctotl,nrowa,infoqp,mode,iw,leniw,x0,d0,
     *               xl,xu,f,fmax,gradf,grdpsf,g,gradg,cvec,a,bl,bu,
     *               clamda,bj,hess1,x,w,lenw)
c     implicit real*8(a-h,o-z)
      integer nparam,nqpram,nob,nobL,nineqn,neq,neqn,ncnstr,nclin,
     *        nctotl,nrowa,infoqp,mode,leniw,lenw,iw(leniw)
      double  precision fmax
      double  precision x0(nparam),d0(nparam),xl(nparam),xu(nparam),
     *        f(1),gradf(nparam,1),grdpsf(nparam),g(1),
     *        gradg(nparam,1),cvec(1),a(nrowa,1),
     *        bl(1),bu(1),clamda(1),bj(1),
     *        hess1(nparam+1,nparam+1),x(1),w(lenw)
c     double  precision x0(nparam),d0(nparam),xl(nparam),xu(nparam),
c    *        f(nob),gradf(nparam,nob+1),grdpsf(nparam),g(ncnstr),
c    *        gradg(nparam,ncnstr),cvec(nqpram),a(nrowa,nqpram),
c    *        bl(nctotl),bu(nctotl),clamda(nctotl+nqpram),bj(nrowa),
c    *        hess1(nparam+1,nparam+1),x(nqpram),w(lenw)
c
      integer io,idum1,idum2,idum3,idum4,idum5,idum6,idum7
      double  precision epsmac,rteps,dumm1,dumm2,bigbnd,dummy
      common  /fsqpp2/io,idum1,idum2,idum3,idum4,idum5,idum6,idum7,
     *        /fsqpp3/epsmac,rteps,dumm1,dumm2,
     *        /fsqpq1/bigbnd,dummy
c
c     bj(1) is equivalent to bl(nparam+3)
c
      integer i,ii,iout,j,mnn
      double  precision x0i,eta
c
      iout=io
      if(mode.eq.0) eta=0.1d0
      if(mode.eq.1) eta=3.d0
      do 100 i=1,nparam
        x0i=x0(i)
        bl(i)=xl(i)-x0i
        bu(i)=xu(i)-x0i
        if(mode.eq.0) cvec(i)=-eta*d0(i)
        if(mode.eq.1) cvec(i)=0.d0
 100  continue
      bl(nqpram)=-bigbnd
      bu(nqpram)=bigbnd
      cvec(nqpram)=1.d0
      ii=ncnstr-nineqn
      do 400 i=1,ncnstr
        bj(i)=-g(ncnstr+1-i)
        do 300 j=1,nparam
 300      a(i,j)=-gradg(j,ncnstr+1-i)
        a(i,nqpram)=0.d0
        if((i.gt.(neq-neqn).and.i.le.neq).or.i.gt.ii) a(i,nqpram)=1.d0
 400  continue
      if(mode.eq.1) goto 610
      i=0
 450  continue
        i=i+1
        ii=ncnstr+i
        if(nob.eq.0) bj(ii)=fmax
        if(nob.gt.0) bj(ii)=fmax-f(i)
        do 500 j=1,nparam
          if(nob.eq.0) a(ii,j)=grdpsf(j)
          if(nob.gt.0) a(ii,j)=-gradf(j,i)+grdpsf(j)
          if(nobL.gt.nob) a(ii+nob,j)=gradf(j,i)+grdpsf(j)
 500    continue
        a(ii,nqpram)=1.d0
        if(nobL.gt.nob) a(ii+nob,nqpram)=1.d0
        if(i.lt.nob) goto 450
 610  call diagnl(nqpram,eta,hess1)
      call nullvc(nqpram,x)
      hess1(nqpram,nqpram)=0.d0
c
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  The following modification is done inside QP0001
c  for the ease of interfacing with QPSOL
c
c     hess1(nqpram,nqpram)=qleps
C
      mnn=nclin+2*nqpram
      iw(1)=1
      call QL0001(nclin,neq-neqn,nrowa,nqpram,nparam+1,mnn,hess1,cvec,A,
     *             bj,bL,bU,X,clamda,iout,infoqp,0,w,lenw,iw,leniw)
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end

      subroutine step(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,ncg,
     *                ncf,indxob,indxcn,iact,iskp,iskip,istore,feasb,
     *                grdftd,f,fmax,fM,fMp,psf,penp,steps,scvneq,xnew,
     *                x,di,d,g,w,backup,signeq,obj,constr)
c
c     FFSQP : Armijo or nonmonotone line search, with
c             some ad hoc strategies to decrease the number
c             of function evaluation as much as possible
c
c     implicit real*8(a-h,o-z)
      integer nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,ncg,ncf,iskp
      integer indxob(1),indxcn(1),iact(1),iskip(1),
     *        istore(1)
c     integer indxob(nob),indxcn(ncnstr),iact(nn+nob),iskip(4),
c    *        istore(nineqn+nob)
      double  precision grdftd,fmax,fM,fMp,steps,scvneq,psf
      double  precision xnew(nparam),x(nparam),di(nparam),d(nparam),
     *        f(1),penp(1),g(1),w(1),backup(1),
     *        signeq(1)
c     double  precision xnew(nparam),x(nparam),di(nparam),d(nparam),
c    *        f(nob),penp(neqn),g(ncnstr),w(1),backup(nob+ncnstr),
c    *        signeq(neqn)
      external obj,constr
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,ipspan,ipyes,info,
     *        idum1,idum2,idum3,nstop,lstype
      double  precision epsmac,bigbnd,tolfea,dum1,dum2,dum3,fii,
     *        objeps,objrep,gLgeps
      logical lqpsl,ldummy,dlfeas,local,update,ldumm2,xisnew
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,idum1,idum2,idum3,
     *        /fsqpp3/epsmac,dum1,dum2,dum3,
     *        /fsqpq1/bigbnd,tolfea,
     *        /fsqplo/dlfeas,local,update,ldumm2,
     *        /fsqpqp/lqpsl,ldummy
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
      integer i,ii,ij,itry,ikeep,j,job,nlin,mnm,ntwo
      double  precision prod1,prod,dummy,tolfe,dmax1,ostep,
     *                  adummy(1),temp,fmaxl
      logical ltem1,ltem2,reform,fbind,cdone,fdone,eqdone
      data ntwo/2/
c
c    The above data statement is to fool some compilers
c    that do not like a hard number as the index of w in
c    "call resign ..."
c
      nlin=nnineq-nineqn
      ii=1
      itry=1
      steps=1.d0
      ostep=steps
      fbind=.false.
      cdone=.false.
      fdone=.false.
      eqdone=.false.
      if(local) dlfeas=.false.
      ikeep=nlin-iskp
      prod1=0.1d0*grdftd
      tolfe=0.d0
      if(lqpsl) tolfe=tolfea
      if(iprint.ge.3.and.ipyes.eq.0)
     *  call sbout1(io,0,'directional deriv.  ',grdftd,adummy,1,2)
 
      w(1)=fM
 100  continue
        reform=.true.
        if(iprint.ge.3.and.ipyes.eq.0) 
     *    write(io,9901) itry
        prod=prod1*steps
        if(.not.feasb.or.nobL.gt.1) prod=prod+tolfe
        do 200 i=1,nparam
          if(local)      xnew(i)=x(i)+steps*di(i)
          if(.not.local) xnew(i)=x(i)+steps*di(i)+d(i)*steps**2
 200    continue
        xisnew=.true.
        if(iprint.lt.3.or.ipyes.gt.0) goto 205
          call sbout1(io,0,'trial step          ',steps,adummy,1,2)
          call sbout1(io,nparam,'trial point         ',
     *                dummy,xnew,2,2)
 205    if(iskp.eq.0) goto 209
          ostep=steps
          do 207 i=ii,iskp
            ij=iskip(i)
            call constr(nparam,ij,xnew,g(ij))
            if(iprint.lt.3.or.ipyes.gt.0) goto 206
              if(i.eq.1) write(io,9900) ij,g(ij)
              if(i.ne.1) write(io,9902) ij,g(ij)
 206        if(g(ij).le.tolfe) goto 207
            ii=i
            goto 1120
 207      continue
          iskp=0
 209    if(nn.eq.0) goto 310
        if(.not.local.and.fbind) goto 315
 210    continue
        do 300 i=1,nn
          ncg=i
          ii=iact(i)
          ij=nnineq+neqn
          if(ii.le.nineqn.and.istore(ii).eq.1) goto 215
          if(ii.gt.nnineq.and.ii.le.ij.and.eqdone) goto 215
            temp=1.d0
            if(ii.gt.nnineq.and.ii.le.ij) temp=signeq(ii-nnineq)
            call constr(nparam,ii,xnew,g(ii))
            g(ii)=g(ii)*temp
            ncallg=ncallg+1
 215      if(iprint.lt.3.or.ipyes.gt.0) goto 220
            if(i.eq.1.and.ikeep.eq.nlin) 
     *        write(io,9900) ii,g(ii)
            if(i.ne.1.or.ikeep.ne.nlin) write(io,9902) ii,g(ii)
 220      if(local.or.g(ii).le.tolfe) goto 230
            call shift(nn,ii,iact)
            goto 1110
 230      if(local.and.g(ii).gt.tolfe) goto 1500
 300    continue
 310    cdone=.true.
        eqdone=.true.
        if(local) dlfeas=.true.
 315    if(fdone) goto 410
        i = 0
        if(nob.gt.0) fmaxl=-bigbnd
 400    continue
          if(i.gt.nob) goto 405
          if(nob.ne.0.and.i.eq.0) i = 1
          ncf=i
          ii=iact(nn+i)
          if(feasb) then
            if(eqdone.or.neqn.eq.0) goto 317
              do 316 j=1,neqn
 316            call constr(nparam,nnineq+j,xnew,g(nnineq+j))
              ncallg=ncallg+neqn
 317        if(neqn.eq.0) goto 318
              if(eqdone)      job=20
              if(.not.eqdone) job=10
              call resign(nparam,neqn,psf,w(ntwo),penp,
     *                    g(nnineq+1),w(ntwo),signeq,job,10)
 318        if(i.eq.0.or.istore(nineqn+ii).eq.1) goto 320
              call obj(nparam,ii,xnew,f(ii))
              ncallf=ncallf+1
 320        if(i.eq.0) fii = 0.d0
            if(i.ne.0) fii = f(ii)
            if(nob.gt.0.and.(i.le.1.and.iprint.ge.3.and.ipyes.eq.0)) 
     *        write(io,9903) ii,fii-psf
            if(nob.eq.0.and.(i.le.1.and.iprint.ge.3.and.ipyes.eq.0)) 
     *        write(io,9904) fii-psf
            if(i.gt.1.and.iprint.ge.3.and.ipyes.eq.0) 
     *        write(io,9902) ii,fii-psf
          else
            if(istore(ii).eq.1) goto 325
              call constr(nparam,indxob(ii),xnew,f(ii))
              ncallg=ncallg+1
 325        if(f(ii).gt.tolfe) reform=.false.
            if(i.eq.1.and.iprint.ge.3.and.ipyes.eq.0) 
     *        write(io,9903) indxob(ii),f(ii)
            if(i.ne.1.and.iprint.ge.3.and.ipyes.eq.0) 
     *        write(io,9902) indxob(ii),f(ii)
            fii=f(ii)
          endif
          fmaxl=dmax1(fmaxl,fii)
          if(nobL.ne.nob) fmaxl=dmax1(fmaxl,-fii)
          if(.not.feasb.and.reform) goto 401
          if(local) goto 340
          if((fii-psf).le.(fMp+prod)) goto 330
            fbind=.true.
            call shift(nob,ii,iact(nn+1))
          goto 1110
 330      if(nobL.eq.nob.or.(-fii-psf).le.(fMp+prod)) goto 401
            fbind=.true.
            call shift(nob,ii,iact(nn+1))
          goto 1110
 340      ltem1=(fii-psf).gt.(fMp+prod)
          ltem2=nobL.ne.nob.and.(-fii-psf).gt.(fMp+prod)
          if(ltem1.or.ltem2) goto 1500
 401      i = i + 1
        goto 400
 405    fbind=.false.
        fdone=.true.
        eqdone=.true.
        if(.not.cdone) goto 210
 410    if(ostep.eq.steps) mnm=ikeep+neq-neqn
        if(ostep.ne.steps) mnm=ncnstr-nn
        do 500 i=1,mnm
          ii=indxcn(i+nn)
          if(ikeep.ne.nlin.and.ostep.eq.steps) then
            if(i.le.ikeep) ii=iskip(nlin+2-i)
            if(i.gt.ikeep) ii=indxcn(nn+i-ikeep+nlin)
          endif
          call constr(nparam,ii,xnew,g(ii))
 500    continue
        scvneq=0.d0
        do 600 i=1,ncnstr
          if(i.gt.nnineq.and.i.le.(nnineq+neqn)) scvneq=scvneq-g(i)
 600      backup(i)=g(i)
        do 700 i=1,nob
 700      backup(i+ncnstr)=f(i)
        if(feasb.or..not.reform) goto 810
          do 800 i=1,nparam
 800        x(i)=xnew(i)
          nstop=0
          goto 1500
 810    if(local) ncg=ncnstr
        if(local) update=.true.
        fM=fmaxl
        fMp=fmaxl-psf
        fmax=fmaxl
        do 1000 i=1,nn
 1000     iact(i)=indxcn(i)
        do 1100 i=1,nob
 1100     iact(nn+i)=i
        goto 1500
c
 1110   cdone=.false.
        fdone=.false.
        eqdone=.false.
        reform=.false.
        if(lstype.eq.2) fbind=.false.
 1120   itry=itry+1
        if(steps.lt.1.d0) goto 1140
        do 1130 i=1,nob+nineqn
 1130     istore(i)=0
 1140   steps=steps*.5d0
        if(steps.lt.epsmac) goto 1150
      goto 100
c
 1150 info=4
      nstop=0
 1500 if(steps.lt.6.d-1) goto 9000
        do 1600 i=1,nob+nineqn
 1600     istore(i)=0
 9000 return
 9900 format(1x,t17,17htrial constraints,t37,i7,t45,e22.14)
 9901 format(1x,t17,12htrial number,t45,i22)
 9902 format(1x,t37,i7,t45,e22.14)
 9903 format(1x,t17,16htrial objectives,t37,i7,t45,e22.14)
 9904 format(1x,t17,25htrial penalized objective,t45,e22.14)
      end

      subroutine hesian(nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,
     *                  nctotl,nfs,nstart,feasb,bigbnd,xnew,x,f,fmax,
     *                  fM,fMp,psf,gradf,grdpsf,penp,g,gm,gradg,indxob,
     *                  indxcn,cllamd,delta,eta,gamma,hess,hd,steps,
     *                  nrst,signeq,span,obj,constr,gradob,gradcn,
     *                  phess,psb,psmu,w,lenw,iw,leniw)
c
c     FFSQP : updating the Hessian matrix using BFGS
c             formula with Powell's modification
c
c     implicit real*8(a-h,o-z)
      integer nparam,nob,nobL,nineqn,neq,neqn,nn,ncnstr,nctotl,nfs,
     *        nstart,indxob(1),indxcn(1),nrst,lenw,leniw,iw(leniw)
c    *        nstart,indxob(nob),indxcn(1),nrst,lenw,leniw,iw(leniw)
      double  precision bigbnd,steps,psf,fmax,fM,fMp,
     *        xnew(nparam),x(nparam),f(1),gradf(nparam,1),
     *        grdpsf(nparam),penp(1),g(1),gm(1),
     *        gradg(nparam,1),cllamd(1),delta(nparam),
     *        eta(nparam),gamma(nparam),hess(nparam,nparam),hd(nparam),
     *        signeq(1),span(1),phess(1),psb(1),psmu(1),w(lenw)
c     double  precision bigbnd,steps,psf,fmax,fM,fMp,
c    *        xnew(nparam),x(nparam),f(nob),gradf(nparam,nob),
c    *        grdpsf(nparam),penp(neqn),g(ncnstr),gm(4*neqn),
c    *        gradg(nparam,ncnstr),cllamd(nctotl),delta(nparam),
c    *        eta(nparam),gamma(nparam),hess(nparam,nparam),hd(nparam),
c    *        signeq(neqn),span(1),phess(neq,neq),psb(neq),
c    *        psmu(neq),w(lenw)
      external obj,constr,gradob,gradcn
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,ipspan,ipyes,info,
     *        ipd,iter,nstop,initvl,lstype
      double  precision epsmac,rteps,udelta,valnom,objeps,objrep,gLgeps
      logical rolis1,d0is0,xisnew
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,iter,initvl,
     *        /fsqpp3/epsmac,rteps,udelta,valnom,
     *        /fsqpp4/rolis1,d0is0,
     *        /fsqpus/objeps,objrep,gLgeps,xisnew
c
      integer ng,i,j,ifail,indexs,np,mnm,iout
      double  precision dhd,gammd,etad,scaprd,dummy,theta,signgj,psfnew
      logical done
c
      if(feasb.and.nstop.ne.0.and.neqn.eq.0) then
c
c       check of gLgeps is just after computing d0!
c
        if(dabs(w(1)-fmax).le.objeps) then
          nstop=0
        else if(dabs(w(1)-fmax).le.objrep*dabs(w(1))) then
          nstop=0
        endif
      endif
      if(nstop.eq.0) goto 810
c
      ipd=0
      done=.false.
      psfnew=0.d0
      call nullvc(nparam,delta)
      call nullvc(nparam,eta)
      if(nobL.gt.1) ng=2
      if(nobL.le.1) ng=1
c
 100  continue
        call nullvc(nparam,gamma)
        if(nobL.gt.1) call matrvc(nparam,nob,nparam,nob,gradf,
     *                     cllamd(nparam+ng+ncnstr),hd)
        if(.not.feasb) goto 120
        if(nineqn.eq.0) goto 110
        call matrvc(nparam,nineqn,nparam,nineqn,gradg,cllamd(nparam+ng),
     *              gamma)
 110    if(neqn.eq.0) goto 120
        call matrvc(nparam,neqn,nparam,neqn,gradg(1,nnineq+1),
     *              cllamd(nparam+nnineq+ng),eta)
 120    do 200 i=1,nparam
          if(nobL.gt.1) then
            if(done) psb(i)=hd(i)+cllamd(i)+gamma(i)
            gamma(i)=gamma(i)+hd(i)-grdpsf(i)+eta(i)
          else if(nobL.eq.1) then
            if(done) psb(i)=gradf(i,1)+cllamd(i)+gamma(i)
            gamma(i)=gamma(i)+gradf(i,1)-grdpsf(i)+eta(i)
          else if(nobL.eq.0) then
            if(done) psb(i)=cllamd(i)+gamma(i)
            gamma(i)=gamma(i)-grdpsf(i)+eta(i)
          endif
          if(.not.done) delta(i)=gamma(i)
 200    continue
        if(done) goto 410
        if(d0is0) goto 405
        if(nn.eq.0) goto 310
        do 300 i=1,nn
          if(feasb.and.i.gt.nineqn)     signgj=signeq(i-nineqn)
          if(.not.feasb.or.i.le.nineqn) signgj=1.d0
          valnom=g(indxcn(i))*signgj
          call gradcn(nparam,indxcn(i),xnew,gradg(1,indxcn(i)),constr)
 300    continue
        call resign(nparam,neqn,psf,grdpsf,penp,g(nnineq+1),
     *              gradg(1,nnineq+1),signeq,11,11)
 310    do 400 i=1,nob
          valnom=f(i)
          if(feasb) call gradob(nparam,i,xnew,gradf(1,i),obj)
          if(.not.feasb)
     *      call gradcn(nparam,indxob(i),xnew,gradf(1,i),constr)
 400    continue
 405    done=.true.
      goto 100
c
 410  if(d0is0) goto 910 
      if(nrst.lt.(5*nparam).or.steps.gt.0.1d0) goto 420
        nrst=0
        call diagnl(nparam,1.d0,hess)
        goto 810
 420  nrst=nrst+1
      do 500 i=1,nparam
        gamma(i)=gamma(i)-delta(i)
        delta(i)=xnew(i)-x(i)
 500  continue
      call matrvc(nparam,nparam,nparam,nparam,hess,delta,hd)
      dhd=scaprd(nparam,delta,hd)
      if(sqrt(scaprd(nparam,delta,delta)).le.epsmac) then
        nstop=0
        info=8
        goto 9000
      endif
      gammd=scaprd(nparam,delta,gamma)
      if(gammd.ge.0.2d0*dhd) theta=1.d0
      if(gammd.lt.0.2d0*dhd) theta=.8d0*dhd/(dhd-gammd)
      do 600 i=1,nparam
 600    eta(i)=hd(i)*(1.d0-theta)+theta*gamma(i)
      etad=theta*gammd+(1.d0-theta)*dhd
      do 800  i=1,nparam
        do 700 j=i,nparam
          hess(i,j)=hess(i,j)-hd(i)*hd(j)/dhd+eta(i)*eta(j)/etad
 700    hess(j,i)=hess(i,j)
 800  continue
 810  do 900 i=1,nparam
 900    x(i)=xnew(i)
      xisnew=.true.
 910  if(nstop.eq.0) goto 9000
      if(neqn.eq.0.or..not.feasb) goto 1400
        iout=io
        i=neq-neqn
        if(i.eq.0) goto 940
        call matrvc(nparam,i,nparam,i,gradg(1,nnineq+neqn+1),
     *              cllamd(nparam+ng+nnineq+neqn),gamma)
        do 930 i=1,nparam
 930      psb(i)=psb(i)+gamma(i)
 940    i=nnineq-nineqn
        if(i.eq.0) goto 990
        call matrvc(nparam,i,nparam,i,gradg(1,nineqn+1),
     *              cllamd(nparam+ng+nineqn),gamma)
        do 950 i=1,nparam
 950      psb(i)=psb(i)+gamma(i)
 990    call estlam(nparam,neqn,ifail,iout,bigbnd,phess,delta,eta,gamma,
     *              gradg(1,nnineq+1),psb,hd,xnew,psmu,w,lenw,iw,leniw)
        do 1000 i=1,neqn
          if(ifail.ne.0.or.d0is0) then
            penp(i)=2.d0*penp(i)
          else if(ifail.eq.0) then
            etad=psmu(i)+penp(i)
            if(etad.ge.1.d0) goto 1000
            penp(i)=dmax1(1.0d0-psmu(i),2.0d0*penp(i))
          endif
          if(penp(i).gt.bigbnd) then
            nstop=0
            info=9
            goto 9000
          endif
 1000   continue
        call resign(nparam,neqn,psf,grdpsf,penp,g(nnineq+1),
     *              gradg(1,nnineq+1),signeq,20,12)
        fMp=fM-psf
 1400   if(nfs.eq.0) goto 1430
        nstart=nstart+1
        np=indexs(nstart,nfs)
        span(np)=fmax
        do 1410 i=1,neqn
 1410     gm((np-1)*neqn+i)=g(nnineq+i)
        if(neqn.ne.0) call resign(nparam,neqn,psfnew,grdpsf,penp,
     *                            gm(1),gradg,signeq,20,10)
        fM=span(1)
        fMp=span(1)-psfnew
        mnm=min0(nstart,nfs)
        do 1420 i=2,mnm
          if(neqn.ne.0) call resign(nparam,neqn,psfnew,grdpsf,penp,
     *                           gm((i-1)*neqn+1),gradg,signeq,20,10)
          fM=dmax1(fM,span(i))
          fMp=dmax1(fMp,span(i)-psfnew)
 1420   continue
 1430 if(iprint.lt.3.or.ipyes.gt.0) goto 9000
        do 1700 i=1,nob
          if(.not.feasb) goto 1600
            if(nob.gt.1) call sbout2(io,nparam,i,'gradf(j,',')',
     *                               gradf(1,i))
            if(nob.eq.1) call sbout1(io,nparam,'gradf(j)            ',
     *                               dummy,gradf(1,i),2,2)
          goto 1700
 1600       call sbout2(io,nparam,indxob(i),'gradg(j,',')',
     *                  gradf(1,i))
 1700   continue
        if(ncnstr.eq.0) goto 1900
        do 1800 i=1,ncnstr
 1800     call sbout2(io,nparam,indxcn(i),'gradg(j,',')',
     *                gradg(1,indxcn(i)))
        if(neqn.eq.0) goto 1900
        call sbout1(io,nparam,'grdpsf(j)           ',dummy,grdpsf,2,2)
        call sbout1(io,neqn,'P                   ',dummy,penp,2,2)
c       call sbout1(io,neqn,'psmu                ',dummy,psmu,2,2)
 1900   call sbout1(io,nparam,'multipliers  for  x ',dummy,cllamd,2,2)
        if(ncnstr.ne.0) call sbout1(io,ncnstr,'             for  g ',
     *                              dummy,cllamd(nparam+ng),2,2)
        if(nobL.gt.1) call sbout1(io,nob,'             for  f ',
     *                            dummy,cllamd(nparam+ng+ncnstr),2,2)
        do 2000 i=1,nparam
 2000     call sbout2(io,nparam,i,'hess (j,',')',hess(1,i))
 9000 return 
      end
      subroutine grobfd(nparam,j,x,gradf,obj)
c
c     FFSQP : computation of gradients of objective
c             functions by forward finite differences
c
c     implicit real*8(a-h,o-z)
      integer nparam,j
      double  precision x(nparam),gradf(nparam)
      external obj
c
      integer io,iprint,ipspan,ipyes,info,ipd,idum,idum2
      double  precision epsmac,rteps,udelta,fj,objeps,objrep,gLgeps
      logical xisnew
      common  /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,idum,idum2,
     *        /fsqpp3/epsmac,rteps,udelta,fj
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
c     estimates the gradient of the objective function 
c     by forward finite differences
c
      integer i
      double  precision xi,delta,dmax1
c
      do 10 i=1,nparam
        xi=x(i)
        delta=dmax1(udelta,rteps*dmax1(1.d0,dabs(xi)))
        if (xi.lt.0.d0) delta=-delta
        if (ipd.eq.1.or.j.ne.1.or.iprint.lt.3.or.ipyes.gt.0) goto 9
          if(i.eq.1) write(io,1001) delta
          if(i.ne.1) write(io,1002) delta
  9     x(i)=xi+delta
        xisnew=.true.
        call obj(nparam,j,x,gradf(i))
        gradf(i)=(gradf(i)-fj)/delta
        x(i)=xi
        xisnew=.true.
 10   continue
      return
 1001 format(1x,t17,8hdelta(i),t45,e22.14)
 1002 format(1x,t45,e22.14)
      end
c
      subroutine grcnfd(nparam,j,x,gradg,constr)
c
c     FFSQP : computation of gradients of constraint
c             functions by forward finite differences
c
c     implicit real*8(a-h,o-z)
      integer nparam,j
      double  precision x(nparam),gradg(nparam)
      external constr
c
      integer io,iprint,ipspan,ipyes,info,ipd,idum,idum2
      double  precision epsmac,rteps,udelta,gj,objeps,objrep,gLgeps
      logical xisnew
      common  /fsqpp2/io,iprint,ipspan,ipyes,info,ipd,idum,idum2,
     *        /fsqpp3/epsmac,rteps,udelta,gj
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
c     estimate the gradient of the ith constraint 
c     by forward finite differences
c
      integer i
      double  precision xi,delta,dmax1
c
      do 10 i=1,nparam
        xi=x(i)
        delta=dmax1(udelta,rteps*dmax1(1.d0,dabs(xi)))
        if (xi.lt.0.d0) delta=-delta
        if (j.ne.1.or.iprint.lt.3) goto 9
        if (ipspan.ge.10.and.ipyes.gt.0) goto 9
          if(i.eq.1) write(io,1001) delta
          if(i.ne.1) write(io,1002) delta
        ipd=1
  9     x(i)=xi+delta
        xisnew=.true.
        call constr(nparam,j,x,gradg(i))
        gradg(i)=(gradg(i)-gj)/delta
        x(i)=xi
        xisnew=.true.
 10   continue
      return
 1001 format(1x,t17,8hdelta(i),t45,e22.14)
 1002 format(1x,t45,e22.14)
      end
      subroutine out(miter,nparam,nob,nineqn,nn,neqn,ncnstr,x,g,f,
     *               fmax,fM,psf,steps,sktnom,d0norm,feasb)
c
c     FFSQP : output for different value of iprint
c
c     implicit real*8(a-h,o-z)
      integer miter,nparam,nob,nineqn,nn,neqn,ncnstr
      double  precision fmax,fM,steps,sktnom,d0norm,psf
      double  precision x(nparam),g(1),f(1)
c     double  precision x(nparam),g(ncnstr),f(nob)
      logical feasb
c
      integer nnineq,M,ncallg,ncallf,mode,io,iprint,ipspan,ipyes,
     *        info,idum1,iter,nstop,initvl,lstype
      common  /fsqpp1/nnineq,M,ncallg,ncallf,mode,lstype,nstop,
     *        /fsqpp2/io,iprint,ipspan,ipyes,info,idum1,iter,initvl
c
      integer i
      double precision SNECV,dummy,adummy(1)
c
      if(nstop.eq.0) ipyes=0
      if (iter.ge.miter.and.nstop.ne.0) then
        info=3
        nstop=0
      endif
 10   if(iprint.eq.0.or.ipyes.gt.0) then
        iter=iter+1
        goto 9000
      endif
      if(info.gt.0.and.info.lt.3) goto 120
      if(iprint.ne.1.or.nstop.eq.0) goto 20
        iter=iter+1
        if(initvl.eq.0) goto 9000
        if(feasb.and.nob.gt.0)
     *    call sbout1(io,nob,'objectives          ',dummy,f,2,1)
        if (mode.eq.1.and.iter.gt.1.and.feasb)
     *  call sbout1(io,0,'objective max4      ',fM,adummy,1,1)
        if(nob.gt.1) call sbout1(io,0,'objmax              ',
     *                           fmax,adummy,1,1)
        if(ncnstr.eq.0) write(io,9909)
        call sbout1(io,ncnstr,'constraints         ',dummy,g,2,1)
        if(ncnstr.ne.0) write(io,9909)
        goto 9000
 20   if(iprint.eq.1.and.nstop.eq.0) write(io,9900) iter
      if(iprint.ge.2.and.nstop.eq.0.and.ipspan.ge.10) 
     *  write(io,9900) iter
      iter=iter+1
      if(initvl.eq.0) 
     *  call sbout1(io,nparam,'x                   ',dummy,x,2,1)
      if(nob.gt.0)
     *  call sbout1(io,nob,'objectives          ',dummy,f,2,1)
      if (mode.eq.1.and.iter.gt.1)
     *  call sbout1(io,0,'objective max4      ',fM,adummy,1,1)
      if(nob.gt.1) call sbout1(io,0,'objmax              ',
     *                          fmax,adummy,1,1)
      if(ncnstr.eq.0) go to 110
      call sbout1(io,ncnstr,'constraints         ',dummy,g,2,1)
      if(.not.feasb) goto 110
        SNECV=0.d0
        do 100 i=nnineq+1,nnineq+nn-nineqn
          SNECV=SNECV+dabs(g(i))
 100    continue
        if(initvl.eq.0.and.(nn-nineqn).ne.0)
     *    call sbout1(io,0,'SNECV               ',SNECV,adummy,1,1)
 110  continue
      if(iter.le.1) write(io,9909)
      if(iter.le.1.and.ipspan.lt.10) write(io,9900) iter
      if(iter.le.1) goto 9000
      if(iprint.ge.2.and.initvl.eq.0)
     *  call sbout1(io,0,'step                ',steps,adummy,1,1)
      if(initvl.eq.0.and.(nstop.eq.0.or.info.ne.0.or.iprint.eq.2)) then
        call sbout1(io,0,'d0norm              ',d0norm,adummy,1,1)
        call sbout1(io,0,'ktnorm              ',sktnom,adummy,1,1)
      endif
      if(initvl.eq.0.and.feasb) write(io,9902) ncallf
      if(initvl.eq.0.and.(nn.ne.0.or..not.feasb)) write(io,9903) ncallg
      if(nstop.ne.0) write(io,9909)
      if(nstop.ne.0.and.iter.le.miter.and.ipspan.lt.10) 
     *  write(io,9900) iter
 120  if(nstop.ne.0.or.iprint.eq.0) goto 9000
      write(io,9909)
      write(io,9901) info
      if(info.eq.0) write(io,9904)
      if(info.eq.0.and.sktnom.gt.0.1d0) write(io,9910)
      if(info.eq.3) write(io,9905)
      if(info.eq.4) write(io,9906)
      if(info.eq.5) write(io,9907)
      if(info.eq.6) write(io,9908)
      if(info.eq.8) write(io,9911)
      if(info.eq.9) write(io,9912)
      write(io,9909)
 9000 initvl=0
      if(ipspan.ge.10) ipyes=mod(iter,ipspan)
c      if(iter.le.miter) return
c        nstop=0
c        info=3
c        write(io,9905)
      return
 9900 format(1x,9hiteration,t22,i22)
 9901 format(1x,6hinform,t22,i22)
 9902 format(1x,6hncallf,t22,i22)
 9903 format(1x,6hncallg,t22,i22)
 9904 format(1x,'Normal termination: You have obtained a solution !!')
 9905 format(1x,'Warning : Maximum iterations have been reached ',
     *          'before obtaining a solution !!'/)
 9906 format(1x,'Error : Step size has been smaller than ',
     *          'the computed machine precision !!'/)
 9907 format(1x,'Error : Failure of the QP solver ',
     *          'in constructing d0 !!',
     *      /1x,'        A more robust QP solver may succeed.'/)
 9908 format(1x,'Error : Failure of the QP solver ',
     *          'in constructing d1 !!',
     *      /1x,'        A more robust QP solver may succeed.'/)
 9909 format(1x,/)
 9910 format(1x,'Warning: Norm of Kuhn-Tucker vector is large !!'/)
 9911 format(1x,'Error : The new iterate is numerically equivalent to ',
     *      /1x,'the current iterate, though the stopping criterion',
     *      /1x,'is not satisfied.'/)
 9912 format(1x,'Error : Could not satisfy nonlinear equality',
     *      /1x,'constraints - penalty parameter too large.'/)
      end
c=== subroutines used in FFSQP ====================================c
c                                                                  c
c  diagnl  error   estlam  fool    indexs  lfuscp  matrcp  matrvc  c
c  nullvc  resign  sbout1  sbout2  scaprd  shift   slope   small   c 
c                                                                  c
c==================================================================c
c
      subroutine diagnl(nrowa,diag,a)
c     implicit real*8(a-h,o-z)
      integer nrowa,i,j
      double  precision a(nrowa,1),diag
c     double  precision a(nrowa,nrowa),diag
c
c     set a=diag*I, the diagonal matrix
c
      do 200 i=1,nrowa
        do 100 j=i,nrowa
          a(i,j)=0.d0
 100      a(j,i)=0.d0
 200    a(i,i)=diag
      return
      end
c
      subroutine error(string,inform,io)
c     implicit real*8 (a-h,o-z)
      integer inform,io
      character*40 string
c
      write(io,9900) string
 9900 format(1x,a40)
      inform=7
      return
      end
c
      subroutine estlam(nparam,neqn,ifail,iout,bigbnd,hess,cvec,a,b,
     *                  gradh,psb,bl,bu,x,w,lenw,iw,leniw)
      integer nparam,neqn,ifail,iout,lenw,leniw,iw(leniw)
      double precision bigbnd,hess(neqn,1),cvec(1),a(1),b(1),
     *                 gradh(nparam,1),psb(1),bl(1),bu(1),
     *                 x(1),w(lenw)
c     double precision bigbnd,hess(neqn,neqn),cvec(neqn),a(1),b(1),
c    *                 gradh(nparam,neqn),psb(nparam),bl(1),bu(1),
c    *                 x(neqn),w(lenw)
c
c     compute an estimate of multipliers for updating penalty parameter
c
      integer i,j
      double precision scaprd
c
      do 200 i=1,neqn
        bl(i)=-bigbnd
        bu(i)=bigbnd
        cvec(i)=scaprd(nparam,gradh(1,i),psb)
        x(i)=0.d0
        do 100 j=i,neqn 
          hess(i,j)=scaprd(nparam,gradh(1,i),gradh(1,j))
 100      hess(j,i)=hess(i,j)
 200  continue
      iw(1)=1
      call ql0001(0,0,1,neqn,neqn,2*neqn,hess,cvec,a,b,bl,bu,x,w,
     c            iout,ifail,0,w(2),lenw-1,iw,leniw)
c
      return
      end
c
      subroutine fool(x,y,z)
      double precision x,y,z
c
      z=x*y+y
      return
      end
c
      double precision function lfuscp(val,thrshd)
c     implicit real*8(a-h,o-z)
      double precision val,thrshd
c
      if(dabs(val).le.thrshd) lfuscp=0
      if(dabs(val).gt.thrshd) lfuscp=1
      return
      end
c
      integer function indexs(i,nfs)
c     implicit real*8(a-h,o-z)
      integer i,nfs,mm
c
c     find the residue of i with respect to nfs
c
      mm=i
      if(mm.le.nfs) goto 120
 110  mm=mm-nfs
      if(mm.gt.nfs) goto 110
 120  indexs=mm
      return
      end
c
      subroutine matrcp(ndima,a,ndimb,b)
c     implicit real*8(a-h,o-z)
      integer ndima,ndimb,i,j
      double  precision a(ndima,1),b(ndimb,1)
c     double  precision a(ndima,ndima),b(ndimb,ndimb)
c
      do 100 i=1,ndima
        do 100 j=1,ndima
 100      b(i,j)=a(i,j)
      if(ndimb.le.ndima) goto 9000
        do 200 i=1,ndimb
          b(ndimb,i)=0.d0
 200      b(i,ndimb)=0.d0
 9000 return
      end
c 
      subroutine matrvc(l,n,la,na,a,x,y)
c     implicit real*8(a-h,o-z)
      integer l,n,la,na,i,j
      double  precision a(l,n),x(n),y(l),yi
c     double  precision a(l,1),x(1),y(1),yi
c
c     computation of y=ax
c
      do 200 i=1,la
        yi=0.d0
        do 100 j=1,na
 100      yi=yi+a(i,j)*x(j)
 200      y(i)=yi
      return
      end       
c
      subroutine nullvc(nparam,x)
c     implicit real*8(a-h,o-z)
      integer nparam,i
      double  precision x(nparam)
c
c     set x=0
c
      do 100 i=1,nparam
 100    x(i)=0.d0
      return
      end
c
      subroutine resign(n,neqn,psf,grdpsf,penp,g,gradg,signeq,job1,job2)
      integer i,j,job1,job2,n,neqn
      double precision psf,grdpsf(1),penp(1),g(1),gradg(n,1),
     *                 signeq(1)
c     double precision psf,grdpsf(n),penp(neqn),g(neqn),gradg(n,neqn),
c    *                 signeq(neqn)
c
c     job1=10: g*signeq, job1=11: gradg*signeq, job1=12: job1=10&11
c     job1=20: do not change sign
c     job2=10: psf,      job2=11: grdpsf,       job2=12: job2=10&11
c     job2=20: do not compute psf or grdpsf
c
      if(job2.eq.10.or.job2.eq.12) psf=0.d0
      do 100 i=1,neqn
        if(job1.eq.10.or.job1.eq.12) g(i)=signeq(i)*g(i)
        if(job2.eq.10.or.job2.eq.12) psf=psf+g(i)*penp(i)
        if(job1.eq.10.or.job1.eq.20) goto 100
          do 50 j=1,n
            gradg(j,i)=gradg(j,i)*signeq(i)
  50      continue
 100  continue
      if(job2.eq.10.or.job2.eq.20) goto 9000
      call nullvc(n,grdpsf)
      do 120 i=1,n
        do 110 j=1,neqn
 110      grdpsf(i)=grdpsf(i)+gradg(i,j)*penp(j)
 120  continue
c
 9000 return
      end
c
      subroutine sbout1(io,n,s1,z,z1,job,level)
c     implicit real*8(a-h,o-z)
      integer io,n,job,level,j
      double precision z,z1(1)
      character*20 s1
c
      if (job.eq.2) goto 10
      if (level.eq.1)write(io,9900) s1,z
      if (level.eq.2)write(io,9901) s1,z
      return
 10   if(n.eq.0) goto 101
      if (level.eq.1)write(io,9900) s1,z1(1)
      if (level.eq.2)write(io,9901) s1,z1(1)
      do 100 j=2,n
        if (level.eq.1) write(io,9902) z1(j)
        if (level.eq.2) write(io,9903) z1(j)
 100  continue
 101  return
 9900 format(1x,a20,e22.14)
 9901 format(1x,t17,a20,t45,e22.14)
 9902 format(1x,t22,e22.14)
 9903 format(1x,t45,e22.14)
      end
c
      subroutine sbout2(io,n,i,s1,s2,z)
c     implicit real*8(a-h,o-z)
      integer io,n,i,j
      double precision z(n)
      character*8 s1
      character*1 s2
c
      write(io,9900) s1,i,s2,z(1)
      do 100 j=2,n
 100    write(io,9901) z(j)
      return
 9900 format(1x,t17,a8,i5,a1,t45,e22.14)
 9901 format(1x,t45,e22.14)
      end
c
      double  precision function scaprd(n,x,y)
c     implicit real*8(a-h,o-z)
      integer n,i
      double  precision x(1),y(1),z
c     double  precision x(n),y(n),z
c
c     compute z=x'y
c
      z=0.d0
      do 100 i=1,n
        z=x(i)*y(i)+z 
 100  continue
      scaprd=z
      return
      end
c
      subroutine shift(n,ii,iact)
      integer n,ii,iact(1),j,k
c
      if(ii.eq.iact(1)) return
      do 200 j=1,n
        if(ii.ne.iact(j)) goto 200
        do 100 k=j,2,-1
 100      iact(k)=iact(k-1)
        goto 210
 200  continue
 210  if(n.ne.0) iact(1)=ii
      return
      end
c
      double precision function slope(nob,nobL,neqn,nparam,feasb,
     *                         f,gradf,grdpsf,x,y,fmax,theta,job)
c     implicit real*8(a-h,o-z)
      integer nob,nobL,neqn,nparam,job,i
      double  precision fmax,theta,slope1,dmax1,dmin1,rhs,rhog,
     *        grdftx,grdfty,diff,scaprd,grpstx,grpsty
      double  precision f(nob),gradf(nparam,nob),grdpsf(nparam),
     *        x(nparam),y(nparam)
c     double  precision f(1),gradf(nparam,1),grdpsf(nparam),
c    *        x(nparam),y(nparam)
      logical feasb
c
      double  precision bigbnd,dummy
      common  /fsqpq1/bigbnd,dummy
c
c     job=0 : compute the generalized gradient of the minimax
c     job=1 : compute rhog in mode = 1
c
      slope=-bigbnd
      if(feasb.and.nob.eq.0) slope=0.d0
      if(neqn.eq.0.or..not.feasb) then
        grpstx=0.d0
        grpsty=0.d0
      else
        grpstx=scaprd(nparam,grdpsf,x)
        grpsty=scaprd(nparam,grdpsf,y)
      endif
      do 100 i=1,nob
        slope1=f(i)+scaprd(nparam,gradf(1,i),x)
        slope=dmax1(slope,slope1)
        if(nobL.ne.nob) slope=dmax1(slope,-slope1)
 100  continue
      slope=slope-fmax-grpstx
      if (job.eq.0) goto 9000
      rhs=theta*slope+fmax
      rhog=1.d0
      do 200 i=1,nob
        grdftx=scaprd(nparam,gradf(1,i),x)-grpstx
        grdfty=scaprd(nparam,gradf(1,i),y)-grpsty
        diff=grdfty-grdftx
        if (diff.le.0.d0) goto 200
        rhog=dmin1(rhog,(rhs-f(i)-grdftx)/diff)
        if(nobL.ne.nob) rhog=dmin1(rhog,-(rhs+f(i)+grdftx)/diff)
 200  continue
      slope=rhog
 9000 return
      end
c
      double precision function small()
c     implicit real*8(a-h,o-z)
      double precision one, two, z
c
      one=1.d0
      two=2.d0
      small=one
10    small=small/two
      call fool(small,one,z)
      if(z.gt.one) goto 10
      small=small*two*two
c
c The simpler sequence commented out below fails on some machines that use
c extra-length registers for internal computation.  This was pointed out
c to us by Roque Donizete de Oliveira (Michigan) who suggested to sequence
c used now.
c
c     small=1.d0
c100  if ((small+1.d0).eq.1.d0) goto 110
c     small=small/2.d0
c     goto 100
c110  small=small*4.d0
      return
      end
c
      block data
      double  precision objeps,objrep,gLgeps
      logical xisnew
      common  /fsqpus/objeps,objrep,gLgeps,xisnew
c
      data objeps,objrep,gLgeps/-1.d0,-1.d0,-1.d0/
      end

      SUBROUTINE QL0001(M,ME,MMAX,N,NMAX,MNN,C,D,A,B,XL,XU,
     1           X,U,IOUT,IFAIL,IPRINT,WAR,LWAR,IWAR,LIWAR)
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c                     !!!! NOTICE !!!!
c
c 1. The routines contained in this file are due to Prof. K.Schittkowski
c    of the University of Bayreuth, Germany (modification of routines
c    due to Prof. MJD Powell at the University of Cambridge).  They can
c    be freely distributed.
c
c 2. A minor modification was performed at the University of Maryland. 
c    It is marked in the code by "c umd".
c
c                                      A.L. Tits and J.L. Zhou
c                                      University of Maryland
C
C***********************************************************************
C
C
C             SOLUTION OF QUADRATIC PROGRAMMING PROBLEMS
C
C
C
C   QL0001 SOLVES THE QUADRATIC PROGRAMMING PROBLEM
C
C   MINIMIZE        .5*X'*C*X + D'*X
C   SUBJECT TO      A(J)*X  +  B(J)   =  0  ,  J=1,...,ME
C                   A(J)*X  +  B(J)  >=  0  ,  J=ME+1,...,M
C                   XL  <=  X  <=  XU
C   
C   HERE C MUST BE AN N BY N SYMMETRIC AND POSITIVE MATRIX, D AN N-DIMENSIONAL
C   VECTOR, A AN M BY N MATRIX AND B AN M-DIMENSIONAL VECTOR. THE ABOVE 
C   SITUATION IS INDICATED BY IWAR(1)=1. ALTERNATIVELY, I.E. IF IWAR(1)=0,
C   THE OBJECTIVE FUNCTION MATRIX CAN ALSO BE PROVIDED IN FACTORIZED FORM.
C   IN THIS CASE, C IS AN UPPER TRIANGULAR MATRIX.
C
C   THE SUBROUTINE REORGANIZES SOME DATA SO THAT THE PROBLEM CAN BE SOLVED
C   BY A MODIFICATION OF AN ALGORITHM PROPOSED BY POWELL (1983).
C
C
C   USAGE:
C
C      QL0001(M,ME,MMAX,N,NMAX,MNN,C,D,A,B,XL,XU,X,U,IOUT,IFAIL,IPRINT,
C             WAR,LWAR,IWAR,LIWAR)
C
C
C   DEFINITION OF THE PARAMETERS:
C
C   M :        TOTAL NUMBER OF CONSTRAINTS.
C   ME :       NUMBER OF EQUALITY CONSTRAINTS.
C   MMAX :     ROW DIMENSION OF A. MMAX MUST BE AT LEAST ONE AND GREATER
C              THAN M.
C   N :        NUMBER OF VARIABLES.
C   NMAX :     ROW DIMENSION OF C. NMAX MUST BE GREATER OR EQUAL TO N.
C   MNN :      MUST BE EQUAL TO M + N + N.
C   C(NMAX,NMAX): OBJECTIVE FUNCTION MATRIX WHICH SHOULD BE SYMMETRIC AND
C              POSITIVE DEFINITE. IF IWAR(1) = 0, C IS SUPPOSED TO BE THE
C              CHOLESKEY-FACTOR OF ANOTHER MATRIX, I.E. C IS UPPER
C              TRIANGULAR.
C   D(NMAX) :  CONTAINS THE CONSTANT VECTOR OF THE OBJECTIVE FUNCTION.
C   A(MMAX,NMAX): CONTAINS THE DATA MATRIX OF THE LINEAR CONSTRAINTS.
C   B(MMAX) :  CONTAINS THE CONSTANT DATA OF THE LINEAR CONSTRAINTS.
C   XL(N),XU(N): CONTAIN THE LOWER AND UPPER BOUNDS FOR THE VARIABLES.
C   X(N) :     ON RETURN, X CONTAINS THE OPTIMAL SOLUTION VECTOR.
C   U(MNN) :   ON RETURN, U CONTAINS THE LAGRANGE MULTIPLIERS. THE FIRST
C              M POSITIONS ARE RESERVED FOR THE MULTIPLIERS OF THE M
C              LINEAR CONSTRAINTS AND THE SUBSEQUENT ONES FOR THE 
C              MULTIPLIERS OF THE LOWER AND UPPER BOUNDS. ON SUCCESSFUL
C              TERMINATION, ALL VALUES OF U WITH RESPECT TO INEQUALITIES 
C              AND BOUNDS SHOULD BE GREATER OR EQUAL TO ZERO.
C   IOUT :     INTEGER INDICATING THE DESIRED OUTPUT UNIT NUMBER, I.E.
C              ALL WRITE-STATEMENTS START WITH 'WRITE(IOUT,... '.
C   IFAIL :    SHOWS THE TERMINATION REASON.
C      IFAIL = 0 :   SUCCESSFUL RETURN.
C      IFAIL = 1 :   TOO MANY ITERATIONS (MORE THAN 40*(N+M)).
C      IFAIL = 2 :   ACCURACY INSUFFICIENT TO SATISFY CONVERGENCE
C                    CRITERION.
C      IFAIL = 5 :   LENGTH OF A WORKING ARRAY IS TOO SHORT.
C      IFAIL > 10 :  THE CONSTRAINTS ARE INCONSISTENT.
C   IPRINT :   OUTPUT CONTROL.
C      IPRINT = 0 :  NO OUTPUT OF QL0001.
C      IPRINT > 0 :  BRIEF OUTPUT IN ERROR CASES.
C   WAR(LWAR) : REAL WORKING ARRAY. THE LENGTH LWAR SHOULD BE GRATER THAN
C               3*NMAX*NMAX/2 + 10*NMAX + 2*MMAX.
C   IWAR(LIWAR): INTEGER WORKING ARRAY. THE LENGTH LIWAR SHOULD BE AT
C              LEAST N.
C              IF IWAR(1)=0 INITIALLY, THEN THE CHOLESKY DECOMPOSITION
C              WHICH IS REQUIRED BY THE DUAL ALGORITHM TO GET THE FIRST
C              UNCONSTRAINED MINIMUM OF THE OBJECTIVE FUNCTION, IS
C              PERFORMED INTERNALLY. OTHERWISE, I.E. IF IWAR(1)=1, THEN
C              IT IS ASSUMED THAT THE USER PROVIDES THE INITIAL FAC-
C              TORIZATION BY HIMSELF AND STORES IT IN THE UPPER TRIAN-
C              GULAR PART OF THE ARRAY C.
C
C   A NAMED COMMON-BLOCK  /CMACHE/EPS   MUST BE PROVIDED BY THE USER,
C   WHERE EPS DEFINES A GUESS FOR THE UNDERLYING MACHINE PRECISION.
C
C*********************************************************************
C
C
      INTEGER NMAX,MMAX,N,MNN,LWAR,LIWAR
      DIMENSION C(NMAX,NMAX),D(NMAX),A(MMAX,NMAX),B(MMAX),
     1      XL(N),XU(N),X(N),U(MNN),WAR(LWAR),IWAR(LIWAR)
      DOUBLE PRECISION C,D,A,B,X,XL,XU,U,WAR,DIAG,ZERO,
     1      EPS,QPEPS,TEN
      INTEGER M,ME,IOUT,IFAIL,IPRINT,IWAR,INW1,INW2,IN,J,LW,MN,I,         
     1      IDIAG,INFO,NACT,MAXIT
      LOGICAL LQL
C
C     INTRINSIC FUNCTIONS:  DSQRT
C
      COMMON /CMACHE/EPS
C
C     CONSTANT DATA
C
c#################################################################
c

      if(dabs(c(nmax,nmax)).lt.eps) c(nmax,nmax)=eps
c
c umd
c  This prevents a subsequent more major modification of the Hessian
c  matrix in the important case when a minmax problem (yielding a 
c  singular Hessian matrix) is being solved.
c                                 ----UMCP, April 1991, Jian L. Zhou
c#################################################################
c
      LQL=.FALSE.
      IF (IWAR(1).EQ.1) LQL=.TRUE.
      ZERO=0.0D+0
      TEN=1.D+1
      MAXIT=40*(M+N)
      QPEPS=EPS
      INW1=1
      INW2=INW1+MMAX
C
C     PREPARE PROBLEM DATA FOR EXECUTION
C
      IF (M.LE.0) GOTO 20
      IN=INW1
      DO 10 J=1,M
      WAR(IN)=-B(J)
   10 IN=IN+1
   20 LW=3*NMAX*NMAX/2 + 10*NMAX + M
      IF ((INW2+LW).GT.LWAR) GOTO 80
      IF (LIWAR.LT.N) GOTO 81
      IF (MNN.LT.M+N+N) GOTO 82
      MN=M+N
C
C     CALL OF QL0002
C
      CALL QL0002(N,M,ME,MMAX,MN,MNN,NMAX,LQL,A,WAR(INW1),
     1    D,C,XL,XU,X,NACT,IWAR,MAXIT,QPEPS,INFO,DIAG,
     2    WAR(INW2),LW)
C
C     TEST OF MATRIX CORRECTIONS
C
      IFAIL=0
      IF (INFO.EQ.1) GOTO 40
      IF (INFO.EQ.2) GOTO 90
      IDIAG=0
      IF ((DIAG.GT.ZERO).AND.(DIAG.LT.1000.0)) IDIAG=DIAG
      IF ((IPRINT.GT.0).AND.(IDIAG.GT.0))
     1   WRITE(IOUT,1000) IDIAG
      IF (INFO .LT. 0) GOTO  70
C
C     REORDER MULTIPLIER
C
      DO  50 J=1,MNN
   50 U(J)=ZERO
      IN=INW2-1
      IF (NACT.EQ.0) GOTO 30
      DO  60 I=1,NACT
      J=IWAR(I)
      U(J)=WAR(IN+I)
   60 CONTINUE
   30 CONTINUE
      RETURN
C
C     ERROR MESSAGES
C
   70 IFAIL=-INFO+10
      IF ((IPRINT.GT.0).AND.(NACT.GT.0))
     1     WRITE(IOUT,1100) -INFO,(IWAR(I),I=1,NACT)
      RETURN
   80 IFAIL=5
      IF (IPRINT .GT. 0) WRITE(IOUT,1200)
      RETURN
   81 IFAIL=5
      IF (IPRINT .GT. 0) WRITE(IOUT,1210)
      RETURN
   82 IFAIL=5
      IF (IPRINT .GT. 0) WRITE(IOUT,1220)
      RETURN
   40 IFAIL=1
      IF (IPRINT.GT.0) WRITE(IOUT,1300) MAXIT
      RETURN
   90 IFAIL=2
      IF (IPRINT.GT.0) WRITE(IOUT,1400) 
      RETURN
C
C     FORMAT-INSTRUCTIONS
C
 1000 FORMAT(/8X,28H***QL: MATRIX G WAS ENLARGED,I3,
     1        20H-TIMES BY UNITMATRIX)
 1100 FORMAT(/8X,18H***QL: CONSTRAINT ,I5,
     1        19H NOT CONSISTENT TO ,/,(10X,10I5))
 1200 FORMAT(/8X,21H***QL: LWAR TOO SMALL)
 1210 FORMAT(/8X,22H***QL: LIWAR TOO SMALL)
 1220 FORMAT(/8X,20H***QL: MNN TOO SMALL)
 1300 FORMAT(/8X,37H***QL: TOO MANY ITERATIONS (MORE THAN,I6,1H))
 1400 FORMAT(/8X,50H***QL: ACCURACY INSUFFICIENT TO ATTAIN CONVERGENCE) 
      END
C
      SUBROUTINE QL0002(N,M,MEQ,MMAX,MN,MNN,NMAX,LQL,A,B,GRAD,G,
     1      XL,XU,X,NACT,IACT,MAXIT,VSMALL,INFO,DIAG,W,LW)
C
C**************************************************************************
C
C
C   THIS SUBROUTINE SOLVES THE QUADRATIC PROGRAMMING PROBLEM 
C
C       MINIMIZE      GRAD'*X  +  0.5 * X*G*X
C       SUBJECT TO    A(K)*X  =  B(K)   K=1,2,...,MEQ,
C                     A(K)*X >=  B(K)   K=MEQ+1,...,M,
C                     XL  <=  X  <=  XU
C
C   THE QUADRATIC PROGRAMMING METHOD PROCEEDS FROM AN INITIAL CHOLESKY-
C   DECOMPOSITION OF THE OBJECTIVE FUNCTION MATRIX, TO CALCULATE THE
C   UNIQUELY DETERMINED MINIMIZER OF THE UNCONSTRAINED PROBLEM. 
C   SUCCESSIVELY ALL VIOLATED CONSTRAINTS ARE ADDED TO A WORKING SET 
C   AND A MINIMIZER OF THE OBJECTIVE FUNCTION SUBJECT TO ALL CONSTRAINTS 
C   IN THIS WORKING SET IS COMPUTED. IT IS POSSIBLE THAT CONSTRAINTS
C   HAVE TO LEAVE THE WORKING SET.
C
C
C   DESCRIPTION OF PARAMETERS:    
C
C     N        : IS THE NUMBER OF VARIABLES.
C     M        : TOTAL NUMBER OF CONSTRAINTS.
C     MEQ      : NUMBER OF EQUALITY CONTRAINTS.
C     MMAX     : ROW DIMENSION OF A, DIMENSION OF B. MMAX MUST BE AT
C                LEAST ONE AND GREATER OR EQUAL TO M.
C     MN       : MUST BE EQUAL M + N.
C     MNN      : MUST BE EQUAL M + N + N.
C     NMAX     : ROW DIEMSION OF G. MUST BE AT LEAST N.
C     LQL      : DETERMINES INITIAL DECOMPOSITION.
C        LQL = .FALSE.  : THE UPPER TRIANGULAR PART OF THE MATRIX G
C                         CONTAINS INITIALLY THE CHOLESKY-FACTOR OF A SUITABLE 
C                         DECOMPOSITION.
C        LQL = .TRUE.   : THE INITIAL CHOLESKY-FACTORISATION OF G IS TO BE
C                         PERFORMED BY THE ALGORITHM.
C     A(MMAX,NMAX) : A IS A MATRIX WHOSE COLUMNS ARE THE CONSTRAINTS NORMALS.
C     B(MMAX)  : CONTAINS THE RIGHT HAND SIDES OF THE CONSTRAINTS.
C     GRAD(N)  : CONTAINS THE OBJECTIVE FUNCTION VECTOR GRAD.
C     G(NMAX,N): CONTAINS THE SYMMETRIC OBJECTIVE FUNCTION MATRIX.
C     XL(N), XU(N): CONTAIN THE LOWER AND UPPER BOUNDS FOR X.
C     X(N)     : VECTOR OF VARIABLES.
C     NACT     : FINAL NUMBER OF ACTIVE CONSTRAINTS.
C     IACT(K) (K=1,2,...,NACT): INDICES OF THE FINAL ACTIVE CONSTRAINTS.
C     INFO     : REASON FOR THE RETURN FROM THE SUBROUTINE.
C         INFO = 0 : CALCULATION WAS TERMINATED SUCCESSFULLY.
C         INFO = 1 : MAXIMUM NUMBER OF ITERATIONS ATTAINED.
C         INFO = 2 : ACCURACY IS INSUFFICIENT TO MAINTAIN INCREASING
C                    FUNCTION VALUES.
C         INFO < 0 : THE CONSTRAINT WITH INDEX ABS(INFO) AND THE CON-
C                    STRAINTS WHOSE INDICES ARE IACT(K), K=1,2,...,NACT,
C                    ARE INCONSISTENT.
C     MAXIT    : MAXIMUM NUMBER OF ITERATIONS.
C     VSMALL   : REQUIRED ACCURACY TO BE ACHIEVED (E.G. IN THE ORDER OF THE 
C                MACHINE PRECISION FOR SMALL AND WELL-CONDITIONED PROBLEMS).
C     DIAG     : ON RETURN DIAG IS EQUAL TO THE MULTIPLE OF THE UNIT MATRIX
C                THAT WAS ADDED TO G TO ACHIEVE POSITIVE DEFINITENESS.
C     W(LW)    : THE ELEMENTS OF W(.) ARE USED FOR WORKING SPACE. THE LENGTH
C                OF W MUST NOT BE LESS THAN (1.5*NMAX*NMAX + 10*NMAX + M).
C                WHEN INFO = 0 ON RETURN, THE LAGRANGE MULTIPLIERS OF THE
C                FINAL ACTIVE CONSTRAINTS ARE HELD IN W(K), K=1,2,...,NACT.
C   THE VALUES OF N, M, MEQ, MMAX, MN, MNN AND NMAX AND THE ELEMENTS OF
C   A, B, GRAD AND G ARE NOT ALTERED.
C
C   THE FOLLOWING INTEGERS ARE USED TO PARTITION W:
C     THE FIRST N ELEMENTS OF W HOLD LAGRANGE MULTIPLIER ESTIMATES.
C     W(IWZ+I+(N-1)*J) HOLDS THE MATRIX ELEMENT Z(I,J).
C     W(IWR+I+0.5*J*(J-1)) HOLDS THE UPPER TRIANGULAR MATRIX
C       ELEMENT R(I,J). THE SUBSEQUENT N COMPONENTS OF W MAY BE
C       TREATED AS AN EXTRA COLUMN OF R(.,.).
C     W(IWW-N+I) (I=1,2,...,N) ARE USED FOR TEMPORARY STORAGE.
C     W(IWW+I) (I=1,2,...,N) ARE USED FOR TEMPORARY STORAGE.
C     W(IWD+I) (I=1,2,...,N) HOLDS G(I,I) DURING THE CALCULATION.
C     W(IWX+I) (I=1,2,...,N) HOLDS VARIABLES THAT WILL BE USED TO
C       TEST THAT THE ITERATIONS INCREASE THE OBJECTIVE FUNCTION.
C     W(IWA+K) (K=1,2,...,M) USUALLY HOLDS THE RECIPROCAL OF THE
C       LENGTH OF THE K-TH CONSTRAINT, BUT ITS SIGN INDICATES
C       WHETHER THE CONSTRAINT IS ACTIVE.
C
C*************************************************************************
C
      INTEGER MMAX,NMAX,N,LW
      DIMENSION A(MMAX,NMAX),B(MMAX),GRAD(N),G(NMAX,N),X(N),IACT(N),
     1      W(LW),XL(N),XU(N)
      INTEGER M,MEQ,MN,MNN,NACT,IACT,INFO,MAXIT
      DOUBLE PRECISION CVMAX,DIAG,DIAGR,FDIFF,FDIFFA,GA,GB,PARINC,PARNEW       
     1      ,RATIO,RES,STEP,SUM,SUMX,SUMY,SUMA,SUMB,SUMC,TEMP,TEMPA,
     2       VSMALL,XMAG,XMAGR,ZERO,ONE,TWO,ONHA,VFACT
      DOUBLE PRECISION A,B,G,GRAD,W,X,XL,XU
C
C   INTRINSIC FUNCTIONS:   DMAX1,DSQRT,DABS,DMIN1
C
      INTEGER IWZ,IWR,IWW,IWD,IWA,IFINC,KFINC,K,I,IA,ID,II,IR,IRA,
     1     IRB,J,NM,IZ,IZA,ITERC,ITREF,JFINC,IFLAG,IWS,IS,K1,IW,KK,IP,
     2     IPP,IL,IU,JU,KFLAG,LFLAG,JFLAG,KDROP,NU,MFLAG,KNEXT,IX,IWX,
     3     IWY,IY,JL
      LOGICAL LQL,LOWER
C
C   INITIAL ADDRESSES
C
      IWZ=NMAX
      IWR=IWZ+NMAX*NMAX
      IWW=IWR+(NMAX*(NMAX+3))/2
      IWD=IWW+NMAX
      IWX=IWD+NMAX
      IWA=IWX+NMAX
C
C     SET SOME CONSTANTS.
C
      ZERO=0.D+0
      ONE=1.D+0
      TWO=2.D+0
      ONHA=1.5D+0
      VFACT=1.D+0
C
C     SET SOME PARAMETERS.
C     NUMBER LESS THAN VSMALL ARE ASSUMED TO BE NEGLIGIBLE.
C     THE MULTIPLE OF I THAT IS ADDED TO G IS AT MOST DIAGR TIMES
C       THE LEAST MULTIPLE OF I THAT GIVES POSITIVE DEFINITENESS.
C     X IS RE-INITIALISED IF ITS MAGNITUDE IS REDUCED BY THE
C       FACTOR XMAGR.
C     A CHECK IS MADE FOR AN INCREASE IN F EVERY IFINC ITERATIONS,
C       AFTER KFINC ITERATIONS ARE COMPLETED.
C
      DIAGR=TWO
      XMAGR=1.0D-2
      IFINC=3
      KFINC=MAX0(10,N)
C
C     FIND THE RECIPROCALS OF THE LENGTHS OF THE CONSTRAINT NORMALS.
C     RETURN IF A CONSTRAINT IS INFEASIBLE DUE TO A ZERO NORMAL.
C
      NACT=0
      IF (M .LE. 0) GOTO 45
      DO 40 K=1,M
      SUM=ZERO
      DO 10 I=1,N
   10 SUM=SUM+A(K,I)**2
      IF (SUM .GT. ZERO) GOTO 20
      IF (B(K) .EQ. ZERO) GOTO 30
      INFO=-K
      IF (K .LE. MEQ) GOTO 730
      IF (B(K)) 30,30,730
   20 SUM=ONE/DSQRT(SUM)
   30 IA=IWA+K
   40 W(IA)=SUM
   45 DO 50 K=1,N
      IA=IWA+M+K
   50 W(IA)=ONE
C
C     IF NECESSARY INCREASE THE DIAGONAL ELEMENTS OF G.
C
      IF (.NOT. LQL) GOTO 165
      DIAG=ZERO
      DO 60 I=1,N
      ID=IWD+I
      W(ID)=G(I,I)
      DIAG=DMAX1(DIAG,VSMALL-W(ID))
      IF (I .EQ. N) GOTO 60
      II=I+1
      DO 55 J=II,N
      GA=-DMIN1(W(ID),G(J,J))
      GB=DABS(W(ID)-G(J,J))+DABS(G(I,J))
      IF (GB .GT. ZERO) GA=GA+G(I,J)**2/GB
   55 DIAG=DMAX1(DIAG,GA)
   60 CONTINUE
      IF (DIAG .LE. ZERO) GOTO 90
   70 DIAG=DIAGR*DIAG
      DO 80 I=1,N
      ID=IWD+I
   80 G(I,I)=DIAG+W(ID)
C
C     FORM THE CHOLESKY FACTORISATION OF G. THE TRANSPOSE
C     OF THE FACTOR WILL BE PLACED IN THE R-PARTITION OF W.
C
   90 IR=IWR
      DO 130 J=1,N
      IRA=IWR
      IRB=IR+1
      DO 120 I=1,J
      TEMP=G(I,J)
      IF (I .EQ. 1) GOTO 110
      DO 100 K=IRB,IR
      IRA=IRA+1
  100 TEMP=TEMP-W(K)*W(IRA)
  110 IR=IR+1
      IRA=IRA+1
      IF (I .LT. J) W(IR)=TEMP/W(IRA)
  120 CONTINUE
      IF (TEMP .LT. VSMALL) GOTO 140
  130 W(IR)=DSQRT(TEMP)
      GOTO 170
C
C     INCREASE FURTHER THE DIAGONAL ELEMENT OF G.
C
  140 W(J)=ONE
      SUMX=ONE
      K=J
  150 SUM=ZERO
      IRA=IR-1
      DO 160 I=K,J
      SUM=SUM-W(IRA)*W(I)
  160 IRA=IRA+I
      IR=IR-K
      K=K-1
      W(K)=SUM/W(IR)
      SUMX=SUMX+W(K)**2
      IF (K .GE. 2) GOTO 150
      DIAG=DIAG+VSMALL-TEMP/SUMX
      GOTO 70
C
C     STORE THE CHOLESKY FACTORISATION IN THE R-PARTITION
C     OF W.
C
  165 IR=IWR
      DO 166 I=1,N
      DO 166 J=1,I
      IR=IR+1
  166 W(IR)=G(J,I)
C
C     SET Z THE INVERSE OF THE MATRIX IN R.
C
  170 NM=N-1
      DO 220 I=1,N
      IZ=IWZ+I
      IF (I .EQ. 1) GOTO 190
      DO 180 J=2,I
      W(IZ)=ZERO
  180 IZ=IZ+N
  190 IR=IWR+(I+I*I)/2
      W(IZ)=ONE/W(IR)
      IF (I .EQ. N) GOTO 220
      IZA=IZ
      DO 210 J=I,NM
      IR=IR+I
      SUM=ZERO
      DO 200 K=IZA,IZ,N
      SUM=SUM+W(K)*W(IR)
  200 IR=IR+1
      IZ=IZ+N
  210 W(IZ)=-SUM/W(IR)
  220 CONTINUE
C
C     SET THE INITIAL VALUES OF SOME VARIABLES.
C     ITERC COUNTS THE NUMBER OF ITERATIONS.
C     ITREF IS SET TO ONE WHEN ITERATIVE REFINEMENT IS REQUIRED.
C     JFINC INDICATES WHEN TO TEST FOR AN INCREASE IN F.
C
      ITERC=1
      ITREF=0
      JFINC=-KFINC
C
C     SET X TO ZERO AND SET THE CORRESPONDING RESIDUALS OF THE
C     KUHN-TUCKER CONDITIONS.
C
  230 IFLAG=1
      IWS=IWW-N
      DO 240 I=1,N
      X(I)=ZERO
      IW=IWW+I
      W(IW)=GRAD(I)
      IF (I .GT. NACT) GOTO 240
      W(I)=ZERO
      IS=IWS+I
      K=IACT(I)
      IF (K .LE. M) GOTO 235
      IF (K .GT. MN) GOTO 234
      K1=K-M
      W(IS)=XL(K1)
      GOTO 240
  234 K1=K-MN
      W(IS)=-XU(K1)
      GOTO 240
  235 W(IS)=B(K)
  240 CONTINUE
      XMAG=ZERO
      VFACT=1.D+0
      IF (NACT) 340,340,280
C
C     SET THE RESIDUALS OF THE KUHN-TUCKER CONDITIONS FOR GENERAL X.
C
  250 IFLAG=2
      IWS=IWW-N
      DO 260 I=1,N
      IW=IWW+I
      W(IW)=GRAD(I)
      IF (LQL) GOTO 259
      ID=IWD+I
      W(ID)=ZERO
      DO 251 J=I,N
  251 W(ID)=W(ID)+G(I,J)*X(J)
      DO 252 J=1,I
      ID=IWD+J
  252 W(IW)=W(IW)+G(J,I)*W(ID)
      GOTO 260
  259 DO 261 J=1,N
  261 W(IW)=W(IW)+G(I,J)*X(J)
  260 CONTINUE
      IF (NACT .EQ. 0) GOTO 340
      DO 270 K=1,NACT
      KK=IACT(K)
      IS=IWS+K
      IF (KK .GT. M) GOTO 265
      W(IS)=B(KK)
      DO 264 I=1,N
      IW=IWW+I
      W(IW)=W(IW)-W(K)*A(KK,I)
  264 W(IS)=W(IS)-X(I)*A(KK,I)
      GOTO 270
  265 IF (KK .GT. MN) GOTO 266
      K1=KK-M
      IW=IWW+K1
      W(IW)=W(IW)-W(K)
      W(IS)=XL(K1)-X(K1)
      GOTO 270
  266 K1=KK-MN
      IW=IWW+K1
      W(IW)=W(IW)+W(K)
      W(IS)=-XU(K1)+X(K1)
  270 CONTINUE
C
C     PRE-MULTIPLY THE VECTOR IN THE S-PARTITION OF W BY THE
C     INVERS OF R TRANSPOSE.
C
  280 IR=IWR
      IP=IWW+1
      IPP=IWW+N
      IL=IWS+1
      IU=IWS+NACT
      DO 310 I=IL,IU
      SUM=ZERO
      IF (I .EQ. IL) GOTO 300
      JU=I-1
      DO 290 J=IL,JU
      IR=IR+1
  290 SUM=SUM+W(IR)*W(J)
  300 IR=IR+1
  310 W(I)=(W(I)-SUM)/W(IR)
C
C     SHIFT X TO SATISFY THE ACTIVE CONSTRAINTS AND MAKE THE
C     CORRESPONDING CHANGE TO THE GRADIENT RESIDUALS.
C
      DO 330 I=1,N
      IZ=IWZ+I
      SUM=ZERO
      DO 320 J=IL,IU
      SUM=SUM+W(J)*W(IZ)
  320 IZ=IZ+N
      X(I)=X(I)+SUM
      IF (LQL) GOTO 329
      ID=IWD+I
      W(ID)=ZERO
      DO 321 J=I,N
  321 W(ID)=W(ID)+G(I,J)*SUM
      IW=IWW+I
      DO 322 J=1,I
      ID=IWD+J
  322 W(IW)=W(IW)+G(J,I)*W(ID)
      GOTO 330
  329 DO 331 J=1,N
      IW=IWW+J
  331 W(IW)=W(IW)+SUM*G(I,J)
  330 CONTINUE
C
C     FORM THE SCALAR PRODUCT OF THE CURRENT GRADIENT RESIDUALS
C     WITH EACH COLUMN OF Z.
C
  340 KFLAG=1
      GOTO 930
  350 IF (NACT .EQ. N) GOTO 380
C
C     SHIFT X SO THAT IT SATISFIES THE REMAINING KUHN-TUCKER
C     CONDITIONS.
C
      IL=IWS+NACT+1
      IZA=IWZ+NACT*N
      DO 370 I=1,N
      SUM=ZERO
      IZ=IZA+I
      DO 360 J=IL,IWW
      SUM=SUM+W(IZ)*W(J)
  360 IZ=IZ+N
  370 X(I)=X(I)-SUM
      INFO=0
      IF (NACT .EQ. 0) GOTO 410
C
C     UPDATE THE LAGRANGE MULTIPLIERS.
C
  380 LFLAG=3
      GOTO 740
  390 DO 400 K=1,NACT
      IW=IWW+K
  400 W(K)=W(K)+W(IW)
C
C     REVISE THE VALUES OF XMAG.
C     BRANCH IF ITERATIVE REFINEMENT IS REQUIRED.
C
  410 JFLAG=1
      GOTO 910
  420 IF (IFLAG .EQ. ITREF) GOTO 250
C
C     DELETE A CONSTRAINT IF A LAGRANGE MULTIPLIER OF AN
C     INEQUALITY CONSTRAINT IS NEGATIVE.
C
      KDROP=0
      GOTO 440
  430 KDROP=KDROP+1
      IF (W(KDROP) .GE. ZERO) GOTO 440
      IF (IACT(KDROP) .LE. MEQ) GOTO 440
      NU=NACT
      MFLAG=1
      GOTO 800
  440 IF (KDROP .LT. NACT) GOTO 430
C
C     SEEK THE GREATEAST NORMALISED CONSTRAINT VIOLATION, DISREGARDING
C     ANY THAT MAY BE DUE TO COMPUTER ROUNDING ERRORS.
C
  450 CVMAX=ZERO
      IF (M .LE. 0) GOTO 481
      DO 480 K=1,M
      IA=IWA+K
      IF (W(IA) .LE. ZERO) GOTO 480
      SUM=-B(K)
      DO 460 I=1,N
  460 SUM=SUM+X(I)*A(K,I)
      SUMX=-SUM*W(IA)
      IF (K .LE. MEQ) SUMX=DABS(SUMX)
      IF (SUMX .LE. CVMAX) GOTO 480
      TEMP=DABS(B(K))
      DO 470 I=1,N
  470 TEMP=TEMP+DABS(X(I)*A(K,I))
      TEMPA=TEMP+DABS(SUM)
      IF (TEMPA .LE. TEMP) GOTO 480
      TEMP=TEMP+ONHA*DABS(SUM)
      IF (TEMP .LE. TEMPA) GOTO 480
      CVMAX=SUMX
      RES=SUM
      KNEXT=K
  480 CONTINUE
  481 DO 485 K=1,N
      LOWER=.TRUE.
      IA=IWA+M+K
      IF (W(IA) .LE. ZERO) GOTO 485
      SUM=XL(K)-X(K)
      IF (SUM) 482,485,483
  482 SUM=X(K)-XU(K)
      LOWER=.FALSE.
  483 IF (SUM .LE. CVMAX) GOTO 485
      CVMAX=SUM
      RES=-SUM
      KNEXT=K+M
      IF (LOWER) GOTO 485
      KNEXT=K+MN
  485 CONTINUE
C
C     TEST FOR CONVERGENCE
C
      INFO=0
      IF (CVMAX .LE. VSMALL) GOTO 700
C
C     RETURN IF, DUE TO ROUNDING ERRORS, THE ACTUAL CHANGE IN
C     X MAY NOT INCREASE THE OBJECTIVE FUNCTION
C
      JFINC=JFINC+1
      IF (JFINC .EQ. 0) GOTO 510
      IF (JFINC .NE. IFINC) GOTO 530
      FDIFF=ZERO
      FDIFFA=ZERO
      DO 500 I=1,N
      SUM=TWO*GRAD(I)
      SUMX=DABS(SUM)
      IF (LQL) GOTO 489
      ID=IWD+I
      W(ID)=ZERO
      DO 486 J=I,N
      IX=IWX+J
  486 W(ID)=W(ID)+G(I,J)*(W(IX)+X(J))
      DO 487 J=1,I
      ID=IWD+J
      TEMP=G(J,I)*W(ID)
      SUM=SUM+TEMP
  487 SUMX=SUMX+DABS(TEMP)
      GOTO 495
  489 DO 490 J=1,N
      IX=IWX+J
      TEMP=G(I,J)*(W(IX)+X(J))
      SUM=SUM+TEMP
  490 SUMX=SUMX+DABS(TEMP)
  495 IX=IWX+I
      FDIFF=FDIFF+SUM*(X(I)-W(IX))
  500 FDIFFA=FDIFFA+SUMX*DABS(X(I)-W(IX))
      INFO=2
      SUM=FDIFFA+FDIFF
      IF (SUM .LE. FDIFFA) GOTO 700
      TEMP=FDIFFA+ONHA*FDIFF
      IF (TEMP .LE. SUM) GOTO 700
      JFINC=0
      INFO=0
  510 DO 520 I=1,N
      IX=IWX+I
  520 W(IX)=X(I)
C
C     FORM THE SCALAR PRODUCT OF THE NEW CONSTRAINT NORMAL WITH EACH
C     COLUMN OF Z. PARNEW WILL BECOME THE LAGRANGE MULTIPLIER OF
C     THE NEW CONSTRAINT.
C
  530 ITERC=ITERC+1
      IF (ITERC.LE.MAXIT) GOTO 531
      INFO=1
      GOTO 710
  531 CONTINUE
      IWS=IWR+(NACT+NACT*NACT)/2
      IF (KNEXT .GT. M) GOTO 541
      DO 540 I=1,N
      IW=IWW+I
  540 W(IW)=A(KNEXT,I)
      GOTO 549
  541 DO 542 I=1,N
      IW=IWW+I
  542 W(IW)=ZERO
      K1=KNEXT-M
      IF (K1 .GT. N) GOTO 545
      IW=IWW+K1
      W(IW)=ONE
      IZ=IWZ+K1
      DO 543 I=1,N
      IS=IWS+I
      W(IS)=W(IZ)
  543 IZ=IZ+N
      GOTO 550
  545 K1=KNEXT-MN
      IW=IWW+K1
      W(IW)=-ONE
      IZ=IWZ+K1
      DO 546 I=1,N
      IS=IWS+I
      W(IS)=-W(IZ)
  546 IZ=IZ+N
      GOTO 550
  549 KFLAG=2
      GOTO 930
  550 PARNEW=ZERO
C
C     APPLY GIVENS ROTATIONS TO MAKE THE LAST (N-NACT-2) SCALAR
C     PRODUCTS EQUAL TO ZERO.
C
      IF (NACT .EQ. N) GOTO 570
      NU=N
      NFLAG=1
      GOTO 860
C
C     BRANCH IF THERE IS NO NEED TO DELETE A CONSTRAINT.
C
  560 IS=IWS+NACT
      IF (NACT .EQ. 0) GOTO 640
      SUMA=ZERO
      SUMB=ZERO
      SUMC=ZERO
      IZ=IWZ+NACT*N
      DO 563 I=1,N
      IZ=IZ+1
      IW=IWW+I
      SUMA=SUMA+W(IW)*W(IZ)
      SUMB=SUMB+DABS(W(IW)*W(IZ))
  563 SUMC=SUMC+W(IZ)**2
      TEMP=SUMB+.1D+0*DABS(SUMA)
      TEMPA=SUMB+.2D+0*DABS(SUMA)
      IF (TEMP .LE. SUMB) GOTO 570
      IF (TEMPA .LE. TEMP) GOTO 570
      IF (SUMB .GT. VSMALL) GOTO 5
      GOTO 570
    5 SUMC=DSQRT(SUMC)
      IA=IWA+KNEXT
      IF (KNEXT .LE. M) SUMC=SUMC/W(IA)
      TEMP=SUMC+.1D+0*DABS(SUMA)
      TEMPA=SUMC+.2D+0*DABS(SUMA)
      IF (TEMP .LE. SUMC) GOTO 567
      IF (TEMPA .LE. TEMP) GOTO 567
      GOTO 640
C
C     CALCULATE THE MULTIPLIERS FOR THE NEW CONSTRAINT NORMAL
C     EXPRESSED IN TERMS OF THE ACTIVE CONSTRAINT NORMALS.
C     THEN WORK OUT WHICH CONTRAINT TO DROP.
C
  567 LFLAG=4
      GOTO 740
  570 LFLAG=1
      GOTO 740
C
C     COMPLETE THE TEST FOR LINEARLY DEPENDENT CONSTRAINTS.
C
  571 IF (KNEXT .GT. M) GOTO 574
      DO 573 I=1,N
      SUMA=A(KNEXT,I)
      SUMB=DABS(SUMA)
      IF (NACT.EQ.0) GOTO 581
      DO 572 K=1,NACT
      KK=IACT(K)
      IF (KK.LE.M) GOTO 568
      KK=KK-M
      TEMP=ZERO
      IF (KK.EQ.I) TEMP=W(IWW+KK)
      KK=KK-N
      IF (KK.EQ.I) TEMP=-W(IWW+KK)
      GOTO 569
  568 CONTINUE
      IW=IWW+K
      TEMP=W(IW)*A(KK,I)
  569 CONTINUE
      SUMA=SUMA-TEMP
  572 SUMB=SUMB+DABS(TEMP)
  581 IF (SUMA .LE. VSMALL) GOTO 573
      TEMP=SUMB+.1D+0*DABS(SUMA)
      TEMPA=SUMB+.2D+0*DABS(SUMA)
      IF (TEMP .LE. SUMB) GOTO 573
      IF (TEMPA .LE. TEMP) GOTO 573
      GOTO 630
  573 CONTINUE
      LFLAG=1
      GOTO 775
  574 K1=KNEXT-M
      IF (K1 .GT. N) K1=K1-N
      DO 578 I=1,N
      SUMA=ZERO
      IF (I .NE. K1) GOTO 575
      SUMA=ONE
      IF (KNEXT .GT. MN) SUMA=-ONE
  575 SUMB=DABS(SUMA)
      IF (NACT.EQ.0) GOTO 582
      DO 577 K=1,NACT
      KK=IACT(K)
      IF (KK .LE. M) GOTO 579
      KK=KK-M
      TEMP=ZERO
      IF (KK.EQ.I) TEMP=W(IWW+KK)
      KK=KK-N
      IF (KK.EQ.I) TEMP=-W(IWW+KK)
      GOTO 576
  579 IW=IWW+K
      TEMP=W(IW)*A(KK,I)
  576 SUMA=SUMA-TEMP
  577 SUMB=SUMB+DABS(TEMP)
  582 TEMP=SUMB+.1D+0*DABS(SUMA)
      TEMPA=SUMB+.2D+0*DABS(SUMA)
      IF (TEMP .LE. SUMB) GOTO 578
      IF (TEMPA .LE. TEMP) GOTO 578
      GOTO 630
  578 CONTINUE
      LFLAG=1
      GOTO 775
C
C     BRANCH IF THE CONTRAINTS ARE INCONSISTENT.
C
  580 INFO=-KNEXT
      IF (KDROP .EQ. 0) GOTO 700
      PARINC=RATIO
      PARNEW=PARINC
C
C     REVISE THE LAGRANGE MULTIPLIERS OF THE ACTIVE CONSTRAINTS.
C
  590 IF (NACT.EQ.0) GOTO 601
      DO 600 K=1,NACT
      IW=IWW+K
      W(K)=W(K)-PARINC*W(IW)
      IF (IACT(K) .GT. MEQ) W(K)=DMAX1(ZERO,W(K))
  600 CONTINUE
  601 IF (KDROP .EQ. 0) GOTO 680
C
C     DELETE THE CONSTRAINT TO BE DROPPED.
C     SHIFT THE VECTOR OF SCALAR PRODUCTS.
C     THEN, IF APPROPRIATE, MAKE ONE MORE SCALAR PRODUCT ZERO.
C
      NU=NACT+1
      MFLAG=2
      GOTO 800
  610 IWS=IWS-NACT-1
      NU=MIN0(N,NU)
      DO 620 I=1,NU
      IS=IWS+I
      J=IS+NACT
  620 W(IS)=W(J+1)
      NFLAG=2
      GOTO 860
C
C     CALCULATE THE STEP TO THE VIOLATED CONSTRAINT.
C
  630 IS=IWS+NACT
  640 SUMY=W(IS+1)
      STEP=-RES/SUMY
      PARINC=STEP/SUMY
      IF (NACT .EQ. 0) GOTO 660
C
C     CALCULATE THE CHANGES TO THE LAGRANGE MULTIPLIERS, AND REDUCE
C     THE STEP ALONG THE NEW SEARCH DIRECTION IF NECESSARY.
C
      LFLAG=2
      GOTO 740
  650 IF (KDROP .EQ. 0) GOTO 660
      TEMP=ONE-RATIO/PARINC
      IF (TEMP .LE. ZERO) KDROP=0
      IF (KDROP .EQ. 0) GOTO 660
      STEP=RATIO*SUMY
      PARINC=RATIO
      RES=TEMP*RES
C
C     UPDATE X AND THE LAGRANGE MULTIPIERS.
C     DROP A CONSTRAINT IF THE FULL STEP IS NOT TAKEN.
C
  660 IWY=IWZ+NACT*N
      DO 670 I=1,N
      IY=IWY+I
  670 X(I)=X(I)+STEP*W(IY)
      PARNEW=PARNEW+PARINC
      IF (NACT .GE. 1) GOTO 590
C
C     ADD THE NEW CONSTRAINT TO THE ACTIVE SET.
C
  680 NACT=NACT+1
      W(NACT)=PARNEW
      IACT(NACT)=KNEXT
      IA=IWA+KNEXT
      IF (KNEXT .GT. MN) IA=IA-N
      W(IA)=-W(IA)
C
C     ESTIMATE THE MAGNITUDE OF X. THEN BEGIN A NEW ITERATION,
C     RE-INITILISING X IF THIS MAGNITUDE IS SMALL.
C
      JFLAG=2
      GOTO 910
  690 IF (SUM .LT. (XMAGR*XMAG)) GOTO 230
      IF (ITREF) 450,450,250
C
C     INITIATE ITERATIVE REFINEMENT IF IT HAS NOT YET BEEN USED,
C     OR RETURN AFTER RESTORING THE DIAGONAL ELEMENTS OF G.
C
  700 IF (ITERC .EQ. 0) GOTO 710
      ITREF=ITREF+1
      JFINC=-1
      IF (ITREF .EQ. 1) GOTO 250
  710 IF (.NOT. LQL) RETURN
      DO 720 I=1,N
      ID=IWD+I
  720 G(I,I)=W(ID)
  730 RETURN
C
C
C     THE REMAINIG INSTRUCTIONS ARE USED AS SUBROUTINES.
C
C
C********************************************************************
C
C
C     CALCULATE THE LAGRANGE MULTIPLIERS BY PRE-MULTIPLYING THE
C     VECTOR IN THE S-PARTITION OF W BY THE INVERSE OF R.
C
  740 IR=IWR+(NACT+NACT*NACT)/2
      I=NACT
      SUM=ZERO
      GOTO 770
  750 IRA=IR-1
      SUM=ZERO
      IF (NACT.EQ.0) GOTO 761
      DO 760 J=I,NACT
      IW=IWW+J
      SUM=SUM+W(IRA)*W(IW)
  760 IRA=IRA+J
  761 IR=IR-I
      I=I-1
  770 IW=IWW+I
      IS=IWS+I
      W(IW)=(W(IS)-SUM)/W(IR)
      IF (I .GT. 1) GOTO 750
      IF (LFLAG .EQ. 3) GOTO 390
      IF (LFLAG .EQ. 4) GOTO 571
C
C     CALCULATE THE NEXT CONSTRAINT TO DROP.
C
  775 IP=IWW+1
      IPP=IWW+NACT
      KDROP=0
      IF (NACT.EQ.0) GOTO 791
      DO 790 K=1,NACT
      IF (IACT(K) .LE. MEQ) GOTO 790
      IW=IWW+K
      IF ((RES*W(IW)) .GE. ZERO) GOTO 790
      TEMP=W(K)/W(IW)
      IF (KDROP .EQ. 0) GOTO 780
      IF (DABS(TEMP) .GE. DABS(RATIO)) GOTO 790
  780 KDROP=K
      RATIO=TEMP
  790 CONTINUE
  791 GOTO (580,650), LFLAG
C
C
C********************************************************************
C
C
C     DROP THE CONSTRAINT IN POSITION KDROP IN THE ACTIVE SET.
C
  800 IA=IWA+IACT(KDROP)
      IF (IACT(KDROP) .GT. MN) IA=IA-N
      W(IA)=-W(IA)
      IF (KDROP .EQ. NACT) GOTO 850
C
C     SET SOME INDICES AND CALCULATE THE ELEMENTS OF THE NEXT
C     GIVENS ROTATION.
C
      IZ=IWZ+KDROP*N
      IR=IWR+(KDROP+KDROP*KDROP)/2
  810 IRA=IR
      IR=IR+KDROP+1
      TEMP=DMAX1(DABS(W(IR-1)),DABS(W(IR)))
      SUM=TEMP*DSQRT((W(IR-1)/TEMP)**2+(W(IR)/TEMP)**2)
      GA=W(IR-1)/SUM
      GB=W(IR)/SUM
C
C     EXCHANGE THE COLUMNS OF R.
C
      DO 820 I=1,KDROP
      IRA=IRA+1
      J=IRA-KDROP
      TEMP=W(IRA)
      W(IRA)=W(J)
  820 W(J)=TEMP
      W(IR)=ZERO
C
C     APPLY THE ROTATION TO THE ROWS OF R.
C
      W(J)=SUM
      KDROP=KDROP+1
      DO 830 I=KDROP,NU
      TEMP=GA*W(IRA)+GB*W(IRA+1)
      W(IRA+1)=GA*W(IRA+1)-GB*W(IRA)
      W(IRA)=TEMP
  830 IRA=IRA+I
C
C     APPLY THE ROTATION TO THE COLUMNS OF Z.
C
      DO 840 I=1,N
      IZ=IZ+1
      J=IZ-N
      TEMP=GA*W(J)+GB*W(IZ)
      W(IZ)=GA*W(IZ)-GB*W(J)
  840 W(J)=TEMP
C
C     REVISE IACT AND THE LAGRANGE MULTIPLIERS.
C
      IACT(KDROP-1)=IACT(KDROP)
      W(KDROP-1)=W(KDROP)
      IF (KDROP .LT. NACT) GOTO 810
  850 NACT=NACT-1
      GOTO (250,610), MFLAG
C
C
C********************************************************************
C
C
C     APPLY GIVENS ROTATION TO REDUCE SOME OF THE SCALAR
C     PRODUCTS IN THE S-PARTITION OF W TO ZERO.
C
  860 IZ=IWZ+NU*N
  870 IZ=IZ-N
  880 IS=IWS+NU
      NU=NU-1
      IF (NU .EQ. NACT) GOTO 900
      IF (W(IS) .EQ. ZERO) GOTO 870
      TEMP=DMAX1(DABS(W(IS-1)),DABS(W(IS)))
      SUM=TEMP*DSQRT((W(IS-1)/TEMP)**2+(W(IS)/TEMP)**2)
      GA=W(IS-1)/SUM
      GB=W(IS)/SUM
      W(IS-1)=SUM
      DO 890 I=1,N
      K=IZ+N
      TEMP=GA*W(IZ)+GB*W(K)
      W(K)=GA*W(K)-GB*W(IZ)
      W(IZ)=TEMP
  890 IZ=IZ-1
      GOTO 880
  900 GOTO (560,630), NFLAG
C
C
C********************************************************************
C
C
C     CALCULATE THE MAGNITUDE OF X AN REVISE XMAG.
C
  910 SUM=ZERO
      DO 920 I=1,N
      SUM=SUM+DABS(X(I))*VFACT*(DABS(GRAD(I))+DABS(G(I,I)*X(I)))
      IF (LQL) GOTO 920
      IF (SUM .LT. 1.D-30) GOTO 920
      VFACT=1.D-10*VFACT
      SUM=1.D-10*SUM
      XMAG=1.D-10*XMAG
  920 CONTINUE
  925 XMAG=DMAX1(XMAG,SUM)
      GOTO (420,690), JFLAG
C
C
C********************************************************************
C
C
C     PRE-MULTIPLY THE VECTOR IN THE W-PARTITION OF W BY Z TRANSPOSE.
C
  930 JL=IWW+1
      IZ=IWZ
      DO 940 I=1,N
      IS=IWS+I
      W(IS)=ZERO
      IWWN=IWW+N
      DO 940 J=JL,IWWN
      IZ=IZ+1
  940 W(IS)=W(IS)+W(IZ)*W(J)
      GOTO (350,550), KFLAG
      RETURN
      END
