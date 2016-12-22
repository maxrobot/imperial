!START OF THE MAIN PROGRAM
!
program navierstokes
!
  implicit none   !-->all the variables MUST be declared
!
  integer,parameter :: nx=129,ny=129,nt=2,ns=3,nf=3,mx=nf*nx,my=nf*ny
  !size of the computational domain (nx x ny) 
  !size of the exchanger (mx x my)
  !number of time step for the simulation
!
  !Declaration of variables
  real(8),dimension(nx,ny) :: uuu,vvv,rho,eee,pre,tmp,rou,rov,wz,tuu,tvv
  real(8),dimension(nx,ny) :: roe,tb1,tb2,tb3,tb4,tb5,tb6,tb7,tb8,tb9
  real(8),dimension(nx,ny) :: tba,tbb,fro,fru,frv,fre,gro,gru,grv,gre,&
        rot,eps,ftp,gtp,scp
  real(8),dimension(mx) :: xx
  real(8),dimension(my) :: yy
  real(8),dimension(mx,my) :: tf
  real(8),dimension(2,ns) :: coef
  integer :: i,j,itemp,k,n,nxm,iread,ni,nj,isave,longueur,imodulo
  real(8) :: xlx,yly,CFL,dlx,dx,xmu,xkt,um0,vm0,tm0
  real(8) :: xba,gma,chp,eta,uu0,dlt,um,vm,tm,x,y,dy
!***********************************************************************
  character(len=80) path_network,nom_network,path_files,name_file,&
        name_script,name_film,nchamp
  character(len=4) suffix
  character(len=20) nfile
!***********************************************************************
  !Name of the file for visualisation:
990 format('./data/output',I4.4)
  imodulo=500 !snapshots to be saved every imodulo time steps

  ! AB2 temporal scheme itemp=1
  ! RK3 temporal scheme itemp=2
  itemp=1


  ! Subroutine for the initialisation of the variables 
  call initl(uuu,vvv,rho,eee,pre,tmp,rou,rov,roe,nx,ny,xlx,yly, &
       xmu,xba,gma,chp,dlx,eta,eps,scp,xkt,uu0)
  !we need to define the time step
  dx=xlx/nx !mesh size in x
  dy=yly/ny !mesh sixe in y
  CFL=0.75  !CFL number for time step
  dlt=CFL*dlx
  print *,'The time step of the simulation is',dlt
  
  !Computation of the average velocity and temperature at t=0
  call average(uuu,um0,nx,ny)
  call average(vvv,vm0,nx,ny)
  call average(scp,tm0,nx,ny)
  write(*,*) 'Average values at t=0', um0,vm0,tm0


!BEGINNING OF TIME LOOP
  do n=1,nt
     if (itemp.eq.1) then   !TEMPORAL SCHEME AB2
      
        call fluxx(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,tb1,tb2,tb3,tb4, &
             tb5,tb6,tb7,tb8,tb9,tba,tbb,fro,fru,frv,fre,xlx,yly,xmu,xba,eps, &
             eta,ftp,scp,xkt)

        call adams(rho,rou,rov,roe,fro,gro,fru,gru,frv,grv,&
             fre,gre,ftp,gtp,scp,nx,ny,dlt)
        
        call udt(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,gma,chp)
        
     endif
        
     if (itemp.eq.2) then !TEMPORAL SCHEME RK3
        
        !loop for sub-time steps
        do k=1,ns

           call fluxx(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,tb1,tb2,tb3,tb4,&
                tb5,tb6,tb7,tb8,tb9,tba,tbb,fro,fru,frv,fre,xlx,yly,xmu,xba,eps,&
                eta,ftp,scp,xkt)
       
           call rkutta(rho,rou,rov,roe,fro,gro,fru,gru,frv,grv,&
                fre,gre,ftp,gtp,nx,ny,ns,dlt,coef,scp,k)
     
           call udt(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,gma,chp)
           
        enddo
     endif
     !loop for the snapshots, to be save every imodulo
     if (mod(n,imodulo).eq.0) then
        !this is design for Gnuplot but feel free to implement your
        !own code if you want to use Matlab or Paraview
        write(nfile, 990) n/imodulo
        open(21,file=nfile,form='formatted',status='unknown')
        ! write(21,*) nx,ny
        x=0.
        do i=1,mx
           xx(i)=x
           x=x+dx 
        enddo
        y=0.
        do j=1,my
           yy(j)=y
           y=y+dy
        enddo

        !computation of the vorticity
        call derix(vvv,nx,ny,tvv,xlx)
        call deriy(uuu,nx,ny,tuu,yly)
        do j=1,ny
          do i=1,nx
             wz(i,j)=tvv(i,j)-tuu(i,j)
          enddo
        enddo
        
        !using periodicity we copy the vorticity for the heat exchanger
        do ni=1,nf
          do nj=1,nf
             do j=1,ny
               do i=1,nx
                  tf(i+(ni-1)*nx,j+(nj-1)*ny)=wz(i,j)
               enddo
             enddo
          enddo
        enddo
        !this file will be used by gnuplot for visualisations
        ! open(21,file=nfile,form='formatted',status='unknown')
        do j=1,my
          do i=1,mx
             write(21,"(F10.5)",advance='no')  tf(i,j)
          enddo
          write(21,*) 
        enddo
        close (21)  

       !Computation of average values
       call average(uuu,um,nx,ny)
       call average(vvv,vm,nx,ny)
       call average(scp,tm,nx,ny)

       !we write the average values for velocity and temperature
       write(*,*) n,um,vm,tm
  !    
     endif

  enddo
  !END OF THE TIME LOOP
!
end program navierstokes
!
!END OF THE MAIN PROGRAMME