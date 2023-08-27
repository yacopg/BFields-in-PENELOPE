!*******************************************************************
!*                          TALLY                                  *
!*             EMF - uniform EM fields initiation                  *
!*                                                                 *
!* Short description:                                              *
!*   Allow user to input EM file into penEasy                      *
!*   This file must be included in penEasy.F additionally this     *
!*   subroutine must be called in the initally subroutine in       *
!*   penEasy.F                                                     *
!*   e.g: call EMFinitally(EMFactive)                              *
!*                                                                 *        
!*                                                                 *
!* Rev:                                                            *
!*   2022-10-05                                                    *
!*******************************************************************

      subroutine EMFinitally(active)
!*******************************************************************
!*    Initializes. To be called when reading input file            *
!*                                                                 *
!*    Output:                                                      *
!*      active: TRUE if the tally is ON.                           *
!*                                                                 *
!*******************************************************************
      use EMFmod

      logical active
      character*(*) secid,eos
      parameter (secid ='[SECTION EMF]')
      parameter (eos ='[END OF EMF SECTION]')
      character*80 buffer, upcase
      integer errorcode
      real*8 b
      integer i, j, k, l
      character*80 emfname
      integer ufile, finduf, error
      real*8 x0, y0, z0, Dx, Dy, Dz
      integer Nx, Ny, Nz

! Find Section EMF in penEasy.
      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'EMF: ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

! Bookkeeping for if section is ON or OFF.
      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        BFieldFlag = .false.

        write(*, '(a)')
     &    '>>>> EBEMF is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'EBEMF: ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'EBEMF: ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop !terminates execution of the program.
      endif

! penEasy.in File format....
![SECTION EMF]
!  ON                             STATUS (ON or OFF)
!  Bfilename                      BF filename
!  0.0 0.0 0.0                    EX, EY, EZ
! [END OF EMF SECTION]

! magnetic field file format.
c x0, y0, z0 (meters)
c Dx, Dy, Dz (meters)
c Nx, Ny, Nz
c Bx(yi = 1, xi = 1, zi = 1)
c Bx(yi = 2, xi = 1, zi = 1)
c ...
c Bx(yi = Ny, xi = 1, zi = 1)
c Bx(yi = 1, xi = 2, zi = 1)
c ...
c Bx(yi = 1, xi = Nx, zi = 1)
c Bx(yi = 1, xi = 1, zi = 2)
c ...
c Bx(yi = 1, xi = 1, zi = Nz)
c By(yi = 1, xi = 1, zi = 1)
c ...
c By(yi = Ny, xi = 1, zi = 1)
c ... ect...
c.... Note that y index is first because Matlab is default column major.
c.... We fix that here to avoid confusion.

      ! Read filename
      write(*, '(a)') "Bfield Filename:"
      read(*, *) emfname
      write(*, '(a)') emfname

      write(*,'(a)') 'EFX, EFY, EFZ:'
      read(*,*) EFX, EFY, EFZ
      write(*,'(3(1x,es12.5),1x)') EFX, EFY, EFZ

      ! Read Field Vectors from txt file into a common LookUpTable.
      write(*,'(a)') 'Starting to read BField file...'
      ufile = finduf()
      open(ufile,file=emfname,status='old',iostat=error)
      if (error.ne.0) then
        write(*,'(a)') 'EMFinisrc:ERROR: cannot open the EMF file'
        stop
      endif

      ! B field format
      write(*, '(a)') 'y0, x0, z0 (from Bfile in m converted to cm)'
      read(ufile, *)  y0, x0, z0
      ! Convert to cm
      x0 = x0*100
      y0 = y0*100
      z0 = z0*100
      write(*, '(3(1x, es12.2), 1x)') y0, x0, z0

      write(*, '(a)') 'dy, dx, dz (from Bfile in m converted to cm)'
      read(ufile, *) Dy, Dx, Dz
      Dx = Dx*100
      Dy = Dy*100
      Dz = Dz*100
      write(*, '(3(1x, es12.2), 1x)') Dy, Dx, Dz

      write(*, '(a)') 'Ny, Nx, Nz'
      read(ufile, *) Ny, Nx, Nz
      write(*, '(3(1x, I15), 1x)') Ny, Nx, Nz

      write(*, '(a)') 'Reading in Magnetic field strength, assuming T'
      write(*, '(a)') 'Converting to Gauss'

      ! Build common variable X, Y, Z, BFX, BFY, BFZ
      ALLOCATE(X_arr(Nx))
      ALLOCATE(Y_arr(Ny))
      ALLOCATE(Z_arr(Nz))
      allocate(BFX(Nx,Ny,Nz))
      allocate(BFY(Nx,Ny,Nz))
      allocate(BFZ(Nx,Ny,Nz), stat = errorcode)

      if (errorcode.ne.0) then
          write(*, '(a)') "BField:ERROR: not enough memory."
          stop
      endif

      ! Build X, Y, Z arrays
      do i = 1, Nx
        X_arr(i) = (i-1)*Dx + x0
      end do

      do i = 1, Ny
        Y_arr(i) = (i-1)*Dy + y0
      end do

      do i = 1, Nz
        Z_arr(i) = (i-1)*Dz + z0
      end do

      ! Initialize common variables that help with efficiency of the interpolation.
      ilox = 1
      iloy = 1
      iloz = 1

      ! Read in the BX, BY, BZ
      ! Should throw checks to make sure that entire file was read.
      ! that size of arrays is what is expected.
      BF_loop: do l = 1, 3
          Z_loop: do k = 1, Nz
              X_loop: do i = 1, Nx
                  Y_loop: do j = 1, Ny
                              read(ufile, *, iostat = errorcode) b
                              if (errorcode.eq.0) then
                                  !write(*, *) b
                                  !fixing the order, such that BFX(xi, yj, xk),
                                  !this is why the y_loop is done first.
                                  !Converting from T to Gauss 1T = 1e4 G
                                  if (l.eq.1) BFX(i, j, k) = b*10000.0
                                  if (l.eq.2) BFY(i, j, k) = b*10000.0
                                  if (l.eq.3) BFZ(i, j, k) = b*10000.0
                              else
                                  write(*,*) errorcode
                                  exit
                              endif
                          end do Y_loop
                      end do X_loop
                  end do Z_loop
              end do BF_loop

      ! Close Bfield file
      close(ufile)

      ! Check integrity:
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'EMF: ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif



      ! Sim restart file:
      call EMFdump(0)

      ! Set BFieldFlag to True.
      BFieldFlag = .true.


      write(*,'(a)') '>>>> EMF initialization finished >>>>'
      end


      subroutine EMFdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode: 1 to write dump file, else to read from it.          *
!*******************************************************************
      use EMFmod
      use dumpMod
      implicit none
      integer mode
      ! dump file not configured to work, this should throw an error. 
      ! Needs to be fixed in order to restart a simulation. 

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) 0
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) EFX,EFY,EFZ,BFX,BFY,BFZ
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
