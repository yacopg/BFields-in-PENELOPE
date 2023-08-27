CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2018)                                     C
C  Copyright (c) 2001-2018                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE PACKAGE PENFIELD
C  *********************************************************************
C
C                                       Francesc Salvat. December, 2018.
C
C  This subroutine package generates electron/positron trajectories in
C  an external static electromagnetic field. The material structure is
C  described by means of the geometry package PENGEOM.
C
C  The user must provide the following subroutines (see the examples at
C  the end of this file);
C     SUBROUTINE GETEMF(X,Y,Z,EX,EY,EZ,BX,BY,BZ)
C     SUBROUTINE GETSCP(X,Y,Z,PHI)
C  which deliver the electric field (EX,EY,EZ) in V/cm, the magnetic
C  field (BX,BY,BZ) in G, and the scalar potential PHI in V, at an
C  arbitrary point with coordinates (X,Y,Z) cm.
C
C  Notice that 1 G = 1.0E-4 T and 1 statV = 299.792458 V.
C
C  The sequence of calls to generate a particle track is the following,
C  1: CALL START  ! Initialises the simulation of the track.
C  2: CALL TPEMF0(ULDV,ULDE,ULEM,DSMAX)  ! Determines the maximum
C       ! step length DSMAX consistent with the adopted delta values.
C  3: CALL JUMP(DSMAX,DS)  ! Determines the path length DS to the next
C       ! interaction event.
C  4: CALL TPEMF1(DS,DSEF,NCROSS)  ! Moves the particle to the
C       ! end of the step and determines its final direction and energy.
C  5: CALL KNOCK(DE,ICOL)  ! Simulates next interaction. Interaction
C       !  energy losses are assumed to be deposited at the hinge.
C  6: GO TO 1 or 2.
C  This sequence has to be discontinued when the particle leaves the
C  material system or is absorbed.
C
C  Although the transport of photons is not affected by the field, the
C  present subroutines can also be called to generate photon histories.
C  This is useful to simplify the structure of the main program (all
C  kinds of particles can be tracked by using the same sequence of calls
C  to the simulation subroutines).
C
C  IMPORTANT NOTE: One of the practical features of PENGEOM is that
C  particles are automatically transported through void regions (i.e.
C  those with MAT=0) following straight trajectories without loosing
C  energy. Of course, in the presence of fields this would give
C  erroneous results. The easiest method to avoid this conflict is to
C  fill up the voids with a fake material of very small density.
C
C  J. Groeneveld:
C       NOTE: this is not the complete penfield package, but rather the
C           subroutines inside penfield.f that have been edited to enable
C           arbitrary magnetic fields within PENELOPE
C  *********************************************************************
C                       module EMFmod
C  *********************************************************************
      module EMFmod
C
C  This module contains the 'global' variables for the electromagentic field
C  components. Used in GETEMF, GETSCP, and EMF.F
C
      implicit none
      save
      real*8 EFX,EFY,EFZ ! constants
      real*8, allocatable :: BFX(:,:,:) ! magnetic field components.
      real*8, allocatable :: BFY(:,:,:)
      real*8, allocatable :: BFZ(:,:,:)
      real*8, allocatable :: X_arr(:) ! grid locations.
      real*8, allocatable :: Y_arr(:)
      real*8, allocatable :: Z_arr(:)
      integer ilox, iloy, iloz ! efficiency parameter for binary search.
      logical BFieldFlag ! bookkeeping

      end module

C  *********************************************************************
C                       SUBROUTINE GETEMF
C  *********************************************************************
      SUBROUTINE GETEMF(XP,YP,ZP,EX,EY,EZ,BX,BY,BZ)
C
C  This subroutine returns the electric and magnetic field components
C  at the input point (XP,YP,ZP).
C
C  Input arguments:
C  (XP,YP,ZP) ...  position coordinates (cm).
C
C  Output arguments:
C  (EX,EY,EZ) ...  components of the E-field, in V/cm.
C  (BX,BY,BZ) ...  components of the B-field, in G.
C
C  Notice that 1 G = 1.0E-4 T and 1 statV = 299.792458 V.
C
C  It is assumed that users will provide an equivalent subroutine, with
C  the same name and arguments, to define their fields.
C
C  For example:
C  For a uniform E-field: (EX,EY,EZ) three constants.
C                         (BX,BY,BZ)=(0,0,0).
C  For a uniform B-field: (EX,EY,EZ)=(0,0,0).
C                         (BX,BY,BZ) three constants.
C  For a point charge at rest at the origin of coordinates:
C                         EX=C*XP/(SQRT(XP**2+YP**2+ZP**2))**3,
C                         EY=C*YP/(SQRT(XP**2+YP**2+ZP**2))**3,
C                         EZ=C*ZP/(SQRT(XP**2+YP**2+ZP**2))**3,
C                         (BX,BY,BZ)=(0,0,0).
C                         where C is a constant.
C
C  This example is for a uniform E field, the components of the E-field
C  and of the B-field are entered through the module EMFmod.
C
C *********************************************************************
C An arbitrary magnetic field is introduced with an external text file.
C
C This subroutine uses a trilinear interpolation algorthim as described
C in wikipedia
C
      USE TRACK_mod
      use EMFmod

      implicit none
      real*8 :: XP, YP, ZP, EX, EY, EZ, BX, BY, BZ

      ! Internal variables
      integer, DIMENSION(2) :: ix, iy, iz
      INTEGER, DIMENSION(3) :: mflag
      real*8 :: p1, p2, p3  ! differences between each of x,y,z and smaller related coordinates
      real*8 :: q1, q2, q3  ! 1-q: for interpolation
      real*8 :: fx11, fx21, fx12, fx22, fxy1, fxy2 ! values between eight corners of our cube

      if (.NOT.(BFieldFlag)) THEN !No magnetic field initialized, set to zero.
        BX = 0.0
        BY = 0.0
        BZ = 0.0
        EX = 0.0
        EY = 0.0
        EZ = 0.0
        RETURN
      end if

      ! Use X, Y, Z to determine the indices
      ! mflag indicates whether a requested point is outside the map.
      ! DINTRV enables extrapolation, but might not be appropriate, or necessary.
      call DINTRV(X_arr, XP, ilox, ix(1), ix(2), mflag(1))
      call DINTRV(Y_arr, YP, iloy, iy(1), iy(2), mflag(2))
      call DINTRV(Z_arr, ZP, iloz, iz(1), iz(2), mflag(3))

C
C  ****  The field is set equal to zero in regions where MAT=0, i.e. in
C        the outer vacuum, because otherwise particles that leave the
C        system could return to it under the action of the field.
C
      IF(ALL(mflag.eq.0).and.(MAT.ne.0)) then
        !Begin a trilinear interpolation
        q1 = (XP - X_arr(ix(1)))/(X_arr(ix(2)) - X_arr(ix(1)))
        q2 = (YP - Y_arr(iy(1)))/(Y_arr(iy(2)) - Y_arr(iy(1)))
        q3 = (ZP - Z_arr(iz(1)))/(Z_arr(iz(2)) - Z_arr(iz(1)))
        p1 = 1 - q1
        p2 = 1 - q2
        p3 = 1 - q3

C*********************Interpolate in BFX******************************
        ! Interpolate along x
        fx11 = p1*BFX(ix(1),iy(1),iz(1)) + q1*BFX(ix(2),iy(1),iz(1))
        fx21 = p1*BFX(ix(1),iy(2),iz(1)) + q1*BFX(ix(2),iy(2),iz(1))
        fx12 = p1*BFX(ix(1),iy(1),iz(2)) + q1*BFX(ix(2),iy(1),iz(2))
        fx22 = p1*BFX(ix(1),iy(2),iz(2)) + q1*BFX(ix(2),iy(2),iz(2))

        ! Interpolate these values along y
        fxy1 = p2*fx11 + q2*fx21
        fxy2 = p2*fx12 + q2*fx22

        ! Interpolate these values along z
        BX = p3*fxy1 + q3*fxy2

C*********************Interpolate in BFY******************************
        ! Interpolate along x
        fx11 = p1*BFY(ix(1),iy(1),iz(1)) + q1*BFY(ix(2),iy(1),iz(1))
        fx21 = p1*BFY(ix(1),iy(2),iz(1)) + q1*BFY(ix(2),iy(2),iz(1))
        fx12 = p1*BFY(ix(1),iy(1),iz(2)) + q1*BFY(ix(2),iy(1),iz(2))
        fx22 = p1*BFY(ix(1),iy(2),iz(2)) + q1*BFY(ix(2),iy(2),iz(2))

        ! Interpolate these values along y
        fxy1 = p2*fx11 + q2*fx21
        fxy2 = p2*fx12 + q2*fx22

        ! Interpolate these values along z
        BY = p3*fxy1 + q3*fxy2

C*********************Interpolate in BFZ******************************
        ! Interpolate along x
        fx11 = p1*BFZ(ix(1),iy(1),iz(1)) + q1*BFZ(ix(2),iy(1),iz(1))
        fx21 = p1*BFZ(ix(1),iy(2),iz(1)) + q1*BFZ(ix(2),iy(2),iz(1))
        fx12 = p1*BFZ(ix(1),iy(1),iz(2)) + q1*BFZ(ix(2),iy(1),iz(2))
        fx22 = p1*BFZ(ix(1),iy(2),iz(2)) + q1*BFZ(ix(2),iy(2),iz(2))

        ! Interpolate these values along y
        fxy1 = p2*fx11 + q2*fx21
        fxy2 = p2*fx12 + q2*fx22

        ! Interpolate these values along z to calculate BZ
        BZ = p3*fxy1 + q3*fxy2

      else ! NOT(ALL(mflag.eq.0).and.(MAT.neq.0))
        BX = 0.0
        BY = 0.0
        BZ = 0.0
        EX = 0.0
        EY = 0.0
        EZ = 0.0
        RETURN
      end if

      ! Set electric field constant value in emfmod
      EX=EFX
      EY=EFY
      EZ=EFZ

      RETURN


      contains

C  ********************************************************************
C  Function to return the indices in `xt` that bound `x`, to use for interpolation.
!  If outside the range, then the indices are returned that can
!  be used for extrapolation.
!  Precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   iright=2,    mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   iright=i+1,  mflag=0
!         if   xt(n) <= x           then ileft=n-1, iright=n,    mflag=1
!```
!   Code adapted from code posted by Jacob Williams, 10/9/2019 which was itself
!   based on original code written by carl de boor and amos, d. e.in 1980...
!   https://github.com/jacobwilliams/finterp/blob/master/src/linear_interpolation_module.F90
!
      SUBROUTINE DINTRV(xt, x, ilo, ileft, iright, mflag, inearest)

        implicit none
        !note: wp in source code was a way of changing the precision
        real*8,dimension(:),intent(in)  :: xt       !! a knot or break point vector
        real*8,intent(in)               :: x        !! argument
        integer,intent(inout)           :: ilo      !! an initialization parameter which must be set
                                                    !! to 1 the first time the array `xt` is
                                                    !! processed by dintrv. `ilo` contains information for
                                                    !! efficient processing after the initial call and `ilo`
                                                    !! must not be changed by the user.  each dimension
                                                    !! requires a distinct `ilo` parameter.
        integer,intent(out)              :: ileft    !! left index
        integer,intent(out)              :: iright   !! right index
        integer,intent(out)              :: mflag    !! signals when `x` lies out of bounds
        integer,intent(out),optional     :: inearest !! nearest index

        integer :: ihi, istep, imid, n

        n = size(xt)

        if (n==1) then
        ! this is only allowed for nearest interpolation
          if (present(inearest)) then
              inearest = 1
              return
          end if
        end if

        ihi = ilo + 1
        if ( ihi>=n ) then
            if ( x >= xt(n) ) then
                if (x == xt(n)) THEN
                  mflag = 0 !don't want to flag if its not outside bounds -
                else
                  mflag = 1
                end if

                ileft = n-1
                iright= n
                if (present(inearest)) inearest = n
                return
            end if
            if ( n<=1 ) then
                mflag = -1
                ileft = 1
                iright= 2
                if (present(inearest)) inearest = 1
                return
            end if
            ilo = n - 1
            ihi = n
        endif

        if ( x>=xt(ihi) ) then

            ! now x >= xt(ilo). find upper bound
            istep = 1
            do
                ilo = ihi
                ihi = ilo + istep
                if ( ihi>=n ) then
                    if ( x>=xt(n) ) then
                        if (x == xt(n)) THEN
                          mflag = 0 !don't want to flag if its not outside bounds -
                        else
                          mflag = 1
                        end if

                        ileft = n-1
                        iright= n
                        if (present(inearest)) inearest = n
                        return
                    end if
                    ihi = n
                elseif ( x>=xt(ihi) ) then
                    istep = istep*2
                    cycle
                endif
                exit
            end do

        else

            if ( x>=xt(ilo) ) then
                mflag = 0
                ileft = ilo
                iright= ilo+1
                if (present(inearest)) then
                    if ( abs(x-xt(ileft)) <= abs(x-xt(iright)) ) then
                        inearest = ileft
                    else
                        inearest = iright
                    end if
                end if
                return
            end if
            ! now x <= xt(ihi). find lower bound
            istep = 1
            do
                ihi = ilo
                ilo = ihi - istep
                if ( ilo<=1 ) then
                    ilo = 1
                    if ( x<xt(1) ) then
                        mflag = -1
                        ileft = 1
                        iright= 2
                        if (present(inearest)) inearest = 1
                        return
                    end if
                elseif ( x<xt(ilo) ) then
                    istep = istep*2
                    cycle
                endif
                exit
            end do

        endif

        ! now xt(ilo) <= x < xt(ihi). narrow the interval
        do
            imid = (ilo+ihi)/2
            if ( imid==ilo ) then
                mflag = 0
                ileft = ilo
                iright= ilo+1
                if (present(inearest)) then
                    if ( abs(x-xt(ileft)) <= abs(x-xt(iright)) ) then
                        inearest = ileft
                    else
                        inearest = iright
                    end if
                end if
                return
            end if
            ! note. it is assumed that imid = ilo in case ihi = ilo+1
            if ( x<xt(imid) ) then
                ihi = imid
            else
                ilo = imid
            endif
        end do

      END SUBROUTINE DINTRV

      end
C  *********************************************************************
C                       SUBROUTINE GETSCP
C  *********************************************************************
      SUBROUTINE GETSCP(XP,YP,ZP,PHI)
C
C  This routine returns the scalar potential at the point (XP,YP,ZP).
C
C  Input arguments:
C  (XP,YP,ZP) ...  position coordinates (cm).
C
C  Output argument:
C  PHI ..........  scalar potential, in V.
C
C  Notice that 1 statV = 299.792458 V.
C
C  It is assumed that users will provide an equivalent subroutine, with
C  the same name and arguments, to define their fields.
C
C  For example:
C  For a uniform E-field: (EX,EY,EZ) three constants.
C                         PHI=-(EX*XP+EY*YP+EZ*ZP).
C  For a point charge at rest at the origin of coordinates:
C                         PHI=C/SQRT(XP**2+YP**2+ZP**2),
C                         where C is a constant.
C
C  This example is for a uniform E field. The components of the E-field
C  are entered through the module EMFmod.
C
      USE TRACK_mod
      USE EMFmod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
C  ****  The potential is set equal to zero in regions where MAT=0,
C        i.e. in the outer vacuum, because otherwise particles that
C        leave the system could return to it under the action of the
C        field.
C
      IF(MAT.EQ.0) THEN
        PHI=0.0D0
        RETURN
      ENDIF
C
      PHI=-(EFX*XP+EFY*YP+EFZ*ZP)
      RETURN
      END
