# BFields-in-PENELOPE
Implementation of Magnetic Fields into PENELOPE/penEasy Monte Carlo Codes for radiation transport.

## Preamble
The PENELOPE distribution includes a library of Fortran routines in the file `penfield.f` to allow the simulation of charged particle transport in the presence of magnetic fields. This `README` outlines the methods I took to use this library with penEasy and the modifications made to the PENELOPE codes to introduce arbitrary non-uniform magnetic fields into the code without recompilation from a user-defined text document.

I have included here the specific modifications I have made to the PENELOPE source code. I share this for your convenience, but it is provided "as is" without express or implied warranty. This is not a stand-alone project. You will need to obtain the source files for PENELOPE from the Radiation Safety Information Computational Center of the Oak Ridge National Laboratory. 

## Specifics
There are three main files which were edited or added in this framework. They are:
1. `penfield.f`: the main package which generates the electron/positron trajectories in an external static electromagnetic field. I have included a subset of the larger `penfield.f` subroutines included in the PENELOPE code that contain my edits to obtain the magnetic field at any given point in space.
2. `PenEasy.F`: an all-purpose main steering program.
3. `EMF.F`: An additional routine to initialize the magnetic field from a text file into PENELOPE/penEasy.

I have also included `WriteBFile.m`: A small MatLab script that writes a magnetic field map into a text document with the format expected by the code in `tallyEMF.F`.

Following instructions in the penEasy manual, the file `PenEasy.F` was edited such that the following lines in the subroutine `jumpx`
```
      if (isforcing) then       ! Interaction forcing is active
        call jumpf(dsmax,ds)    ! Get distance DS until next interaction
      else
        call jump(dsmax,ds)     ! Same without interaction forcing
      endif
```
are substituted with
```
      uldv = 0.01               ! maximum fractional changes in direction
      ulde = 0.01               ! energy
      ulem = 0.01               ! and fields
      call tpemf0(uldv,ulde,ulem,dsmaxem) ! routine in penfield.f
      dsmaxeff = min(dsmaxem,dsmax)
      if (isforcing) then       ! Interaction forcing is active
        call jumpf(dsmaxeff,ds) ! Get distance DS until next interaction
      else
        call jump(dsmaxeff,ds)  ! Same without interaction forcing
      endif
```
where `uldv`, `ulde` and `ulem` are declared as real double variables. The `tpemf0` subroutine in `penfield.f` calls the subroutine `GETEMF`, which returns the electric and magnetic field components at the particle's position. With this, `tpemf0` determines the maximum step length in the electromagnetic field based on the change in direction, energy, and EM field over the step. This steplength is used by `tpemf1` in its calculation of the particle's trajectory. **Additionally, in the file `penvox.F`, which handles the geometry, all occurrences of `call step(...)` were substituted by `call tpemf1(...)` with the same arguments.**

The magnetic field is introduced through a text file with an added section in the penEasy input file:
```
[SECTION EMF]
  ON                             STATUS (ON or OFF)
  Bfilename.txt                  BF filename
  0.0 0.0 0.0                    EX, EY, EZ
[END OF EMF SECTION]
```

The initialization routine `EMF.F` (included in this repository) looks for this section in the input file and reads in the magnetic field components, storing them as matrices in a global variable accessible to the `GETEMF` subroutine. The following diagram illustrates the function calls used by PENELOPE and penEasy to simulate particle tracks through one history. Significant additions to the code to implement magnetic fields are in the `GETEMF` and `EMFinitally` subroutines. The logic and functions have been colour-coded by file structure.

![pen_elope_easy_field](https://github.com/yacopg/BFields-in-PENELOPE/assets/56735216/e00abb5e-7905-464b-bc6a-c2c2f0bc7761)

There are two main components to the `GETEMF` subroutine (included) to return the magnetic field based on a point `Xp, Yp, Zp`. The first is a binary search algorithm to determine where in the grid this point lies. Next, it uses the surrounding eight values to perform a trilinear interpolation to determine the magnetic field components `Bx, By, Bz`. If the point is outside the magnetic field volume, the magnetic field is returned as 0. The code borrows heavily from the multidimensional linear interpolation module written in modern Fortran posted on [GitHub](https://github.com/jacobwilliams/finterp).

One notable limitation with the current implementation is that I've lost the ability to restart simulations, a handy feature of penEasy that can save an experiment from a power outage or an underestimation of the required statistical accuracy (J. Sempau, “PENELOPE/penEasy User Manual,” Mar. 2020). The code can be generalized to enable non-uniform electric fields but has not yet been implemented. It is also important to note that PENELOPE is not expected to yield a reliable result in the presence of a strong electric field (F. Salvat, “PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and Photon Transport,” OECD, 2019).
