MODULE gort_lidar

  ! This module contains Wenge Ni-Meister's Geometric Optical and Radiative
  ! Transfer (GORT) lidar model. Modified Ping and Wenze's old code on May 19, 2010 by WNM
  
  ! TO USE:
  !
  ! call run_gort(all_gort_in, height_levels, G, gamma, fp, rdfp, pgap, dpdz, wvfm) 
  !
  !   where,
  ! Input:
  ! all_gort_in: Input for gort model. It is of type gort_input. If it is an array, then convoluted outputs
  !              will be generated. If it is a scalar variable of type gort_input, then outputs
  !              for single profile will be generated.
  ! 
  ! Output:
  ! height_levels : height levels for the convoluted or single output profile (depending on the type of input)
  ! fp            : foliage profile
  ! efp           : effective foliage profile: clumping factor * foliage profile*G 
  ! pgap          : gap probability, clumped and non-clumped.
  ! dpdz          : dp/dz, clumped and nonclumped 
  ! wvfm          : waveform, clumped and nonclumped
  ! NOTE:
  ! In the case of single profile, lower and higher end of the output are determined by h1, h2, and dz 
  ! of the gort input. In the case of multi-layer input, the model will take the lowest end of all the
  ! profiles as the low end for the convoluted profile, the high end of all the profiles as the high
  ! end for the convoluted profile, and the smalled dz of all the profiles as the delta_z for the 
  ! convoluted profile.
  ! Modified profile_params type, procedure run_single_gort, run_convolute_gort, converlute and interpolate_profile
  
  
  implicit none
  
  private
  public :: run_gort
  public :: gort_input
  
  ! type definition
  type gort_input
     real :: dens_tree         ! Crown count density (m-2)
     real :: dens_foliage      ! foliage area volume density of a single crown(m-1)
     real :: h1, h2            ! lower and upper bound of crown centers (m)
     real :: dz                ! height increase for each height level
     real :: horz_radius       ! horizontal crown radius (m)
     real :: vert_radius       ! vertical crown radius (m)
     real :: zenith            ! solar zenith angle (degree)
     real :: rhol              ! leaf reflectance
     real :: taul              ! leaf transmittance
     real :: rhos              ! background albedo 
     integer :: pft            ! plant function type
  end type gort_input
  
  real, parameter :: PI=3.14159265
  
  ! interface definition
  interface run_gort
     module procedure run_convolute_gort
     module procedure run_single_gort
  end interface
  
  ! temporary type definition
  type profile_params
     real, dimension(:), pointer :: height_levels ! height level for each profile
     real, dimension(:), pointer :: fp            ! foliage profile, no clump involved
     real, dimension(:), pointer :: efp           ! effective foliage profile:g*gamma*fp
     real, dimension(:), pointer :: pgap(:,:)     ! gap probability
     real, dimension(:), pointer :: dpdz(:,:)     ! dp/dz
     real, dimension(:), pointer :: wvfm(:,:)     ! waveform 
     real, dimension(:), pointer :: G             ! leaf orientation factor
     real, dimension(:), pointer :: gamma         ! clumping factor
  end type profile_params

contains
  
  ! --------------------------------------------------------------------
  
  subroutine run_single_gort(gin, height_levels, G, gamma, fp, efp, pgap, dpdz, wvfm) 
    
    ! Input and Output
    real:: G, gamma

    ! Input 
    type(gort_input) ::  gin  ! input structure for gort model
    
    ! Output
    real, dimension(:), pointer   :: height_levels  ! the height levels for the profile
    real, dimension(:), pointer   :: fp             ! foliage profile
    real, dimension(:), pointer   :: efp            ! effective foliage profile: G*gamma*foliage profile
    real, dimension(:), pointer   :: pgap(:,:)           ! gap probability 
    real, dimension(:), pointer   :: dpdz(:,:)           ! dp/dz
    real, dimension(:), pointer   :: wvfm(:,:)           ! waveform
    
    ! optional output
    
    ! local variables
    real :: ratio_rhocg                             ! ratio of canopy and background reflectivity  WNM 1/19/2011
    integer :: N_height_level, i,j
    
    ! ----------------------------------------------------------
    ! DEBUG
    print *, 'I am in run_single_gort ..........'
    
    ! get G, leaf orientation factor
    G = get_G(gin%zenith, gin%pft)
    
    ! Check the validity of the inputs
    call check_inputs(gin)
    
    ! calculate the ratio of canopy and background reflectivities (Equation (7) in Ni-Meister et al. 2011 paper) by WNM 1/19/2011

    ratio_rhocg=2.*gin%rhol/(3.*gin%rhos)

    ! define the vertical layers
    call get_height_level(gin, height_levels)
    N_height_level = size(height_levels)
    
    ! get variables, clump, fp and efp
    allocate(fp(N_height_level))
    allocate(efp(N_height_level))

    call get_foliage_profile(gin,height_levels,fp)

    !computer clumping factor: gamma
    gamma = get_analytical_clump(gin)

    ! get effective foliage profile
    fp = G *fp
    efp =  gamma * fp
    
    ! compute Pgap
    allocate(Pgap(N_height_level,2))
    Pgap(:,1) = T(fp, gin%dz, gin%zenith)
    Pgap(:,2) = T(efp, gin%dz, gin%zenith)

    ! compute dp/dz
    allocate(dpdz(N_height_level,2))
    allocate(wvfm(N_height_level,2))    
    do j=1,2
      dpdz(1,j)= pgap(1,j)
      dpdz(N_height_level,j)= - (pgap(N_height_level-1,j)-pgap(N_height_level,j))/gin%dz ! last layer
      do i=2, N_height_level-1
        dpdz(i,j)= -(pgap(i-1,j)-pgap(i+1,j))/(2*gin%dz)
      enddo

    !compute waveform  (Equation (3) in Ni-Meister et al., 2011) by WNM 1/19/2011
      wvfm(1,j)=dpdz(1,j)/((1-dpdz(1,j))*ratio_rhocg+dpdz(1,j))
      wvfm(1,j)=wvfm(1,j)/gin%dz

      do i=2, N_height_level-1
        wvfm(i,j)= dpdz(i,j)/((1-dpdz(1,j))+dpdz(1,j)/ratio_rhocg)
      enddo
    enddo

 end subroutine run_single_gort

 !----------------------------------------------------------------------------------------------------------

 subroutine run_convolute_gort(gin, height_levels, G, gamma, fp, efp, pgap, dpdz, wvfm) 

   ! Input
   type(gort_input), dimension(:) :: gin

   ! Output

   real, dimension(:), pointer :: G             ! leaf orientation factor
   real, dimension(:), pointer :: gamma         ! clumping factor
   real, dimension(:), pointer   :: height_levels  ! the height levels for the profile
   real, dimension(:), pointer   :: fp             ! foliage profile
   real, dimension(:), pointer   :: efp            ! effective foliage profile: G*gamma*foliage profile
   real, dimension(:), pointer   :: pgap(:,:)           ! gap probability 
   real, dimension(:), pointer   :: dpdz(:,:)           ! dp/dz
   real, dimension(:), pointer   :: wvfm(:,:)           ! waveform

   ! local variables
   integer :: num_profiles, iprofile, N_height_level, i,j
   integer :: N_convolute_height_level, ilevel, iszn
   real    :: szn, dz
   real :: ratio_rhocg,rhol_ave=0, rhos_ave=0      ! ratio of canopy and background reflectivity  WNM 1/19/2011

   type(profile_params), dimension(:), pointer :: all_convolute_input

   real, dimension(:), pointer   :: tmp_fp, tmp_efp
   real, dimension(:), pointer   :: tmp_height_levels

   real, dimension(:), pointer   :: d_fp, d_rdfp, d_rifp, d_efp, d_height_levels
   real, dimension(:,:), pointer :: convoluted_szn_transmit
   !-----------------------------------------------------------------

   ! DEBUG
   print *, 'I am in run_convolute_gort.........'

   num_profiles = size (gin)

   ! check to verify that all profiles in all_gort_in have the same zenith angle
   szn = gin(1)%zenith
   do iprofile=2, num_profiles
     if (gin(iprofile)%zenith .ne. szn) then
        print *, 'The sun zenith angle of all the participating profiles should be the same. Exit.....'
        exit
     end if
   end do

   print *, 'The sun zenith angle for the convoluted profile is ', szn

   ! Get effective efp for each profile 

   nullify(all_convolute_input)
   allocate(all_convolute_input(num_profiles))
   allocate(G(num_profiles))
   allocate(gamma(num_profiles))

   do iprofile=1, num_profiles

     call check_inputs(gin(iprofile))

     G(iprofile) = get_G(gin(iprofile)%zenith, gin(iprofile)%pft)

     ! define the vertical layers
     call get_height_level(gin(iprofile), tmp_height_levels) !tmp_height_levels is allocated in the subroutine
     N_height_level = size(tmp_height_levels)

     ! get intermediate variables, clump and fp
     allocate(tmp_fp(N_height_level))
     allocate(tmp_efp(N_height_level))
   
     ! consider putting the following two subroutines as function instead
     call get_foliage_profile(gin(iprofile), tmp_height_levels, tmp_fp) 
     gamma(iprofile) = get_analytical_clump(gin(iprofile))

     ! add clump factor into foliage profile data 
     tmp_fp  = G(iprofile) * tmp_fp
     tmp_efp = gamma(iprofile) * tmp_fp 

     !DEBUG
     !Change by ping, print *, tmp_fp(10), clumpd, tmp_efp(10)
     print *, tmp_fp(10), gamma(iprofile), tmp_efp(10)
     print *, minval(tmp_height_levels), maxval(tmp_height_levels)
     print *, 'Assinging single_gort output to convolute_input'

     ! fill in convolute_input
     all_convolute_input(iprofile)%fp => tmp_fp
     all_convolute_input(iprofile)%efp => tmp_efp
     all_convolute_input(iprofile)%height_levels => tmp_height_levels
     
     !DEBUG
     print *, all_convolute_input(iprofile)%efp(10)

     ! deallocate and null the tmp pointers
     nullify(tmp_fp)
     nullify(tmp_efp)
     nullify(tmp_height_levels)

   end do

   ! Convolute fp, efp by calling convolute, but here only fp is final output for run_gort
   ! efp will be used for getting pgap

   ! calculate the ratio of canopy and background reflectivities (Equation (7) in Ni-Meister et al. 2011 paper) 
   ! use the averaged rhol and rhos. by WNM 1/19/2011

   do iprofile=1, num_profiles
      rhol_ave=rhol_ave+gin(iprofile)%rhol
      rhos_ave=rhos_ave+gin(iprofile)%rhos
   enddo
   rhol_ave=rhol_ave/num_profiles
   rhos_ave=rhos_ave/num_profiles
   ratio_rhocg=2.*rhol_ave/(3.*rhos_ave)

   print *, 'Calling convolute ......' 
   call convolute(all_convolute_input,         &                    ! input
        fp, efp, height_levels)    ! convoluted output, also run_gort output

   ! Compute Pgap, dp/dz and waveform 
   N_convolute_height_level = size(height_levels)
   dz = height_levels(2)-height_levels(1)

   ! compute Pgap
   allocate(pgap(N_convolute_height_level,2))
   pgap(:,1) =  T(fp, dz, szn) 
   pgap(:,2) =  T(efp, dz, szn) 
    
   ! compute dp/dz
   allocate(dpdz(N_convolute_height_level,2))
   allocate(wvfm(N_convolute_height_level,2))    
   do j=1,2
     dpdz(1,j)=  pgap(1,j) ! first layer, ground return 
     dpdz(N_convolute_height_level,j)= - (pgap(N_convolute_height_level-1,j)-pgap(N_convolute_height_level,j))/dz ! last layer
     do i=2, N_convolute_height_level-1
       dpdz(i,j)=-(pgap(i-1,j)-pgap(i+1,j))/(2*dz)
     enddo
   
   !compute waveform  (Equation (7) in Ni-Meister et al., 2017) by WNM 1/19/2011
     wvfm(1,j)=dpdz(1,j)/((1-dpdz(1,j))*ratio_rhocg+dpdz(1,j))
     wvfm(1,j)=wvfm(1,j)/dz                                           !by WNM 1/19/2011

   !  do i=2, N_convolute_height_level-1
     do i=2, N_convolute_height_level
       wvfm(i,j)= dpdz(i,j)/((1-dpdz(1,j))+dpdz(1,j)/ratio_rhocg)
    enddo
   enddo

   do iprofile=1, num_profiles
     deallocate(all_convolute_input(iprofile)%fp)
     deallocate(all_convolute_input(iprofile)%efp)
     deallocate(all_convolute_input(iprofile)%height_levels)
   end do
   deallocate(all_convolute_input)

 end subroutine run_convolute_gort

! ----------------------------------------------------------------------------------------------------
 function T(fp, delta_z, szn)

   ! return a transmittance array based on the size of input fp

   real, dimension(:), pointer :: fp
   real :: delta_z, szn

   real, dimension(size(fp)) :: T
   
   ! local vars
   integer :: N_height_level, ilevel

   N_height_level = size(fp)
   do ilevel=1, N_height_level
     T(ilevel) = exp((-sum(fp(ilevel:N_height_level))*delta_z) /cos(radius(szn)))
   end do

end function 

! -------------------------------------------------------------------------------------------------
   
   
 subroutine check_inputs(gin)

   ! check the validity of inputs
   type(gort_input) :: gin                  ! gort input

   real :: dens_tree            ! Crown count density (m-2)
   real :: dens_foliage         ! foliage area volume density of a single crown(m-1)
   real :: h1, h2               ! lower and upper bound of crown centers (m)
   real :: delta_z              ! height increase for each height level
   real :: horz_radius          ! horizontal crown radius (m)
   real :: vert_radius          ! vertical crown radius (m)
   real :: zenith               ! solar zenith angle (degree)

   dens_tree = gin%dens_tree
   dens_foliage = gin%dens_foliage
   h1 = gin%h1
   h2 = gin%h2
   delta_z = gin%dz
   horz_radius = gin%horz_radius
   vert_radius = gin%vert_radius
   zenith = gin%zenith

   ! check density
   if ( (dens_tree < 0.) .or. (dens_foliage < 0.)) then
     write(*,*) 'Density must be greater than 0. STOPPING ...'
     stop
   end if

   ! check the height
   if ((h1 <= 0.) .or. (h2 <= 0.)) then
     write(*,*) 'Center of crown height must be greater than 0. STOPPING ...'
     stop
   end if

   if (h1 >= h2) then
     write(*,*) 'h1 must be less than h2. STOPPING  ...'
     stop
   end if

 !  if ( h1 < vert_radius ) then
 !    write(*,*) 'h1 must be greater than the vert_radius. STOPPING...'
 !    stop
 !  end if

   ! check radius
   if ((horz_radius <= 0.) .or. (vert_radius <= 0.)) then
     write(*,*) 'Center of crown height must be greater than 0. STOPPING ...'
     stop
   end if

    !????? not necessary????
   !if (horz_radius >= vert_radius) then
   !  write(*,*) 'horz_radius must be less than vert_radius.'
   !  write(*,*) 'STOPPING ...'
   !  stop
   !end if

   ! check angle
   if ((zenith < 0.) .or. (zenith > 90.)) then
     write (*,*) 'Solar zenith angle, zenith, is in degree, and must be within(0.,90.).'
     write (*,*) 'STOPPING ...'
     stop
   end if

 end subroutine check_inputs

! ----------------------------------------------------------
 subroutine get_height_level(gin, height_levels)

   ! compute the height levels between h1 and h2
   ! NOTE : here last layer could be less than delta_z.
   !        It may be necessary to return the dz at the last layer to
   !        ensure the accuracy of calculation, otherwise, we can make
   !        sure (h2-h1) can be integer divided by delta_z.

   type(gort_input), intent(in) :: gin
   real, dimension(:), pointer :: height_levels

   ! local
   real :: h1, h2, vert_radius, delta_z
   integer :: n_levels, i, temp1
   real :: jh1, jh2, temp2, dif

   ! -------------------------------------------


   h1 = gin%h1
   h2 = gin%h2
   vert_radius = gin%vert_radius
   delta_z = gin%dz


   jh1 = h1 - vert_radius
   jh2 = h2 + vert_radius

   ! based on wenge, calculate all levels to the ground
   jh1 = 0.0

   ! Added the rounding checking for different machines.

   temp1 = ceiling((jh2-jh1)/delta_z)
   temp2 = (jh2-jh1)/delta_z

   if ((temp1-temp2) > 0.999) then
     n_levels = ceiling((jh2-jh1)/delta_z)
   else
     n_levels = ceiling((jh2-jh1)/delta_z) + 1
   end if


   allocate(height_levels(n_levels))

   do i =1, n_levels-1
     height_levels(i) = jh1 + delta_z * (i-1)
   end do

   height_levels(n_levels) = jh2

 end subroutine get_height_level

! ----------------------------------------------------------

 subroutine get_foliage_profile(gin, height_levels, fp)

   type(gort_input) :: gin
   real, dimension(:), intent(in) :: height_levels
   real, dimension(:), intent(inout) :: fp

   ! local
   real :: h1, h2, vert_radius, tmp
   real :: z, jh11, jh12, jh21, jh22, hdif
   integer :: n_levels, i

   real, dimension(:), allocatable :: xx

   ! ---------------------------------

   h1 = gin%h1
   h2 = gin%h2
   vert_radius = gin%vert_radius

   tmp=gin%dens_tree * gin%dens_foliage * PI * (gin%horz_radius ** 2)

   jh11 = h1 - vert_radius
   if (jh11 < 0) then 
      jh11 = 0 ! updated by Wenge on 61/8/2011
   endif

   jh12 = h1 + vert_radius
   jh21 = h2 - vert_radius
   if (jh21 < 0) then 
      jh21 = 0 ! updated by Wenge on 61/8/2011
   endif
   jh22 = h2 + vert_radius
   hdif = h2 - h1

   n_levels = size(height_levels)
   allocate(xx(n_levels))

   ! define xx based on the location of each height level
   do i = 1, n_levels

     z = height_levels(i)
     if (z < jh11) xx(i)=0.0
     if (jh12 < jh21) then
        if ((z >= jh11) .and. (z <= jh12))   &
           xx(i) = height_function1(z, h1, vert_radius, hdif)
        if ((z > jh12) .and. (z <= jh21))   &
           xx(i) = height_function3(hdif, vert_radius)
        if ((z > jh21) .and. (z <= jh22))   &
           xx(i) = height_function2(z, h2, vert_radius, hdif)
     else ! (jh12 >= jh21)
        if ((z >= jh11) .and. (z <= jh21))   &
           xx(i) = height_function1(z, h1, vert_radius, hdif)
        if ((z > jh21) .and. (z <= jh12))   &
           xx(i) = height_function4(z, h1, h2, vert_radius)
        if ((z > jh12) .and. (z <= jh22))   &
           xx(i) = height_function2(z, h2, vert_radius, hdif)
     endif

     fp(i) = tmp * xx(i)

  end do

  deallocate(xx)

 end subroutine get_foliage_profile

! ----------------------------------------------------------

 subroutine convolute(all_convolute_input,                       &      ! input
        fp, efp, height_levelC)         ! convoluted output


   ! NOTE: This subroutine convolute multiple profiles into one.
   !       The algorithm used here is mapping each profile to the desired height levels
   !       and then convolute (either sum or multiply) profiles together.
   !       The height_levelC will be based upon the input for desired height, however,
   !       it will be adjusted to cover all the participating profiles.
   !
   ! note the following params are going to be decided by the input profiles
   !  real :: hc_low       ! Low end for the convoluted profile
   !  real :: hc_high      ! High end for the convoluted profile
   !  real :: delta_z      ! Delta z for the convoluted profile
   ! --------------------------------------------------------------------------


   ! Input
   type(profile_params), dimension(:), pointer :: all_convolute_input


   ! Output
   real, dimension(:), pointer   :: fp         ! foliage profile
   real, dimension(:), pointer   :: efp        ! effective foliage profile
   real, dimension(:), pointer   :: height_levelC  ! Convolute profile height levels


   ! Local variables
   real, dimension(:), allocatable :: arrthl, arrthh, arrdz
   real :: hc_low, hc_high, delta_z
   integer :: n_profiles, iprofile, n_convolute_levels
   type(profile_params), pointer :: interpolate_out


   !---------------------------------------------------------------------

   print *, 'I am in convolute now....'

   n_profiles = size(all_convolute_input)

   print *, all_convolute_input(1)%fp(10)
   print *, all_convolute_input(2)%fp(10)


   ! Get the height levels for the convoluted profile
   !----------------------------------------------------------------------------------------------
   allocate(arrthl(n_profiles))
   allocate(arrthh(n_profiles))
   allocate(arrdz(n_profiles))

   ! compare with all profiles and adjust the low end and high end

   do iprofile = 1, n_profiles
     arrthl(iprofile) = minval(all_convolute_input(iprofile)%height_levels)
     arrthh(iprofile) = maxval(all_convolute_input(iprofile)%height_levels)
     arrdz(iprofile) = all_convolute_input(iprofile)%height_levels(2) - &
                       all_convolute_input(iprofile)%height_levels(1)
   end do

   print *, arrthl
   print *, arrthh
   print *, arrdz

   hc_low=minval(arrthl)
   hc_high=maxval(arrthh)
   delta_z=minval(arrdz)
   print *,hc_low, hc_high, delta_z

   deallocate(arrthl, arrthh, arrdz)

   print *, 'Getting convolute height levels ...'
   call get_convolute_height_level(hc_low, hc_high, delta_z, height_levelC)
   n_convolute_levels=size(height_levelC)


   ! allocate output variables
   print *, 'Allocating for the convoluted outputs ....'

   allocate(fp(n_convolute_levels))
   allocate(efp(n_convolute_levels))

   ! allocate temporary interpolate_out
   print *, 'Allocating for the interpolation output ....'
   allocate(interpolate_out)
   allocate(interpolate_out%fp(n_convolute_levels))
   allocate(interpolate_out%efp(n_convolute_levels))
   allocate(interpolate_out%height_levels(n_convolute_levels))

   interpolate_out%height_levels = height_levelC


   ! interpolate profiles to the convoluted height level
   !---------------------------------------------------------

   ! interpolate profile1 to the convoluted height levels
   do iprofile = 1, n_profiles


      print *, 'Interpolating for profile No. ', iprofile
      call interpolate_profile (all_convolute_input(iprofile),   &
            interpolate_out)
      print *, interpolate_out%fp(10)

      print *, 'Reassign or convolute the output'

      if (iprofile == 1) then
         fp = interpolate_out%fp
         efp = interpolate_out%efp
      else
         fp = fp + interpolate_out%fp
         efp = efp + interpolate_out%efp
      end if
      print *, fp(10)

   end do

   deallocate(interpolate_out%fp)
   deallocate(interpolate_out%efp)
   deallocate(interpolate_out%height_levels)
   deallocate(interpolate_out)


   return


 end subroutine convolute


 ! ----------------------------------------------------------------------------------

 subroutine interpolate_profile (interpolate_in, &
                                interpolate_out)


   ! Input
   type(profile_params), intent(in) :: interpolate_in             ! data for input profile

   ! Output

   type(profile_params), pointer :: interpolate_out            ! interpolated profile

   ! local variables
   integer :: i, j
   integer :: nlevel_in, nlevel_out
   real :: h, hin_low, hin_high
   integer :: id1, id2
   integer, dimension(1) :: id
   real, dimension(:), pointer :: height_in, height_out      ! height levels for input, output profile
   real, dimension(:), pointer :: diff

   ! ---------------------------------------
   print *, interpolate_in%fp(10)


   height_in => interpolate_in%height_levels
   height_out => interpolate_out%height_levels

   nlevel_in = size(height_in)
   nlevel_out = size (height_out)
   hin_low = height_in(1)
   hin_high = height_in(nlevel_in)
   allocate(diff(nlevel_in))


   ! loop through each height level for the convoluted height levels.
   ! Based upon the height comparison to get the interpolated profile

   do i=1, nlevel_out

     h = height_out(i)

     if ( (h > hin_low) .and. (h < hin_high)) then

       ! find the closest two height level
       do j=1, nlevel_in
         diff(j) = (height_in(j) - h)
       end do
       id = minloc ( abs(diff) ) ! find the location that has the height closes to h
       if (diff(id(1)) > 0) then
         id1 = id(1) - 1
         id2 = id(1)
       else
         id1 = id(1)
         id2 = id(1) + 1
       end if

       ! use the data from the above two points to interpolate for the desired height
       

       interpolate_out%fp(i) = interpolate_in%fp(id1) + (h - height_in(id1)) * (interpolate_in%fp(id2)   &
                               - interpolate_in%fp(id1))   / (height_in(id2) - height_in(id1))
       interpolate_out%efp(i) = interpolate_in%efp(id1) + (h - height_in(id1)) * (interpolate_in%efp(id2)   &
                               - interpolate_in%efp(id1))   / (height_in(id2) - height_in(id1))
       

     else

       ! just fill in the output profile with the ceiling or bottom value

       if (h <= hin_low) then
         interpolate_out%fp(i) = interpolate_in%fp(1)
         interpolate_out%efp(i) = interpolate_in%efp(1)
       end if

       if (h >= hin_high) then
         interpolate_out%fp(i) = interpolate_in%fp(nlevel_in)
         interpolate_out%efp(i) = interpolate_in%efp(nlevel_in)
       end if

     end if

   end do

   nullify(height_in)
   nullify(height_out)

   deallocate(diff)
   return

 end subroutine interpolate_profile


 !-------------------------------------------------------------------------------

 subroutine get_convolute_height_level(h1, h2, delta_z, height_levels)

   ! compute the height levels between h1 and h2
   ! NOTE : here last layer could be less than delta_z.
   !        It may be necessary to return the dz at the last layer to
   !        ensure the accuracy of calculation, otherwise, we can make
   !        sure (h2-h1) can be integer divided by delta_z.

   real, intent(in) :: h1, h2, delta_z
   real, dimension(:), pointer :: height_levels

   ! local
   integer :: n_levels, i, temp1
   real :: temp2, dif

   ! -------------------------------------------

   ! Added the rounding checking for different machines.

   temp1 = ceiling((h2-h1)/delta_z)
   temp2 = (h2-h1)/delta_z

   if ((temp1-temp2) > 0.999) then
     n_levels = ceiling((h2-h1)/delta_z)
   else
     n_levels = ceiling((h2-h1)/delta_z) + 1
   end if


   allocate(height_levels(n_levels))

   do i =1, n_levels-1
     height_levels(i) = h1 + delta_z * (i-1)
   end do

   height_levels(n_levels) = h2

 end subroutine get_convolute_height_level

! ----------------------------------------------------

 function get_analytical_clump(gin)

   type(gort_input), intent(in) :: gin
   real :: get_analytical_clump

   ! local variables:
   real :: tr, a, b, c, d
   ! Add by Ping
   real :: G

   a = (tan(radius(gin%zenith))) ** 2
   b = (gin%vert_radius/gin%horz_radius) ** 2
   c = sqrt((1+a)/(1+a*b))
   G = get_G(gin%zenith, gin%pft) ! Add by Ping
   tr = G * gin%dens_foliage * gin%horz_radius * c

   d = (1 - (2*tr+1) * exp(-2*tr)) / (2*tr*tr)

   get_analytical_clump = 3*(1-d)/(4*tr)

 end function get_analytical_clump

! -------------------------------------------------------------


 subroutine swap(x, y)
   ! swaps x and y

   real, intent(inout) :: x, y
   real :: tmp

   tmp = x
   x = y
   y = tmp

 end subroutine swap

! -------------------------------------------------------------------------------------

 function radius(x)
   ! convert degree to radius

   real :: x
   real :: radius

   radius = x*PI/180.
   return

 end function radius

! ------------------------------------------------------

 function height_function1(z, h, b, hdif)

   real :: z, h, b, hdif
   real :: height_function1

   height_function1 = ((b+z-h)**2) * (2*b + h - z) / (3*hdif*b**2)

 end function height_function1

! ------------------------------------------------------
 function height_function2(z, h, b, hdif)

   real :: z, h, b, hdif
   real :: height_function2

   height_function2 = ((b+h-z)**2) * (2*b + z - h) / (3*hdif*b**2)

 end function height_function2
! -------------------------------------------------------

 function height_function3(hdif, b)

   real :: b, hdif
   real :: height_function3

   height_function3 = (4./3.) * b / hdif

 end function height_function3

! -------------------------------------------------------
 function height_function4(z, h1, h2, b)

   real :: z, h1, h2, b
   real :: height_function4

   height_function4 = 1 - ((h2-h1)**2 + 3*(z-h1)*(z-h2)) / (3*b**2)

 end function height_function4

! --------------------------------------------------------
 function get_G(zenith, ivt)

   real :: zenith
   integer :: ivt
   real :: get_G

    integer, parameter :: numpft      =   16   	! number of clm plant function types
    !1: needleleaf_evergreen_temperate_tree   
    !2: needleleaf_evergreen_boreal_tree      

    !3: needleleaf_deciduous_boreal_tree      
    
    !4: broadleaf_evergreen_tropical_tree     
    !5: broadleaf_evergreen_temperate_tree    
    
    !6: broadleaf_deciduous_tropical_tree     
    !7: broadleaf_deciduous_temperate_tree    
    !8: broadleaf_deciduous_boreal_tree       
    
    !9: broadleaf_evergreen_shrub             
    !10:broadleaf_deciduous_temperate_shrub   
    !11:broadleaf_deciduous_boreal_shrub      
    
    !12:c3_arctic_grass                       
    !13:c3_non-arctic_grass                   
    !14:c4_grass                              
    
    !15:crop 1: e.g.corn                                  
    !16:crop 2: e.g.wheat                                 


    integer :: i 
    real :: xl(numpft)         ! ecophys const - leaf/stem orientation index
    data (xl(i),i=1,numpft) /0.01, 0.01, 0.01, 0.10, 0.10, 0.01, 0.25, &
               0.25, 0.01, 0.25, 0.25, -0.30, -0.30, -0.30, -0.30, -0.30/
    real :: phi1,phi2,chil   
                                    
       chil = min( max(xl(ivt), -0.4), 0.6 )
       if (abs(chil) <= 0.01) chil = 0.01
       phi1 = 0.5 - 0.633*chil - 0.330*chil*chil
       phi2 = 0.877 * (1.-2.*phi1)
       get_G = phi1 + phi2*cos(zenith*PI/180)

 end function get_G


END MODULE gort_lidar


