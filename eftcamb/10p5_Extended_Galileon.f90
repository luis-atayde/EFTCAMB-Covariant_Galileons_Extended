!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 10p5_Extended_Galileon.f90
!! This file contains the definition of the Extended Galileon model.
!! Please refer to the numerical notes for details.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Extended Galileon model.
!! Please refer to the numerical notes for details.

!> @author Simone Peirone, Bin Hu, Marco Raveri

module EFTCAMB_full_Extended_Galileon

    use precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_abstract_model_full

    implicit none

    private

    public EFTCAMB_Extended_Galileon

    !----------------------------------------------------------------------------------------
    !> This is the type that contains the definition of the Extended Galileon model.
    type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Extended_Galileon

        ! the model parameters:
        real(dl)  :: ExtendedGalileon_B      !< Extended Galileon model parameter \f$B\f$
        real(dl)  :: ExtendedGalileon_q      !< Extended Galileon model parameter \f$q\f$
        real(dl)  :: csi                   !< Extended Galileon background parameter \f$\xi\f$ deriving from the tracker solution
	real(dl)  :: S
    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBExtendedGalileonReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBExtendedGalileonAllocateModelSelection      !< subroutine that allocates the model selection. For Horava this is a dummy procedure.
        procedure :: init_model_parameters           => EFTCAMBExtendedGalileonInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBExtendedGalileonInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBExtendedGalileonInitBackground              !< subroutine that initializes the background of Extended Galileon.

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBExtendedGalileonComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBExtendedGalileonFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBExtendedGalileonParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBExtendedGalileonParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBExtendedGalileonParameterValues            !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBExtendedGalileonBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBExtendedGalileonSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_adotoa                    => EFTCAMBExtendedGalileonComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBExtendedGalileonComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBExtendedGalileonAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_Extended_Galileon

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBExtendedGalileonReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Extended_Galileon)       :: self   !< the base class
        type(TIniFile)                      :: Ini    !< Input ini file

    end subroutine EFTCAMBExtendedGalileonReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBExtendedGalileonAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_Extended_Galileon) :: self !< the base class

    end subroutine EFTCAMBExtendedGalileonAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    !! Nothing needs to be done but procedure present because it is deferred.
    subroutine EFTCAMBExtendedGalileonInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Extended_Galileon)                          :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

    end subroutine EFTCAMBExtendedGalileonInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file. Nothing needs to be done
    !! but procedure present because it is deferred.
    subroutine EFTCAMBExtendedGalileonInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Extended_Galileon)  :: self   !< the base class
        type(TIniFile)                 :: Ini    !< Input ini file

    end subroutine EFTCAMBExtendedGalileonInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of Extended Galileon.
    subroutine EFTCAMBExtendedGalileonInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                      , intent(out)   :: success        !< wether the background initialization succeded or not

        real(dl) :: Omega_phi0

        Omega_phi0 = params_cache%omegav

        ! Extended -> just c_3
!        self%csi              = sqrt( 6._dl*Omega_phi0 )		!!!!!!!!!
        self%ExtendedGalileon_B = self%S*self%ExtendedGalileon_q			!< B=S*q!
        self%ExtendedGalileon_q = self%ExtendedGalileon_B/self%S			!!!!!!!!!
	self%S			= self%ExtendedGalileon_B/self%ExtendedGalileon_q
        call self%feedback()

        success=.true.

    end subroutine EFTCAMBExtendedGalileonInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBExtendedGalileonComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Extended_Galileon)  :: self   !< the base class

        self%parameter_number = 2								!luis-change::0->2

    end subroutine EFTCAMBExtendedGalileonComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBExtendedGalileonFeedback( self, print_params )

        implicit none

        class(EFTCAMB_Extended_Galileon)  :: self         !< the base class
        logical, optional              :: print_params !< optional flag that decised whether to print numerical values
                                                       !! of the parameters.

        ! print general model informations:
!        if (self%ExtendedGalileon_c2 == -1._dl) then 
!
!          write(*,*)
!          write(*,'(a,a)')       '   Model              =  ', self%name
!          write(*,'(a,I3)')      '   Number of params   ='  , self%parameter_number
!          write(*,'(a,F12.6)')   '                xi    ='  , self%csi
!          write(*,'(a,F12.6)')   '                B    ='  , self%ExtendedGalileon_B
!          write(*,'(a,F12.6)')   '                q    ='  , self%ExtendedGalileon_q
!          write(*,'(a,F12.6)')   '                S    ='  , self%S
!        end if

    end subroutine EFTCAMBExtendedGalileonFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Extended_Galileon) :: self   !< the base class
        integer     , intent(in)      :: i      !< the index of the parameter
        character(*), intent(out)     :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==0 ) then
            name = 'no_name'
            return
        end if

    end subroutine EFTCAMBExtendedGalileonParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Extended_Galileon) :: self       !< the base class
        integer     , intent(in)      :: i          !< The index of the parameter
        character(*), intent(out)     :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_names_latex.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==0 ) then
            latexname = 'noname'
            return
        end if

    end subroutine EFTCAMBExtendedGalileonParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBExtendedGalileonParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Extended_Galileon) :: self   !< the base class
        integer , intent(in)          :: i      !< The index of the parameter
        real(dl), intent(out)         :: value  !< the output value of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'Illegal index for parameter_value.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the appropriate name:
        if ( i==0 ) then
            value = 0._dl
            return
        end if

    end subroutine EFTCAMBExtendedGalileonParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBExtendedGalileonBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache)

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
       real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.


        real(dl)    :: a2, Omega_phi0

        Omega_phi0 = eft_par_cache%omegav
        a2 = a*a

        if(a==0._dl) then
            return
        else if (eft_cache%adotoa==0._dl) then
            call self%compute_adotoa( a, eft_par_cache, eft_cache )
            call self%compute_H_derivs( a, eft_par_cache, eft_cache )
        end if
!
        !! compute the psi field and its derivatives
        !psi           = self%csi*(eft_par_cache%h0_Mpc)**2*a2/eft_cache%adotoa**2
        !psiprime      = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2*a/eft_cache%adotoa**4*( eft_cache%adotoa**2 -eft_cache%Hdot )
        !psiprimeprime = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2/eft_cache%adotoa**4*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot &
        !    & +4._dl*( eft_cache%Hdot/eft_cache%adotoa )**2 -eft_cache%Hdotdot/eft_cache%adotoa )

        ! scalar field conversion
        !phip1 = psi/a
        !phip2 = ( psiprime*a -psi )/a2
        !phip3 = ( psiprimeprime*a2 -2._dl*psiprime*a +2._dl*psi )/a**3

        !c3= -1.0_dl*self%ExtendedGalileon_c3/(eft_par_cache%h0_Mpc)**2

        ! compute the background EFT functions:
        eft_cache%EFTOmegaV    = 0._dl
        eft_cache%EFTOmegaP    = 0._dl
        eft_cache%EFTOmegaPP   = 0._dl
        eft_cache%EFTOmegaPPP  = 0._dl
        eft_cache%EFTc         = 6*a2*self%ExtendedGalileon_B*(Omega_phi0)**2*eft_par_cache%h0_Mpc**(2*&
(self%ExtendedGalileon_B+self%S+1))*eft_cache%adotoa**((-2*self%ExtendedGalileon_B&
*(self%ExtendedGalileon_q+1))/self%ExtendedGalileon_q)*&
((eft_par_cache%h0_Mpc**(2*self%ExtendedGalileon_q+1)*6**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*&
eft_cache%adotoa**(-2*self%ExtendedGalileon_q-3)*a**(4*self%ExtendedGalileon_B+2*self%ExtendedGalileon_B/self%ExtendedGalileon_q+4*&
self%ExtendedGalileon_q+1)*Omega_phi0**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*(a*(2*self%ExtendedGalileon_q+1)*eft_cache%adotoa*&
eft_cache%Hdot+(self%ExtendedGalileon_q-1)*eft_cache%Hdot**2-self%ExtendedGalileon_q*eft_cache%Hdot))/self%ExtendedGalileon_q-3*&
a**(2*self%ExtendedGalileon_B*(1/self%ExtendedGalileon_q+2)))




        eft_cache%EFTLambda    =6*a2*(Omega_phi0)**2*eft_par_cache%h0_Mpc**(2*(self%ExtendedGalileon_B+self%S+1))*&
eft_cache%adotoa**((-2*self%ExtendedGalileon_B*(self%ExtendedGalileon_q+1))/self%ExtendedGalileon_q)*&
(-(self%ExtendedGalileon_B*eft_par_cache%h0_Mpc**(2*self%ExtendedGalileon_q+1)*2**((self%ExtendedGalileon_B+self%ExtendedGalileon_q)&
/self%ExtendedGalileon_B)*3**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*eft_cache%adotoa**(-2*self%ExtendedGalileon_q-3)*&
a**(4*self%ExtendedGalileon_B+2*self%S+4*self%ExtendedGalileon_q+1)*Omega_phi0**(1/self%S)*((2*self%ExtendedGalileon_q+1)&
*eft_cache%adotoa**2-(self%ExtendedGalileon_q+1)*eft_cache%Hdot))/self%ExtendedGalileon_q-3*a**(2*self%ExtendedGalileon_B*(1/self%ExtendedGalileon_q+2)))

       eft_cache%EFTcdot      =(6*self%ExtendedGalileon_B)&
/self%ExtendedGalileon_q**2*&
Omega_phi0**2*&
a**(2*self%ExtendedGalileon_B*(1/self%ExtendedGalileon_q+2)+2)*&
eft_par_cache%h0_Mpc**(2*(self%ExtendedGalileon_B+self%ExtendedGalileon_B/self%ExtendedGalileon_q+1))*&
eft_cache%adotoa**(-2*(self%ExtendedGalileon_B+self%S+self%ExtendedGalileon_q+2))*&
(a**(4*self%ExtendedGalileon_q+1)*&
eft_par_cache%h0_Mpc**(2*self%ExtendedGalileon_q+1)*&
6**(1/self%S)&
*Omega_phi0**(1/self%S)&
*(self%ExtendedGalileon_q*(self%ExtendedGalileon_q+1)*eft_cache%adotoa*eft_cache%Hdotdot+2*eft_cache%adotoa**2*(2*self%ExtendedGalileon_B+self%ExtendedGalileon_q*(3*self%ExtendedGalileon_B+self%ExtendedGalileon_q*(self%ExtendedGalileon_B+self%ExtendedGalileon_q+3)+1))*&
eft_cache%Hdot+(self%ExtendedGalileon_q-1)*eft_cache%adotoa**4*(2*self%ExtendedGalileon_B+4*self%ExtendedGalileon_q**2+4*self%ExtendedGalileon_B*self%ExtendedGalileon_q+self%ExtendedGalileon_q)-(self%ExtendedGalileon_q+1)*(2*self%ExtendedGalileon_B*(self%ExtendedGalileon_q+1)+self%ExtendedGalileon_q*(2*self%ExtendedGalileon_q+3))*eft_cache%Hdot**2)-6*self%ExtendedGalileon_B*&
self%ExtendedGalileon_q*&
eft_cache%adotoa**(2*self%ExtendedGalileon_q+3)*&
((2*self%ExtendedGalileon_q+1)*eft_cache%adotoa**2-(self%ExtendedGalileon_q+1)*eft_cache%Hdot)&
)
            
 
   
         
       eft_cache%EFTLambdadot = ((6**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*a**(1 + 4*self%ExtendedGalileon_q +&
 4*self%ExtendedGalileon_B + (2*self%ExtendedGalileon_B)/self%ExtendedGalileon_q)*(3 + 2*self%ExtendedGalileon_q)*self%ExtendedGalileon_B*eft_par_cache%h0_Mpc**(1 + 2*self%ExtendedGalileon_q)*Omega_phi0**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*eft_cache%Hdot*((-1 - 2*self%ExtendedGalileon_q)*eft_cache%adotoa**2 + (1 + self%ExtendedGalileon_q)*eft_cache%Hdot))/(self%ExtendedGalileon_q*eft_cache%adotoa**(2*self%ExtendedGalileon_q))) - eft_cache%adotoa**5*(-3*a**(2*(2 + 1/self%ExtendedGalileon_q)*self%ExtendedGalileon_B) - (2**((self%ExtendedGalileon_q + self%ExtendedGalileon_B)/self%ExtendedGalileon_B)*3**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*a**(1 + 4*self%ExtendedGalileon_q + 4*self%ExtendedGalileon_B + (2*self%ExtendedGalileon_B)/self%ExtendedGalileon_q)*self%ExtendedGalileon_B*eft_cache%adotoa**(-3 - 2*self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc**(1 + 2*self%ExtendedGalileon_q)*Omega_phi0**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*((1 + 2*self%ExtendedGalileon_q)*eft_cache%adotoa**2 - (1 + self%ExtendedGalileon_q)*eft_cache%Hdot))/self%ExtendedGalileon_q) - ((1 + self%ExtendedGalileon_q)*self%ExtendedGalileon_B*eft_cache%adotoa**3*eft_cache%Hdot*(-3*a**(2*(2 + 1/self%ExtendedGalileon_q)*self%ExtendedGalileon_B) - (2**((self%ExtendedGalileon_q + self%ExtendedGalileon_B)/self%ExtendedGalileon_B)*3**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*a**(1 + 4*self%ExtendedGalileon_q + 4*self%ExtendedGalileon_B + (2*self%ExtendedGalileon_B)/self%ExtendedGalileon_q)*self%ExtendedGalileon_B*eft_cache%adotoa**(-3 - 2*self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc**(1 + 2*self%ExtendedGalileon_q)*Omega_phi0**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*((1 + 2*self%ExtendedGalileon_q)*eft_cache%adotoa**2 - (1 + self%ExtendedGalileon_q)*eft_cache%Hdot))/self%ExtendedGalileon_q))/self%ExtendedGalileon_q - (6**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*a**(1 + 4*self%ExtendedGalileon_q + 4*self%ExtendedGalileon_B + (2*self%ExtendedGalileon_B)/self%ExtendedGalileon_q)*self%ExtendedGalileon_B*eft_cache%adotoa**(1 - 2*self%ExtendedGalileon_q)*eft_par_cache%h0_Mpc**(1 + 2*self%ExtendedGalileon_q)*Omega_phi0**(self%ExtendedGalileon_q/self%ExtendedGalileon_B)*(2*(1 + 2*self%ExtendedGalileon_q)*eft_cache%adotoa*eft_cache%Hdot - (1 + self%ExtendedGalileon_q)*eft_cache%Hdotdot))/self%ExtendedGalileon_q

    end subroutine EFTCAMBExtendedGalileonBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the second order EFT functions at a given time.
    subroutine EFTCAMBExtendedGalileonSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: a2,Omega_phi0

        Omega_phi0 = eft_par_cache%omegav
        a2 = a*a

        if(a*eft_cache%adotoa==0._dl) return

!        ! compute the psi field and its derivatives
!        psi           = self%csi*(eft_par_cache%h0_Mpc)**2*a2/eft_cache%adotoa**2
!        psiprime      = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2*a/eft_cache%adotoa**4*( eft_cache%adotoa**2 -eft_cache%Hdot )
!        psiprimeprime = 2._dl*self%csi*(eft_par_cache%h0_Mpc)**2/eft_cache%adotoa**4*( eft_cache%adotoa**2 -3._dl*eft_cache%Hdot &
!            & +4._dl*( eft_cache%Hdot/eft_cache%adotoa )**2 -eft_cache%Hdotdot/eft_cache%adotoa )
!
!        ! scalar field conversion
!        phip1 = psi/a
!        phip2 = ( psiprime*a -psi )/a2
!        phip3 = ( psiprimeprime*a2 -2._dl*psiprime*a +2._dl*psi )/a**3

!        c3    = -1.0_dl*self%ExtendedGalileon_c3/(eft_par_cache%h0_Mpc)**2

        ! compute the second order EFT functions:
        eft_cache%EFTGamma1V  = self%S/12*a**(2*self%ExtendedGalileon_B*(1/self%ExtendedGalileon_q+2)-1)*eft_par_cache%h0_Mpc**(2*self%ExtendedGalileon_B*(self%ExtendedGalileon_q+1)/self%ExtendedGalileon_q)*eft_cache%adotoa**(-2*self%ExtendedGalileon_B*(self%ExtendedGalileon_q+1)/self%ExtendedGalileon_q-2*self%ExtendedGalileon_q-3)*(6**(1/self%S)*Omega_phi0**(1/self%S))**(2*self%S)*(a**(4*self%ExtendedGalileon_q+1)*eft_par_cache%h0_Mpc**(2*self%ExtendedGalileon_q+1)*6**(1/self%S)*Omega_phi0**(1/self%S)*(eft_cache%adotoa**2*(2*self%ExtendedGalileon_q*(3*self%ExtendedGalileon_B+3*self%ExtendedGalileon_q-2)+1)-(self%ExtendedGalileon_q+1)*eft_cache%Hdot)-6*(self%ExtendedGalileon_B-1)*self%ExtendedGalileon_q*eft_cache%adotoa**(2*self%ExtendedGalileon_q+3))
             
        eft_cache%EFTGamma1P  = Omega_phi0**2*self%S*3/self%ExtendedGalileon_q*a**(2*self%ExtendedGalileon_B*(1/self%ExtendedGalileon_q+2)-1)*eft_par_cache%h0_Mpc**(2*self%ExtendedGalileon_B*(self%ExtendedGalileon_q+1)/self%ExtendedGalileon_q)*eft_cache%adotoa**(-2*&
(self%ExtendedGalileon_B+self%ExtendedGalileon_q)*(self%ExtendedGalileon_q+1)/self%ExtendedGalileon_q-3)*(a**(4*self%ExtendedGalileon_q+1)*eft_par_cache%h0_Mpc**(2*self%ExtendedGalileon_q+1)*6**(1/self%S)*Omega_phi0**(1/self%S)*&
(-self%ExtendedGalileon_q*(self%ExtendedGalileon_q+1)*eft_cache%adotoa*eft_cache%Hdotdot-eft_cache%adotoa**2*(4*self%ExtendedGalileon_B+self%ExtendedGalileon_q*(12*self%ExtendedGalileon_B**2*(self%ExtendedGalileon_q+1)+2*self%ExtendedGalileon_B*&
self%ExtendedGalileon_q*(12*self%ExtendedGalileon_q+7)+self%ExtendedGalileon_q*(2*self%ExtendedGalileon_q*(6*self%ExtendedGalileon_q+1)+3)+2))*eft_cache%Hdot+eft_cache%adotoa**4*(2*self%ExtendedGalileon_B+4*self%ExtendedGalileon_q**2+4*self%ExtendedGalileon_q*&
self%ExtendedGalileon_B+self%ExtendedGalileon_q)*(2*self%ExtendedGalileon_q*(3*self%ExtendedGalileon_B+3*self%ExtendedGalileon_q-2)+1)+(self%ExtendedGalileon_q+1)*(2*self%ExtendedGalileon_B*(self%ExtendedGalileon_q+1)+self%ExtendedGalileon_q*(2*self%ExtendedGalileon_q+3))*eft_cache%Hdot**2)-12*(self%ExtendedGalileon_B-1)*self%ExtendedGalileon_B*self%ExtendedGalileon_q*eft_cache%adotoa**(2*self%ExtendedGalileon_q+3)*((2*self%ExtendedGalileon_q+1)*eft_cache%adotoa**2-(self%ExtendedGalileon_q+1)*eft_cache%Hdot))

        eft_cache%EFTGamma2V  = self%ExtendedGalileon_B*(-2**(1/self%S+2))*3**((self%ExtendedGalileon_B+self%ExtendedGalileon_q)/self%ExtendedGalileon_B)*(a*eft_par_cache%h0_Mpc/eft_cache%adotoa)**((2*(self%ExtendedGalileon_q+1)*(self%ExtendedGalileon_B+self%ExtendedGalileon_q))/self%ExtendedGalileon_q)*Omega_phi0**(self%S+2)
        eft_cache%EFTGamma2P  = 2/3*self%S*(self%ExtendedGalileon_B+self%ExtendedGalileon_q)*a**(4*self%ExtendedGalileon_B+2*self%S+4*self%ExtendedGalileon_q+1)*eft_par_cache%h0_Mpc**((2*(self%ExtendedGalileon_q+1)*(self%ExtendedGalileon_B+self%ExtendedGalileon_q))/self%ExtendedGalileon_q)*eft_cache%adotoa**(-2*(self%ExtendedGalileon_B+self%S+self%ExtendedGalileon_q+2))*(6**(1/self%S)*Omega_phi0**(1/self%S))**(2*self%S+1)*((2*self%ExtendedGalileon_q+1)*eft_cache%adotoa**2-(self%ExtendedGalileon_q+1)*eft_cache%Hdot)
        eft_cache%EFTGamma3V  = 0._dl
        eft_cache%EFTGamma3P  = 0._dl
        eft_cache%EFTGamma4V  = 0._dl
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl
        eft_cache%EFTGamma5V  = 0._dl
        eft_cache%EFTGamma5P  = 0._dl
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBExtendedGalileonSecondOrderEFTFunctions
 	

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H.
    subroutine EFTCAMBExtendedGalileonComputeAdotoa( self, a, eft_par_cache, eft_cache )
	
        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
	
	
       
         	!real(dl), intent(inout)::			first_value
        real(dl)    :: temp, a2, Omega_tot
        integer     :: nu_i , counter 
	real(dl)::limit1, limit2, flimit1, flimit2, dmean, soluction, fsoluction, bolean !soluction=H/H0
!

        real(dl) :: Omega_phi0
	real(dl) :: ATemp1, ATemp2, BTemp1, BTemp2, HorizAsyntB

        Omega_phi0 = eft_par_cache%omegav
        a2 = a*a

        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)



	
		limit1=0
		if (limit1.lt.0) limit1=0
		limit2=10**(9)
		flimit1=Omega_phi0+Omega_tot*(limit1/a)**(self%S)-(limit1/a)**(2+self%S)
		flimit2=Omega_phi0+Omega_tot*(limit2/a)**(self%S)-(limit2/a)**(2+self%S)
		dmean=(limit2-limit1)/2
		soluction=limit2-dmean
		fsoluction=1
		counter=0
		do while(sqrt(fsoluction**2).gt.10**(-1).and.counter.lt.50**1)
			fsoluction=Omega_phi0+Omega_tot*(soluction/a)**(self%S)-(soluction/a)**(2+self%S)
			bolean=fsoluction*flimit1
			if (bolean.gt.0.) then
				limit1=soluction
				flimit1=fsoluction
			endif
			if (bolean.le.0.) then
				limit2=soluction
				flimit2=fsoluction
			endif
			dmean=(limit2-limit1)/2
			soluction=limit1+dmean
			counter=counter+1
			enddo
			
		!first_value=soluction
		temp= soluction*eft_par_cache%h0_Mpc
		
		

        eft_cache%adotoa = temp

    end subroutine EFTCAMBExtendedGalileonComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBExtendedGalileonComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl)    :: temp, a2, Omega_tot, Omega_tot_prime, Omega_tot_primeprime, Omega_phi0
        integer     :: nu_i

        a2 = a*a

        Omega_tot = ( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-3) + ( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-4) +eft_cache%grhonu_tot/(3._dl*eft_par_cache%h0_Mpc**2*a2)
        Omega_tot_prime = -3._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-4) -4._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-5) &
                          & -(eft_cache%grhonu_tot+eft_cache%gpinu_tot)/(eft_par_cache%h0_Mpc**2*a2*a)
        Omega_tot_primeprime = 12._dl*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-5) +20._dl*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-6)&
                          & +(4._dl*(eft_cache%grhonu_tot+eft_cache%gpinu_tot)-eft_cache%gpinudot_tot/eft_cache%adotoa )/(eft_par_cache%h0_Mpc**2*a2**2)
!	Omega_tot_primeprimeprime=-60*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-6)-120*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-7)+1/(eft_par_cache%h0_Mpc**2*a2**2*a)*(-20*(eft_cache%grhonu_tot+eft_cache%gpinu_tot)+eft_cache%gpinudot_tot*(6/eft_cache%adotoa+eft_cache%Hdot/(eft_cache%adotoa**3))-1/(eft_cache%adotoa**2)*eft_cache%gpinudotdot_tot)
!	omega_tot_primeprimeprimeprime=360*( eft_par_cache%omegac +eft_par_cache%omegab )*a**(-7)-840*( eft_par_cache%omegag +eft_par_cache%omegar)*a**(-8)+1/(eft_par_cache%h0_Mpc**2*a2**3)*(120*(eft_cache%grhonu_tot+eft_cache%gpinu_tot)+eft_cache%gpinudot_tot*(-38/eft_cache%adotoa+-9*eft_cache%Hdot/(eft_cache%adotoa**3)+eft_cache%Hdotdot/(eft_cache%adotoa**4)-3*eft_cache%Hdot**2/(eft_cache%adotoa**5))-1/(eft_cache%adotoa**3)*eft_cache%gpinudotdotdot_tot+eft_cache%gpinudotdot_tot*(9/(eft_cache%adotoa**2)+3*eft_cache%Hdot/eft_cache%adotoa**4)
        Omega_phi0 = eft_par_cache%omegav

        eft_cache%Hdot = -(eft_cache%adotoa**2*(a**2*eft_par_cache%h0_Mpc**2*(a*Omega_tot_prime-self%ExtendedGalileon_q*Omega_tot)+&
(self%ExtendedGalileon_q+2)*eft_cache%adotoa**2)/(a**2*eft_par_cache%h0_Mpc**2*self%ExtendedGalileon_q*Omega_tot-&
(self%ExtendedGalileon_q+2)*eft_cache%adotoa**2))
        eft_cache%Hdotdot = eft_cache%adotoa**3/((self%ExtendedGalileon_q+2)*eft_cache%adotoa**2-a**2*eft_par_cache%h0_Mpc**2*self%ExtendedGalileon_q*Omega_tot)**3*(a**2*eft_par_cache%h0_Mpc**2*(self%ExtendedGalileon_q+2)**2*eft_cache%adotoa**4*(a*(5*Omega_tot_prime+a*Omega_tot_primeprime)-6*self%ExtendedGalileon_q*Omega_tot)+a**6*eft_par_cache%h0_Mpc**6*self%ExtendedGalileon_q*Omega_tot*&
(-a**2*(self%ExtendedGalileon_q+2)*Omega_tot_prime**2-2*self%ExtendedGalileon_q**2*Omega_tot**2+a*self%ExtendedGalileon_q*Omega_tot*(5*Omega_tot_prime+a*Omega_tot_primeprime))+a**4*eft_par_cache%h0_Mpc**4*self%ExtendedGalileon_q*(self%ExtendedGalileon_q+2)*eft_cache%adotoa**2*(a**2*Omega_tot_prime**2+6*self%ExtendedGalileon_q*Omega_tot**2-2*a*Omega_tot*(5*Omega_tot_prime+a*Omega_tot_primeprime))+2*(self%ExtendedGalileon_q+2)**3*eft_cache%adotoa**6)
!	eft_cache%Hdotdotdot=  eft_cache%adotoa**4*(5*(2 + self%ExtendedGalileon_q)**5*eft_cache%adotoa**10 + ((2 + self%ExtendedGalileon_q)**5*eft_cache%adotoa**10*((2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*eft_par_cache%h0_Mpc**2*(-(self%ExtendedGalileon_q*Omega_tot) + a*Omega_tot_prime)))/&
!         ((2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2 - a**2*self%ExtendedGalileon_q*eft_par_cache%h0_Mpc**2*Omega_tot) + &
!        (a**10*self%ExtendedGalileon_q**3*eft_par_cache%h0_Mpc**10*Omega_tot**3*(self%ExtendedGalileon_q*Omega_tot - a*Omega_tot_prime)**2*&
!         ((2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*eft_par_cache%h0_Mpc**2*(-(self%ExtendedGalileon_q*Omega_tot) + a*Omega_tot_prime)))/&
!      ((-2 - self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*self%ExtendedGalileon_q*eft_par_cache%h0_Mpc**2*Omega_tot) + &
!    (a**2*(2 + self%ExtendedGalileon_q)**4*eft_cache%adotoa**8*eft_par_cache%h0_Mpc**2*(5*self%ExtendedGalileon_q*Omega_tot + 2*a*Omega_tot_prime)*&
!      ((2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*eft_par_cache%h0_Mpc**2*(-(self%ExtendedGalileon_q*Omega_tot) + a*Omega_tot_prime)))/&
!          ((-2 - self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*self%ExtendedGalileon_q*eft_par_cache%h0_Mpc**2*Omega_tot) - &
!        (a**6*self%ExtendedGalileon_q*(2 + self%ExtendedGalileon_q)**2*eft_cache%adotoa**4*eft_par_cache%h0_Mpc**6*Omega_tot*(-10*self%ExtendedGalileon_q**2*Omega_tot**2 + a**2*Omega_tot_prime**2)*&
!          ((2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*eft_par_cache%h0_Mpc**2*(-(self%ExtendedGalileon_q*Omega_tot) + a*Omega_tot_prime)))/&
!       ((-2 - self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*self%ExtendedGalileon_q*eft_par_cache%h0_Mpc**2*Omega_tot) + &
!     (a**8*self%ExtendedGalileon_q**2*(2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2*eft_par_cache%h0_Mpc**8*Omega_tot**2*&
!       (-5*self%ExtendedGalileon_q**2*Omega_tot**2 + 4*a*self%ExtendedGalileon_q*Omega_tot*Omega_tot_prime + a**2*Omega_tot_prime**2)*&
!      ((2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*eft_par_cache%h0_Mpc**2*(-(self%ExtendedGalileon_q*Omega_tot) + a*Omega_tot_prime)))/&
!   ((-2 - self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*self%ExtendedGalileon_q*eft_par_cache%h0_Mpc**2*Omega_tot) - &
! (a**4*(2 + self%ExtendedGalileon_q)**3*eft_cache%adotoa**6*eft_par_cache%h0_Mpc**4*(10*self%ExtendedGalileon_q**2*Omega_tot**2 + 4*a*self%ExtendedGalileon_q*Omega_tot*Omega_tot_prime + &
!     a**2*Omega_tot_prime**2)*((2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*eft_par_cache%h0_Mpc**2*(-(self%ExtendedGalileon_q*Omega_tot) + a*Omega_tot_prime)))/&
!          ((-2 - self%ExtendedGalileon_q)*eft_cache%adotoa**2 + a**2*self%ExtendedGalileon_q*eft_par_cache%h0_Mpc**2*Omega_tot) + &
!       a**10*self%ExtendedGalileon_q**2*eft_par_cache%h0_Mpc**10*Omega_tot**2*(-5*self%ExtendedGalileon_q**3*Omega_tot**3 + a**3*(5 + 7*self%ExtendedGalileon_q + 2*self%ExtendedGalileon_q**2)*Omega_tot_prime**3 - &
!        a**2*self%ExtendedGalileon_q*Omega_tot*Omega_tot_prime*(2*(11 + 5*self%ExtendedGalileon_q)*Omega_tot_prime + a*(7 + 3*self%ExtendedGalileon_q)*Omega_tot_primeprime) + &
!       a*self%ExtendedGalileon_q**2*Omega_tot**2*(23*Omega_tot_prime + a*(10*Omega_tot_primeprime + a*Omega_tot_primeprimeprime))) + &
!     a**8*self%ExtendedGalileon_q*(2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2*eft_par_cache%h0_Mpc**8*Omega_tot*(25*self%ExtendedGalileon_q**3*Omega_tot**3 - 4*a**3*(-1 + self%ExtendedGalileon_q + self%ExtendedGalileon_q**2)*Omega_tot_prime**3 + &
!      3*a**2*self%ExtendedGalileon_q*Omega_tot*Omega_tot_prime*(2*(9 + 5*self%ExtendedGalileon_q)*Omega_tot_prime + a*(5 + 3*self%ExtendedGalileon_q)*Omega_tot_primeprime) - &
!     4*a*self%ExtendedGalileon_q**2*Omega_tot**2*(24*Omega_tot_prime + a*(10*Omega_tot_primeprime + a*Omega_tot_primeprimeprime))) - &
! a**6*(2 + self%ExtendedGalileon_q)**2*eft_cache%adotoa**4*eft_par_cache%h0_Mpc**6*(50*self%ExtendedGalileon_q**3*Omega_tot**3 + a**3*(1 + 3*self%ExtendedGalileon_q - 2*self%ExtendedGalileon_q**2)*Omega_tot_prime**3 + &
!      a**2*self%ExtendedGalileon_q*Omega_tot*Omega_tot_prime*((38 + 30*self%ExtendedGalileon_q)*Omega_tot_prime + 9*a*(1 + self%ExtendedGalileon_q)*Omega_tot_primeprime) - &
!          6*a*self%ExtendedGalileon_q**2*Omega_tot**2*(25*Omega_tot_prime + a*(10*Omega_tot_primeprime + a*Omega_tot_primeprimeprime))) +1 &
!        a**4*(2 + self%ExtendedGalileon_q)**3*eft_cache%adotoa**6*eft_par_cache%h0_Mpc**4*(50*self%ExtendedGalileon_q**2*Omega_tot**2 + &
!        a**2*Omega_tot_prime*(2*(3 + 5*self%ExtendedGalileon_q)*Omega_tot_prime + a*(1 + 3*self%ExtendedGalileon_q)*Omega_tot_primeprime) - &
!       4*a*self%ExtendedGalileon_q*Omega_tot*(26*Omega_tot_prime + a*(10*Omega_tot_primeprime + a*Omega_tot_primeprimeprime))) + &
!   a**2*(2 + self%ExtendedGalileon_q)**4*eft_cache%adotoa**8*eft_par_cache%h0_Mpc**2*(-25*self%ExtendedGalileon_q*Omega_tot + &
!     a*(27*Omega_tot_prime + a*(10*Omega_tot_primeprime + a*Omega_tot_primeprimeprime))))

!	eft_cache%Hdotdotdotdot=(eft_cache%adotoa**5*(24*(2 + self%ExtendedGalileon_q)**7*eft_cache%adotoa**14 + a**2*(2 + self%ExtendedGalileon_q)**6*eft_cache%adotoa**12*eft_par_cache%h0_Mpc**2*&
!            (-168*self%ExtendedGalileon_q*Omega_tot + a*(154*Omega_tot_prime + &
!                 a*(86*Omega_tot_primeprime + 17*a*Omega_tot_primeprimeprime + a**2*Omega_tot_primeprimeprimeprime))) + &
!           a**4*(2 + self%ExtendedGalileon_q)**5*eft_cache%adotoa**10*eft_par_cache%h0_Mpc**4*(504*self%ExtendedGalileon_q**2*Omega_tot**2 + &
!              a**2*((82 + 86*self%ExtendedGalileon_q)*Omega_tot_prime**2 + a**2*(1 + 3*self%ExtendedGalileon_q)*Omega_tot_primeprime**2 + &
!                 a*Omega_tot_prime*((37 + 51*self%ExtendedGalileon_q)*Omega_tot_primeprime + a*(3 + 4*self%ExtendedGalileon_q)*Omega_tot_primeprimeprime)) - &
!              6*a*self%ExtendedGalileon_q*Omega_tot*(154*Omega_tot_prime + &
!                 a*(86*Omega_tot_primeprime + 17*a*Omega_tot_primeprimeprime + a**2*Omega_tot_primeprimeprimeprime))) - &
!           a**6*self%ExtendedGalileon_q*(2 + self%ExtendedGalileon_q)**4*eft_cache%adotoa**8*eft_par_cache%h0_Mpc**6*(840*self%ExtendedGalileon_q**2*Omega_tot**3 + &
!              a**3*Omega_tot_prime**2*((31 - 34*self%ExtendedGalileon_q)*Omega_tot_prime + a*(13 - 12*self%ExtendedGalileon_q)*Omega_tot_primeprime) + &
!              a**2*Omega_tot*((582 + 430*self%ExtendedGalileon_q)*Omega_tot_prime**2 + a**2*(11 + 15*self%ExtendedGalileon_q)*Omega_tot_primeprime**2 + &
!                 a*Omega_tot_prime*((287 + 255*self%ExtendedGalileon_q)*Omega_tot_primeprime + a*(23 + 20*self%ExtendedGalileon_q)*Omega_tot_primeprimeprime)) - &
!              15*a*self%ExtendedGalileon_q*Omega_tot**2*(154*Omega_tot_prime + &
!                 a*(86*Omega_tot_primeprime + 17*a*Omega_tot_primeprimeprime + a**2*Omega_tot_primeprimeprimeprime))) + &
!           a**8*self%ExtendedGalileon_q*(2 + self%ExtendedGalileon_q)**3*eft_cache%adotoa**6*eft_par_cache%h0_Mpc**8*(840*self%ExtendedGalileon_q**3*Omega_tot**4 + a**4*(12 - 23*self%ExtendedGalileon_q + 6*self%ExtendedGalileon_q**2)*Omega_tot_prime**4 - &
!              2*a**3*Omega_tot*Omega_tot_prime**2*&
!               ((-31 + 23*self%ExtendedGalileon_q + 68*self%ExtendedGalileon_q**2)*Omega_tot_prime + a*(-13 + 4*self%ExtendedGalileon_q + 24*self%ExtendedGalileon_q**2)*Omega_tot_primeprime) + &
!              2*a**2*self%ExtendedGalileon_q*Omega_tot**2*((754 + 430*self%ExtendedGalileon_q)*Omega_tot_prime**2 + a**2*(17 + 15*self%ExtendedGalileon_q)*Omega_tot_primeprime**2 + &
!                 a*Omega_tot_prime*((389 + 255*self%ExtendedGalileon_q)*Omega_tot_primeprime + a*(31 + 20*self%ExtendedGalileon_q)*Omega_tot_primeprimeprime)) - &
!              20*a*self%ExtendedGalileon_q**2*Omega_tot**3*(154*Omega_tot_prime + &
!                 a*(86*Omega_tot_primeprime + 17*a*Omega_tot_primeprimeprime + a**2*Omega_tot_primeprimeprimeprime))) - &
!           a**10*self%ExtendedGalileon_q*(2 + self%ExtendedGalileon_q)**2*eft_cache%adotoa**4*eft_par_cache%h0_Mpc**10*Omega_tot*&
!            (504*self%ExtendedGalileon_q**4*Omega_tot**4 + a**4*(24 - 70*self%ExtendedGalileon_q - 17*self%ExtendedGalileon_q**2 + 18*self%ExtendedGalileon_q**3)*Omega_tot_prime**4 - &
!             6*a**3*self%ExtendedGalileon_q*Omega_tot*Omega_tot_prime**2*&
!             ((3 + 54*self%ExtendedGalileon_q + 34*self%ExtendedGalileon_q**2)*Omega_tot_prime + a*(-1 + 17*self%ExtendedGalileon_q + 12*self%ExtendedGalileon_q**2)*Omega_tot_primeprime) + &
!           2*a**2*self%ExtendedGalileon_q**2*Omega_tot**2*((926 + 430*self%ExtendedGalileon_q)*Omega_tot_prime**2 + a**2*(23 + 15*self%ExtendedGalileon_q)*Omega_tot_primeprime**2 + &
!             a*Omega_tot_prime*((491 + 255*self%ExtendedGalileon_q)*Omega_tot_primeprime + a*(39 + 20*self%ExtendedGalileon_q)*Omega_tot_primeprimeprime)) - &
!         15*a*self%ExtendedGalileon_q**3*Omega_tot**3*(154*Omega_tot_prime + &
!           a*(86*Omega_tot_primeprime + a*(17*Omega_tot_primeprimeprime + a*Omega_tot_primeprimeprimeprime)))) + &
!      a**12*self%ExtendedGalileon_q**2*(2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2*eft_par_cache%h0_Mpc**12*Omega_tot**2*&
!       (168*self%ExtendedGalileon_q**4*Omega_tot**4 + a**4*(-72 - 36*self%ExtendedGalileon_q + 35*self%ExtendedGalileon_q**2 + 18*self%ExtendedGalileon_q**3)*Omega_tot_prime**4 - &
!        2*a**3*self%ExtendedGalileon_q*Omega_tot*Omega_tot_prime**2*&
!        ((111 + 193*self%ExtendedGalileon_q + 68*self%ExtendedGalileon_q**2)*Omega_tot_prime + a*(33 + 64*self%ExtendedGalileon_q + 24*self%ExtendedGalileon_q**2)*Omega_tot_primeprime) + &
!      a**2*self%ExtendedGalileon_q**2*Omega_tot**2*(2*(549 + 215*self%ExtendedGalileon_q)*Omega_tot_prime**2 + a**2*(29 + 15*self%ExtendedGalileon_q)*Omega_tot_primeprime**2 + &
!        a*Omega_tot_prime*((593 + 255*self%ExtendedGalileon_q)*Omega_tot_primeprime + a*(47 + 20*self%ExtendedGalileon_q)*Omega_tot_primeprimeprime)) - &
!      6*a*self%ExtendedGalileon_q**3*Omega_tot**3*(154*Omega_tot_prime + &
!            a*(86*Omega_tot_primeprime + a*(17*Omega_tot_primeprimeprime + a*Omega_tot_primeprimeprimeprime)))) + &
!       a**14*self%ExtendedGalileon_q**3*eft_par_cache%h0_Mpc**14*Omega_tot**3*(-24*self%ExtendedGalileon_q**4*Omega_tot**4 - a**4*(24 + 46*self%ExtendedGalileon_q + 29*self%ExtendedGalileon_q**2 +& 6*self%ExtendedGalileon_q**3)*Omega_tot_prime**4 + &
!       a**3*self%ExtendedGalileon_q*(2 + self%ExtendedGalileon_q)*Omega_tot*Omega_tot_prime**2*&
!       ((71 + 34*self%ExtendedGalileon_q)*Omega_tot_prime + a*(23 + 12*self%ExtendedGalileon_q)*Omega_tot_primeprime) - &
!     a**2*self%ExtendedGalileon_q**2*Omega_tot**2*((254 + 86*self%ExtendedGalileon_q)*Omega_tot_prime**2 + a**2*(7 + 3*self%ExtendedGalileon_q)*Omega_tot_primeprime**2 + &
!       a*Omega_tot_prime*((139 + 51*self%ExtendedGalileon_q)*Omega_tot_primeprime + a*(11 + 4*self%ExtendedGalileon_q)*Omega_tot_primeprimeprime)) + &
!   a*self%ExtendedGalileon_q**3*Omega_tot**3*(154*Omega_tot_prime + &
!     a*(86*Omega_tot_primeprime + a*(17*Omega_tot_primeprimeprime + a*Omega_tot_primeprimeprimeprime))))))/&
!    ((2 + self%ExtendedGalileon_q)*eft_cache%adotoa**2 - a**2*self%ExtendedGalileon_q*eft_par_cache%h0_Mpc**2*Omega_tot)**7
    end subroutine EFTCAMBExtendedGalileonComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBExtendedGalileonAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Extended_Galileon)                :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBExtendedGalileonAdditionalModelStability       !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBExtendedGalileonAdditionalModelStability = .True.

    end function EFTCAMBExtendedGalileonAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_full_Extended_Galileon
