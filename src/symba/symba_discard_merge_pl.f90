!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_merge_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Merge planets
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-planet encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_discard_merge_pl(symba_pl1P, nplplenc, plplenc_list, ldiscard)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_merge_pl(symba_plA, nplplenc, plplenc_list, ldiscard, mergeadd_list, nmergeadd)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_merge_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT)        :: nplplenc, nmergeadd
     TYPE(symba_pl)                  :: symba_plA
     TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
     TYPE(symba_merger), INTENT(INOUT):: mergeadd_list
     LOGICAL(LGT), INTENT(INOUT)     :: ldiscard

! Internals
   INTEGER(I4B)                            :: i, j, nchild, indexchild, enc_big, index1, index2, indexk 
   REAL(DP)                                :: m, mmax, mtot, r, r3, mu, energy, ap, v2, msun, ip_merge
   real(DP)                                :: Mcb, r_child, m_child, r_big, m_big, rmerge, spin_vec_mag 
   REAL(DP), DIMENSION(NDIM)               :: x, v, vbs, xv_child, xv_big, xnew, vnew
   INTEGER(I4B), DIMENSION(NCHILDMAX)      :: array_child
   real(DP), dimension(NDIM)               :: Lspin, xc_child, xc_big, vc_child, vc_big, spin_hat, l_orb_before, l_orb_after
   real(DP), dimension(NDIM)               :: ip_child, ip_big, l_spin_after, l_spin_before, rot_child, rot_big

! Executable code
     msun = symba_plA%helio%swiftest%mass(1)
     vbs(:) = symba_plA%helio%swiftest%vb(:,1)
     DO i = 1, nplplenc
          IF (plplenc_list%status(i) == MERGED) THEN
               index1 = plplenc_list%index1(i)
               index2 = plplenc_list%index2(i)
               ! This IF statement is for handling the merger case 
               IF ((symba_plA%helio%swiftest%status(index1) == ACTIVE) .AND.                                                    &
                  (symba_plA%helio%swiftest%status(index2) == ACTIVE)) THEN

                  enc_big = plplenc_list%index1(i)

                  m = symba_plA%helio%swiftest%mass(enc_big)
                  r = symba_plA%helio%swiftest%radius(enc_big)
                  r3 = r**3
                  mmax = m
                  mtot = m
                  x(:) = m * symba_plA%helio%swiftest%xh(:,enc_big)
                  v(:) = m * symba_plA%helio%swiftest%vb(:,enc_big)
                  indexk = enc_big

                  nchild = symba_plA%nchild(enc_big)
                  array_child(1:NCHILDMAX) = symba_plA%index_child(1:NCHILDMAX,enc_big)
                  DO j = 1, nchild
                     indexchild = array_child(j)
                     m = symba_plA%helio%swiftest%mass(indexchild)
                     r = symba_plA%helio%swiftest%radius(indexchild)
                     r3 = r3 + r**3
                     mtot = mtot + m
                     x(:) = x(:) + m * symba_plA%helio%swiftest%xh(:,indexchild)
                     v(:) = v(:) + m * symba_plA%helio%swiftest%vb(:,indexchild)

                     IF (m > mmax) THEN
                        mmax = m
                        indexk = indexchild
                     END IF
                  END DO

                  x(:) = x(:)/mtot  ! position com of system 
                  v(:) = v(:)/mtot  ! velocity com of system 
                  rmerge = r3**(1.0_DP/3.0_DP)
                  ip_merge = 2.0_DP / 5.0_DP
                  ip_big = 2.0_DP / 5.0_DP
                  spin_vec_mag = 0.0_DP
                  l_orb_before(:) = 0.0_DP
                  l_spin_before(:) = 0.0_DP
                  l_orb_after(:) = 0.0_DP
                  l_spin_after(:) = 0.0_DP
                  
                  xc_big(:) = symba_plA%helio%swiftest%xh(:,enc_big) - x(:)
                  vc_big(:) = symba_plA%helio%swiftest%vb(:,enc_big) - v(:)

                  call util_crossproduct(xc_big, vc_big, xv_big)
                  ip_big = symba_plA%helio%swiftest%Ip(:, enc_big)
                  rot_big = symba_plA%helio%swiftest%rot(:, enc_big)

                  m_big = symba_plA%helio%swiftest%mass(enc_big)
                  r_big = symba_plA%helio%swiftest%radius(enc_big)
                        
                  l_orb_before = l_orb_before + (m_big * xv_big)
                  l_spin_before = l_spin_before + (ip_big * m_big * r_big**2 * rot_big)

                  DO j = 1, nchild
                     indexchild = array_child(j)
                     xc_child(:) = symba_plA%helio%swiftest%xh(:,indexchild) - x(:)
                     vc_child(:) = symba_plA%helio%swiftest%vb(:,indexchild) - v(:)

                     call util_crossproduct(xc_child, vc_child, xv_child)
                     ip_child = symba_plA%helio%swiftest%Ip(:, indexchild)
                     rot_child = symba_plA%helio%swiftest%rot(:, indexchild)

                     m_child = symba_plA%helio%swiftest%mass(indexchild)
                     r_child = symba_plA%helio%swiftest%radius(indexchild)
                        
                     l_orb_before = l_orb_before + (m_child * xv_child)
                     l_spin_before = l_spin_before + (ip_child * m_child * r_child**2 * rot_child) 
                  END DO

                  l_spin_after = l_orb_before + l_spin_before - l_orb_after
                  spin_hat = l_spin_after / NORM2(l_spin_after)
                  spin_vec_mag = NORM2(l_spin_after) / (ip_merge * mtot * rmerge**2)

                  r = rmerge
                  symba_plA%helio%swiftest%mass(indexk) = mtot
                  symba_plA%helio%swiftest%radius(indexk) = r
                  symba_plA%helio%swiftest%xh(:,indexk) = x(:)
                  symba_plA%helio%swiftest%vb(:,indexk) = v(:)
                  symba_plA%helio%swiftest%vh(:,indexk) = v(:) - vbs(:)
                  symba_plA%helio%swiftest%ip(:,indexk) = ip_child
                  symba_plA%helio%swiftest%rot(:,indexk) = spin_vec_mag*spin_hat

                  mu = msun*mtot/(msun + mtot)
                  r = SQRT(DOT_PRODUCT(x(:), x(:)))
                  v(:) = symba_plA%helio%swiftest%vh(:,indexk)
                  v2 = DOT_PRODUCT(v(:), v(:))
                  energy = -1.0_DP*msun*mtot/r + 0.5_DP*mu*v2
                  ap = -1.0_DP*msun*mtot/(2.0_DP*energy)
                  symba_plA%helio%swiftest%rhill(indexk) = ap*(((mu/msun)/3.0_DP)**(1.0_DP/3.0_DP))
                  array_child(1:NCHILDMAX) = symba_plA%index_child(1:NCHILDMAX,enc_big)
                  indexchild = enc_big
                  ldiscard = .TRUE.
                  DO j = 0, nchild
                     IF (indexchild /= indexk) THEN
                        symba_plA%helio%swiftest%status(indexchild) = MERGED
                     END IF
                     indexchild = array_child(j+1)
                  END DO

               ELSE IF ((symba_plA%helio%swiftest%status(index1) == DISRUPTION) .AND.    &                                                
                   (symba_plA%helio%swiftest%status(index2) == DISRUPTION)) THEN 
                    call symba_fragment_calculation(nmergeadd, mergeadd_list, symba_plA, plplenc_list, i) !tentatively
                    enc_big = plplenc_list%index1(i)
                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:NCHILDMAX) = symba_plA%index_child(1:NCHILDMAX,enc_big)
                    DO j = 1, nchild
                         symba_plA%helio%swiftest%status(array_child(j)) = INACTIVE
                    END DO
                    ldiscard = .TRUE.
               ELSE IF ((symba_plA%helio%swiftest%status(index1) == SUPERCATASTROPHIC) .AND.   &                                                 
                   (symba_plA%helio%swiftest%status(index2) == SUPERCATASTROPHIC)) THEN 
                    call symba_fragment_calculation(nmergeadd, mergeadd_list, symba_plA, plplenc_list, i) !tentatively
                    enc_big = plplenc_list%index1(i)
                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:NCHILDMAX) = symba_plA%index_child(1:NCHILDMAX,enc_big)
                    DO j = 1, nchild
                         symba_plA%helio%swiftest%status(array_child(j)) = INACTIVE
                    END DO
                    ldiscard = .TRUE.
               ELSE IF ((symba_plA%helio%swiftest%status(index1) == HIT_AND_RUN) .AND.      &                                              
                   (symba_plA%helio%swiftest%status(index2) == HIT_AND_RUN)) THEN 
                    call symba_fragment_calculation(nmergeadd, mergeadd_list, symba_plA, plplenc_list, i) !tentatively
                    enc_big = plplenc_list%index1(i)
                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:NCHILDMAX) = symba_plA%index_child(1:NCHILDMAX,enc_big)
                    DO j = 1, nchild
                         symba_plA%helio%swiftest%status(array_child(j)) = INACTIVE
                    END DO
                    ldiscard = .TRUE.
               ELSE IF ((symba_plA%helio%swiftest%status(index1) == GRAZE_AND_MERGE) .AND.  &    ! not used in this version, graze and merge are considered pure mergers for now (2021)                                              
                   (symba_plA%helio%swiftest%status(index2) == GRAZE_AND_MERGE)) THEN 

                    enc_big = plplenc_list%index1(i)
                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:NCHILDMAX) = symba_plA%index_child(1:NCHILDMAX,enc_big)
                    DO j = 1, nchild
                         symba_plA%helio%swiftest%status(array_child(j)) = INACTIVE
                    END DO
                    ldiscard = .TRUE.
               END IF

          END IF
     END DO
     RETURN
     
END SUBROUTINE symba_discard_merge_pl
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
