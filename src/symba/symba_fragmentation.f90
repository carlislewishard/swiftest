!**********************************************************************************************************************************
!
!  Unit Name   : symba_fragmentation
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
!  Invocation  : CALL symba_fragmentation(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_fragmentation (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
     encounter_file, npl, symba_plA, nplplenc, plplenc_list, plmaxname, tpmaxname, mtiny, lfragmentation)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_fragmentation
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc
     INTEGER(I4B), INTENT(IN)                         :: npl, nplplenc
     INTEGER(I4B), INTENT(INOUT)                      :: plmaxname, tpmaxname, nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset, mtiny
     REAL(DP), DIMENSION(:), INTENT(IN)               :: vbs
     CHARACTER(*), INTENT(IN)                         :: encounter_file
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     LOGICAL(LGT), INTENT(IN)                         :: lfragmentation

! Internals
 
     INTEGER(I4B)                   :: model, nres, i, itarg, iproj
     REAL(DP), DIMENSION(NDIM)      :: mres, rres
     REAL(DP), DIMENSION(NDIM, 3)   :: pres, vres
     INTEGER(I4B)                   :: regime 
     INTEGER(I4B)                   :: index1, index2, index1_child, index2_child, index1_parent, index2_parent
     INTEGER(I4B)                   :: name1, name2, index_big1, index_big2, stat1, stat2, nchild1, nchild2
     REAL(DP)                       :: r2, rlim, rlim2, vdotr, tcr2, dt2, a, e, q
     REAL(DP)                       :: rad1, rad2, m1, m2, den1, den2, vol1, vol2, vchild, dentarg, denproj, dentot, Mcenter
     REAL(DP)                       :: mass1, mass2, mmax, mtmp, mtot, m1_si, m2_si
     REAL(DP), DIMENSION(NDIM)      :: xr, vr, x1, v1, x2, v2, x1_si, x2_si, v1_si, v2_si, xproj, xtarg, vproj, vtarg, vbs_si
     REAL(DP)                       :: den1_si, den2_si, rad1_si, rad2_si, rproj, rtarg, mtiny_si
     LOGICAL(LGT)                   :: lfrag_add, lmerge
     INTEGER(I4B), DIMENSION(:), allocatable   :: array_index1_child, array_index2_child
     REAL(DP)                       :: Mlr, Mslr, mtarg, mproj
     REAL(DP), DIMENSION(NDIM)      ::  denvec
     

! Executable code

     ! Recalculates vbs 
     CALL coord_vb2vh(npl, symba_plA%helio%swiftest)

     lmerge = .FALSE.
     lfrag_add = .FALSE.
     ! Model 2 is the model for collresolve_resolve (LS12)
     model = 2
     index1 = plplenc_list%index1(index_enc)
     index2 = plplenc_list%index2(index_enc)

     rlim = symba_plA%helio%swiftest%radius(index1) + symba_plA%helio%swiftest%radius(index2)
     xr(:) = symba_plA%helio%swiftest%xh(:,index2) - symba_plA%helio%swiftest%xh(:,index1)
     r2 = DOT_PRODUCT(xr(:), xr(:))
     rlim2 = rlim**2
     ! checks if bodies are actively colliding in this time step
     IF (rlim2 >= r2) THEN 
          lfrag_add = .TRUE.
     ! if they are not actively colliding in  this time step, 
     !checks if they are going to collide next time step based on velocities and q
     ELSE 
          vr(:) = symba_plA%helio%swiftest%vb(:,index2) - symba_plA%helio%swiftest%vb(:,index1)
          vdotr = DOT_PRODUCT(xr(:), vr(:))
          IF (plplenc_list%lvdotr(index_enc) .AND. (vdotr > 0.0_DP)) THEN 
               tcr2 = r2 / DOT_PRODUCT(vr(:), vr(:))
               dt2 = dt**2
               IF (tcr2 <= dt2) THEN
                    mtot = symba_plA%helio%swiftest%mass(index1) + symba_plA%helio%swiftest%mass(index2)
                    CALL orbel_xv2aeq(xr(:), vr(:), mtot, a, e, q)
                    IF (q < rlim) lfrag_add = .TRUE.
               END IF
               ! if no collision is going to happen, write as close encounter, not  merger
               IF (.NOT. lfrag_add) THEN
                    IF (encounter_file /= "") THEN
                         name1 = symba_plA%helio%swiftest%name(index1)
                         m1 = symba_plA%helio%swiftest%mass(index1)
                         rad1 = symba_plA%helio%swiftest%radius(index1)
                         x1(:) = symba_plA%helio%swiftest%xh(:,index1)
                         v1(:) = symba_plA%helio%swiftest%vb(:,index1) - vbs(:)
                         name2 = symba_plA%helio%swiftest%name(index2)
                         m2 = symba_plA%helio%swiftest%mass(index2)
                         rad2 = symba_plA%helio%swiftest%radius(index2)
                         x2(:) = symba_plA%helio%swiftest%xh(:,index2)
                         v2(:) = symba_plA%helio%swiftest%vb(:,index2) - vbs(:)

                         CALL io_write_encounter(t, name1, name2, m1, m2, rad1, rad2, x1(:), x2(:), &
                              v1(:), v2(:), encounter_file)
                    END IF
               END IF
          END IF
     END IF

     nres = 2
     IF (lfrag_add) THEN 
          symba_plA%lmerged(index1) = .TRUE.
          symba_plA%lmerged(index2) = .TRUE.
          index1_parent = symba_plA%index_parent(index1)
          m1 = symba_plA%helio%swiftest%mass(index1_parent)
          mass1 = m1 
          x1(:) = m1 * symba_plA%helio%swiftest%xh(:,index1_parent)
          v1(:) = m1 * symba_plA%helio%swiftest%vb(:,index1_parent)
          nchild1 = symba_plA%nchild(index1_parent)  

          mmax = m1
          name1 = symba_plA%helio%swiftest%name(index1_parent)
          index_big1 = index1_parent
          stat1 = symba_plA%helio%swiftest%status(index1_parent)
          vol1 =  (4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(index1_parent)**3
          if (nchild1 > 0) then
            allocate(array_index1_child(nchild1))
            array_index1_child(:) = symba_plA%index_child(1:nchild1, index1_parent)
            DO i = 1, nchild1 ! initialize an array of children
               index1_child = array_index1_child(i)
               mtmp = symba_plA%helio%swiftest%mass(index1_child)
               vchild = (4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(index1_child)**3
               vol1 = vol1 + vchild
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    name1 = symba_plA%helio%swiftest%name(index1_child)
                    index_big1 = index1_child
                    stat1 = symba_plA%helio%swiftest%status(index1_child)
               END IF
               m1 = m1 + mtmp
               x1(:) = x1(:) + mtmp * symba_plA%helio%swiftest%xh(:, index1_child)
               v1(:) = v1(:) + mtmp * symba_plA%helio%swiftest%vb(:, index1_child)
            END DO
          else
            allocate(array_index1_child(1))
          end if
          den1 =  m1 / vol1
          rad1 = ((3 * m1) / (den1 * 4 * PI))**(1.0_DP / 3.0_DP)
          x1(:) = x1(:) / m1
          v1(:) = v1(:) / m1

          index2_parent = symba_plA%index_parent(index2)
          m2 = symba_plA%helio%swiftest%mass(index2_parent)
          mass2 = m2
          rad2 = symba_plA%helio%swiftest%radius(index2_parent)
          x2(:) = m2 * symba_plA%helio%swiftest%xh(:,index2_parent)
          v2(:) = m2 * symba_plA%helio%swiftest%vb(:,index2_parent)
          mmax = m2
          name2 = symba_plA%helio%swiftest%name(index2_parent)
          index_big2 = index2_parent
          stat2 = symba_plA%helio%swiftest%status(index2_parent)
          nchild2 = symba_plA%nchild(index2_parent)  

          vol2 = (4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(index2_parent)**3

          if (nchild2 > 0) then
            allocate(array_index2_child(nchild2))
            array_index2_child(:) = symba_plA%index_child(1:nchild2, index2_parent)
            DO i = 1, nchild2
               index2_child = array_index2_child(i)
               mtmp = symba_plA%helio%swiftest%mass(index2_child)
               vchild =  (4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(index2_child)**3
               vol2 = vol2 + vchild
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    name2 = symba_plA%helio%swiftest%name(index2_child)
                    index_big2 = index2_child
                    stat2 = symba_plA%helio%swiftest%status(index2_child)
               END IF
               m2 = m2 + mtmp
               x2(:) = x2(:) + mtmp * symba_plA%helio%swiftest%xh(:, index2_child)
               v2(:) = v2(:) + mtmp * symba_plA%helio%swiftest%vb(:, index2_child)
            END DO
          else
            allocate(array_index2_child(1))
          end if 

          if (lfragmentation) then ! Determine the collisional regime and resolve the collision
            den2 =  m2 / vol2
            rad2 = ((3 * m2) / (den2 * 4 * PI))**(1.0_DP / 3.0_DP)
            x2(:) = x2(:) / m2
            v2(:) = v2(:) / m2
   
            m1_si = (m1 / GU) * MU2KG 
            m2_si = (m2 / GU) * MU2KG
            rad1_si = rad1 * DU2M
            rad2_si = rad2 * DU2M
            x1_si(:) = x1(:) * DU2M
            x2_si(:) = x2(:) * DU2M
            v1_si(:) = v1(:) * DU2M / TU2S
            v2_si(:) = v2(:) * DU2M / TU2S
            den1_si = (den1 / GU) * MU2KG / DU2M**3
            den2_si = (den2 / GU) * MU2KG / DU2M**3
            vbs_si = vbs(:) * DU2M / TU2S 
   
            mres(:) = 0.0_DP
            rres(:) = 0.0_DP
            pres(:,:) = 0.0_DP
            vres(:,:) = 0.0_DP
   
            IF (m1_si > m2_si) THEN 
                  itarg = index1
                  iproj = index2
                  dentarg = den1_si
                  denproj = den2_si
                  mtarg = m1_si
                  mproj = m2_si
                  rtarg = rad1_si
                  rproj = rad2_si
                  xtarg(:) = x1_si(:)
                  xproj(:) = x2_si(:)
                  vtarg(:) = v1_si(:)
                  vproj(:) = v2_si(:)
            ELSE
                  itarg = index2
                  iproj = index1
                  dentarg = den2_si
                  denproj = den1_si
                  mtarg = m2_si
                  mproj = m1_si
                  rtarg = rad2_si
                  rproj = rad1_si
                  xtarg(:) = x2_si(:)
                  xproj(:) = x1_si(:)
                  vtarg(:) = v2_si(:)
                  vproj(:) = v1_si(:)
            END IF
            mtot = m1_si + m2_si
            dentot = (m1_si * den1_si + m2_si * den2_si)/ mtot
            Mcenter = symba_plA%helio%swiftest%mass(1) * MU2KG / GU
            mtiny_si = (mtiny / GU) * MU2KG
            !regime = collresolve_resolve(model,mtarg,mproj,rtarg,rproj,xtarg,xproj, vtarg,vproj, nres, mres, rres, pres, vres)
            CALL util_regime(Mcenter, mtarg, mproj, rtarg, rproj, xtarg, xproj, vtarg, vproj, dentarg, denproj, &
                  regime, Mlr, Mslr, mtiny_si)

            mres(1) = min(max(Mlr, 0.0_DP), mtot)
            mres(2) = min(max(Mslr, 0.0_DP), mtot)
            mres(3) = min(max(mtot - Mlr - Mslr, 0.0_DP), mtot)
            denvec(1) = dentarg
            denvec(2) = denproj
            denvec(3) = dentot

            rres(:) = (3 * mres(:)  / (4 * PI * denvec(:)))**(1.0_DP / 3.0_DP)
   
            mres(:) = (mres(:) / MU2KG) * GU
            rres(:) = rres(:) / DU2M
         else ! Only do a pure merger
            regime = COLLRESOLVE_REGIME_MERGE
         end if

         CALL symba_caseresolve(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs_si, & 
                                symba_plA, nplplenc, plplenc_list, regime, plmaxname, tpmaxname, mres, rres, array_index1_child, &
                                array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2, mtiny)
     END IF 
     RETURN

END SUBROUTINE symba_fragmentation
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
