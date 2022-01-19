
subroutine cal_hfacf(acf, xss, nacf, first, last, interval, shift, &
                   & norm, nsample, ndim, nframe, ncom)

   implicit none
   integer, intent(in) :: nacf, first, last, interval
   integer, intent(in) :: shift, nframe, ncom, nsample, ndim
   real(8), intent(in) :: xss(nframe, ndim, ncom)
   logical, optional, intent(in) :: norm
   real(8), intent(out):: acf(nacf, ncom)

   real(8) :: x2_sum_inv(ncom), acf_tmp(nacf, ncom), temp(ncom)
   integer :: iacf, ifrm, ifrm_beg, nfrm_beg, icom, nsam, i
   logical :: use_norm
   real(8) :: a, norm_flg, b(ncom)

   ! determine the default value of use_norm
   ! if ( norm .eqv. .true. ) then
   !    use_norm = norm
   !    norm_flg = 1.0d0
   ! else
   !    use_norm = .false.
   !    norm_flg = 0.0d0
   ! end if

   ! Normalization is disabled
   norm_flg = 0.0d0

   ! determine the number of samples
   if ( nsample <= 0 ) then
      nfrm_beg = nframe - (last-first)
   else if ( nsample >= nframe - (last-first) ) then
      nfrm_beg = nframe - (last-first)
   else
      nfrm_beg = first + shift*(nsample-1)
   end if

   do icom=1, ncom

      do iacf=1, nacf
         a = 0.0d0 
         b(icom)=0.0d0
         do ifrm_beg=first, nfrm_beg, shift

            ! if use_norm then x2_sum_inv = 1/h(0)h(0) else x2_sum_inv = 1.0d0            
            
            x2_sum_inv(icom) = 1.0d0/(xss(ifrm_beg,1,icom)**2 &
                                    + xss(ifrm_beg,2,icom)**2 &
                                    + xss(ifrm_beg,3,icom)**2)
            x2_sum_inv(icom) = x2_sum_inv(icom) - 1.0d0
            x2_sum_inv(icom) = (x2_sum_inv(icom) * norm_flg) + 1.0d0

            ifrm = ifrm_beg+(iacf-1)*interval
            a = a + ( xss(ifrm_beg,1,icom)*xss(ifrm,1,icom) &
                    + xss(ifrm_beg,2,icom)*xss(ifrm,2,icom) &
                    + xss(ifrm_beg,3,icom)*xss(ifrm,3,icom) ) * x2_sum_inv(icom)
            b(icom) = b(icom) + (1.0d0/x2_sum_inv(icom))

         end do
  
         acf(iacf,icom) = a

      end do

   end do ! samples of one file data

   nsam = (nfrm_beg-first)/shift + 1
   ! print*, first, nfrm_beg, nsam
   acf(:,:) = acf(:,:) / real(nsam)
   b(:) = b(:) / real(nsam)
   ! print*, '<h(0)h(0)> = ',b(1)
end subroutine
