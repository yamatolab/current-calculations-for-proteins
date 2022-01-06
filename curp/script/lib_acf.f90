
subroutine cal_acf(acf, xss, nacf, first, last, interval, shift, &
                 & norm, nsample, nframe, ncom)

   implicit none
   integer, intent(in) :: nacf, first, last, interval
   integer, intent(in) :: shift, nframe, ncom, nsample
   real(8), intent(in) :: xss(nframe, ncom)
   logical, optional, intent(in) :: norm
   real(8), intent(out):: acf(nacf, ncom)

   real(8) :: x2_sum_inv(ncom), acf_tmp(nacf, ncom)
   integer :: iacf, ifrm, ifrm_beg, nfrm_beg, icom, nsam
   logical :: use_norm
   real(8) :: a, norm_flg, b(ncom)

   ! determine the number of samples
   if ( nsample <= 0 ) then
      nfrm_beg = nframe - (last-first)
   else if ( nsample >= nframe - (last-first) ) then
      nfrm_beg = nframe - (last-first)
   else
      nfrm_beg = first + shift*(nsample-1)
   end if

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
 
   ! << Calculate auto-correlation function >>
   do icom=1, ncom
      do iacf=1, nacf
         a=0.0d0
         b(icom)=0.0d0
         do ifrm_beg=first, nfrm_beg, shift

            ! if use_norm then x2_sum_inv = 1/J(0)J(0) else x2_sum_inv = 1.0d0
            x2_sum_inv(icom) = 1.0d0/(xss(ifrm_beg,icom)*xss(ifrm_beg,icom))
            x2_sum_inv(icom) = x2_sum_inv(icom) - 1.0d0
            x2_sum_inv(icom) = (x2_sum_inv(icom) * norm_flg) + 1.0d0

            ifrm = ifrm_beg + (iacf-1)*interval
            a = a + ( ( xss(ifrm_beg,icom) * xss(ifrm,icom) ) * x2_sum_inv(icom) )
            b(icom) = b(icom) + (1.0d0/x2_sum_inv(icom))

         end do
         acf(iacf,icom) = a
      end do
   end do ! samples of one file data

   nsam = (nfrm_beg-first)/shift + 1
   ! print*, first, nfrm_beg, nsam
   acf(:,:) = acf(:,:) / real(nsam)
   b(:) = b(:) / real(nsam)
end subroutine

