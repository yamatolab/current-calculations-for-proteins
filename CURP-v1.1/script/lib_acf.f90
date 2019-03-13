
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

   ! determine the default value of use_norm
   if ( present(norm) ) then
      use_norm = norm
   else
      use_norm = .false.
   end if

   ! determine the number of samples
   if ( nsample <= 0 ) then
      nfrm_beg = nframe - (last-first)
   else if ( nsample >= nframe - (last-first) ) then
      nfrm_beg = nframe - (last-first)
   else
      nfrm_beg = first + shift*(nsample-1)
   end if

   ! << Calculate auto-correlation function >>
   do ifrm_beg=first, nfrm_beg, shift

      ! x**2 for norm
      if (use_norm) then
         x2_sum_inv(:) = 1.0d0/(xss(ifrm_beg,:)*xss(ifrm_beg,:))
      else
         x2_sum_inv(:) = 1.0d0
      end if

      ! calculate acf
      acf_tmp(:,:) = 0.0d0
      do iacf=1, nacf
         ifrm = ifrm_beg + (iacf-1)*interval
         acf_tmp(iacf,:) = xss(ifrm_beg,:) * xss(ifrm,:)
      end do

      ! sum up acf by the number of samples
      do icom=1, ncom
         acf(:,icom) = acf(:,icom) + acf_tmp(:,icom)*x2_sum_inv(icom)
      end do

      ! integral acf
      ! print*, '#lib ', sum(acf_tmp(:,1))

   end do ! samples of one file data

   nsam = (nfrm_beg-first)/shift + 1
   ! print*, first, nfrm_beg, nsam
   acf(:,:) = acf(:,:) / real(nsam)

end subroutine

