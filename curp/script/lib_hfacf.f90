
subroutine cal_hfacf(acf, xss, nacf, first, last, interval, shift, &
                   & norm, nsample, ndim, nframe, ncom)

   implicit none
   integer, intent(in) :: nacf, first, last, interval
   integer, intent(in) :: shift, nframe, ncom, nsample, ndim
   real(8), intent(in) :: xss(nframe, ndim, ncom)
   logical, optional, intent(in) :: norm
   real(8), intent(out):: acf(nacf, ncom)

   integer :: iacf, ifrm, ifrm_beg, nfrm_beg, icom, nsam
   real(8) :: a

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

            ifrm = ifrm_beg + (iacf - 1) * interval
            a = a + xss(ifrm_beg,1,icom) * xss(ifrm,1,icom) &
                  + xss(ifrm_beg,2,icom) * xss(ifrm,2,icom) &
                  + xss(ifrm_beg,3,icom) * xss(ifrm,3,icom)

         end do
  
         acf(iacf,icom) = a

      end do

   end do ! samples of one file data

   nsam = (nfrm_beg-first)/shift + 1
   acf(:,:) = acf(:,:) / real(nsam)
end subroutine
