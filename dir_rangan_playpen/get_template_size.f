ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_template_size(nlats,ncur,ntemplatesize,ngridc,
     1           icstart)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Generate quasiuniform grid on concentric circles, with
c     the number of points on each circle at radius k corresponding 
c     to the number for the great circle on the corresponding sphere
c     of radius k (encoded n nlats).
c
c     INPUT:
c     nlats()      number of quarature nodes in theta on sphere
c                  defined by index ncur.
c     ncur         index of highest frequency sphere under consideration
c
c     OUT:         
c     ngridc()       number of output points on successive circles 
c                    in templates
c     ntemplatesize  total number of pts in templates 
c     icstart()      indexing array for points on successive circles 
c
      implicit none
      integer ntemplatesize,nlats(ncur),ncur,ngridc(ncur),icstart(ncur)
      integer i
c
      icstart(1) = 1
      do i = 1,ncur-1
         ngridc(i) = nlats(i)*2     
         icstart(i+1) = icstart(i) + ngridc(i)
      enddo
      ngridc(ncur) = nlats(ncur)*2     
      ntemplatesize = icstart(ncur) + ngridc(ncur) - 1
      return
      end
