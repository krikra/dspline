module ippe_f_interface
   use iso_c_binding
   implicit none
   interface
      function f_ippe_init_2d(n1, n2)
         use, intrinsic :: iso_c_binding, only: C_PTR
         integer :: n1, n2
         type(C_PTR) :: f_ippe_init_2d
      end function

      subroutine f_ippe_choose(ip, para)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: ip
         integer :: para
      end subroutine

      subroutine f_ippe_update(ip, para, val)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: ip
         integer :: para
         double precision :: val
      end subroutine

      function f_ippe_terminated(ip)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: ip
         integer :: f_ippe_terminated
      end function

      function f_ippe_iteration(ip)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: ip
         integer :: f_ippe_iteration
      end function

      subroutine f_ippe_destroy(ip)
         use, intrinsic :: iso_c_binding, only: C_PTR
         type(C_PTR) :: ip
      end subroutine
   end interface
end module
