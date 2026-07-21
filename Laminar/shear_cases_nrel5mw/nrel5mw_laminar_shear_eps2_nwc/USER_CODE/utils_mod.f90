! Module that contains some useful types and procedures
! S. Ouchene
! Date: 19/08/2025
module utils_mod
    implicit none
    
    ! We create a type to be able to use an allocatable array of non-fixed length strings
    ! See details here:
    !https://fortran-lang.discourse.group/t/how-do-i-allocate-an-array-of-strings/3930
    type var_string_type
        character(len=:), allocatable :: filename
    end type

    contains

    subroutine check_file(file_path)
        character(len=*),  intent(in) :: file_path
        logical :: exists_

        inquire(file=trim(file_path), exist=exists_)

        if (.NOT. exists_) then
            error stop "Required file not found: " // trim(file_path)
        end if
            
    end subroutine check_file
end module utils_mod
