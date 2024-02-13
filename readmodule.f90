module readmodule
    implicit none
    
    contains
        subroutine locate_line(fileID, filename, linenum)
            integer, intent(in) :: fileID, linenum
            character(len=200), intent(in):: filename
            
            integer:: i
            
            rewind(fileID)
            
            do i = 1, linenum
                read(fileID, *)
            end do
            
            backspace(fileID)
            
        end subroutine locate_line
    
end module readmodule