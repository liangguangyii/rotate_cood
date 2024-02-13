program main
    use readmodule
    use mathfun     
    implicit none
    
    character(len=200):: filename, temp_str
    character(len=2), allocatable:: atom_list(:)
    integer:: atoms_num, line_num, atomstart, atomend, i, j
    real*8, dimension(3):: zvec, rvec
    real*8, allocatable:: xyz(:,:), newxyz(:,:)
    
    filename = "VinAu12_NMR.gjf"
    atoms_num = 13
    line_num = 10
    
    atomstart = 13
    atomend = 1
    
    zvec=(/0D0, 0D0, 1D0/)
    
    open(10, file=filename, status='old')
    
        allocate(atom_list(atoms_num))
        allocate(xyz(atoms_num, 3))

        call locate_line(10, filename, line_num)
        
        do i = 1, atoms_num
            read(10,*) atom_list(i), xyz(i,1), xyz(i,2), xyz(i,3)
        end do
    close(10)
    
    !write(*, '(1a2)') (atom_list(i), i=1, atoms_num)
    !write(*, '(3f12.7)') ((xyz(i, j), j=1, 3), i=1, atoms_num) 
    
    !write(*, *) vec_cross((/1D0, 0D0, 0D0/), (/0D0, 1D0, 0D0/))
    !write(*, *) vec_norm((/3D0, 4D0, 0D0/))
    allocate(newxyz(atoms_num, 3))
    rvec = xyz(atomstart, :) - xyz(atomend, :)
    
    call rotxyz(atoms_num, zvec, rvec, xyz, newxyz)
    
    open(11, file="newxyz.txt", status='replace')
        write(11, '(1a2,4x,3f12.7)') (atom_list(i), (newxyz(i, j), j=1, 3), i=1, atoms_num)
    close(11)
    
end program main

    
    