module mathfun
    implicit none
    
    contains
    
        function vec_cross(veca, vecb) result(vecc)
            real*8, dimension(3), intent(in) :: veca, vecb
            real*8, dimension(3):: vecc
            
            vecc(1) = veca(2)*vecb(3) - veca(3)*vecb(2)
            vecc(2) = veca(3)*vecb(1) - veca(1)*vecb(3)
            vecc(3) = veca(1)*vecb(2) - veca(2)*vecb(1)
        
        end function vec_cross
        
        subroutine vec_norm(veca)
            real*8, dimension(3), intent(inout):: veca

            
            real*8:: norm
            
            norm = sqrt(dot_product(veca, veca))
            
            veca = veca/norm
        
        end subroutine vec_norm
        
        subroutine rotxyz(atoms_num, zvec, rvec, xyz, newxyz)
            integer, intent(in):: atoms_num
            !zvec z axis, rvec ref vec(new z axis)
            real*8, dimension(3):: zvec, rvec
            real*8, intent(in):: xyz(:, :)
            
            real*8, allocatable, intent(out):: newxyz(:,:)
            
            real*8, dimension(3):: nvec
            real*8:: theta
            real*8, dimension(3,3):: rotmat, temp1, temp2
            
            allocate(newxyz(atoms_num, 3))
            
            nvec = vec_cross(rvec, zvec)
            
            call vec_norm(zvec)
            call vec_norm(rvec)
            call vec_norm(nvec)
            
            
            theta = acos(dot_product(zvec, rvec))
            !vec(1,3) * rotmat(3,3)^T = vec(1,3)
            !rotmat(3,3) *vec(3,1) = vec(3,1)
            rotmat(1,:) = (/ nvec(1)*nvec(1)*(1 - cos(theta)) + cos(theta),&
                             nvec(1)*nvec(2)*(1 - cos(theta)) - nvec(3)*sin(theta),&
                             nvec(1)*nvec(3)*(1 - cos(theta)) + nvec(2)*sin(theta)/)
                             
            rotmat(2,:) = (/ nvec(2)*nvec(1)*(1 - cos(theta)) + nvec(3)*sin(theta),&
                             nvec(2)*nvec(2)*(1 - cos(theta)) + cos(theta),&
                             nvec(2)*nvec(3)*(1 - cos(theta)) - nvec(1)*sin(theta)/)
                             
            rotmat(3,:) = (/ nvec(3)*nvec(1)*(1 - cos(theta)) - nvec(2)*sin(theta),&
                             nvec(3)*nvec(2)*(1 - cos(theta)) + nvec(1)*sin(theta),&
                             nvec(3)*nvec(3)*(1 - cos(theta)) + cos(theta)/)
            
            temp1 = matmul(rotmat, transpose(rotmat))
            temp2 = matmul(transpose(rotmat), rotmat)
            newxyz = matmul(xyz, transpose(rotmat))
        
        end subroutine rotxyz
    
end module mathfun