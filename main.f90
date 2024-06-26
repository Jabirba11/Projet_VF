Program euler
    Use mod_fonctions
    Use mod_schemas  
    Use mod_output
    Use mod_test

    Implicit None

    Character(14), Parameter :: parameters = "parameters.dat"


! Case ('IsentropicVortex')
!     ! Periodic
!     Select Case (axis)
!     Case ('x')
!         If (i == 0 .OR. i == imax) Then
!             ULL = Uvect(:,imax-1,j)
!             UL = Uvect(:,imax,j)
!             UR = Uvect(:,1,j)
!             URR = Uvect(:,2,j)
!         Else If (i == 1) Then
!             ULL = Uvect(:,imax,j)
!             UL = Uvect(:,1,j)
!             UR = Uvect(:,2,j)
!             URR = Uvect(:,3,j)
!         Else If (i == imax-1) Then
!             ULL = Uvect(:,imax-2,j)
!             UL = Uvect(:,imax-1,j)
!             UR = Uvect(:,imax,j)
!             URR = Uvect(:,1,j)
!         Else
!             Write(STDERR, *) "Invalid index for boundary condition ", i
!             Call Exit(1)
!         End If
!     Case ('y')
!         If (j == 0 .OR. j == jmax) Then
!             ULL = Uvect(:,i,jmax-1)
!             UL = Uvect(:,i,jmax)
!             UR = Uvect(:,i,1)
!             URR = Uvect(:,i,2)
!         Else If (j == 1) Then
!             ULL = Uvect(:,i,jmax)
!             UL = Uvect(:,i,1)
!             UR = Uvect(:,i,2)
!             URR = Uvect(:,i,3)
!         Else If (j == jmax-1) Then
!             ULL = Uvect(:,i,jmax-2)
!             UL = Uvect(:,i,jmax-1)
!             UR = Uvect(:,i,jmax)
!             URR = Uvect(:,i,1)
!         Else
!             Write(STDERR, *) "Invalid index for boundary condition ", j
!             Call Exit(1)
!         End If

    ! Parametres
    Real(PR)          :: xmin, xmax, ymin, ymax, time_max, cfl, gamma
    Integer           :: imax, jmax, output_modulo, case
    Character(len=10) :: numflux_name
   
    ! Input related variables
    Character(len=100) :: buffer, label
    Integer            :: pos
    Integer            :: ios = 0
    Integer            :: line_number = 0
    
    ! Arrays
    Real(PR), Dimension(:,:,:), Allocatable :: Uvect, Uvect_exacte, vect_fluxF, vect_fluxG
    Real(PR), Dimension(:,:,:), Allocatable :: U_minus_x,U_plus_x
    Real(PR), Dimension(:,:,:), Allocatable :: U_minus_y,U_plus_y
    Real(PR), Dimension(:,:,:), Allocatable :: U_RK1,U_RK2
 

    Real(PR), Dimension(:), Allocatable     :: x, y, xm, ym
    




    
    ! Loop indices
    Integer  :: i, j, nb_iterations
    ! Other
    Real(PR) :: deltax, deltay, deltat, time


    ! Read parameters
    Open(111, File=parameters)
    Do While (ios == 0)
        Read(111, '(A)', IOstat=ios) buffer
        If (ios == 0) Then
            line_number = line_number + 1

            pos = SCAN(buffer, ' ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            Select Case (label)
            Case ('xmin')
                Read(buffer, *, iostat=ios) xmin
            Case ('xmax')
                Read(buffer, *, iostat=ios) xmax
            Case ('ymin')
                Read(buffer, *, iostat=ios) ymin
            Case ('ymax')
                Read(buffer, *, iostat=ios) ymax
            Case ('Nx')
                Read(buffer, *, iostat=ios) imax
            Case ('Ny')
                Read(buffer, *, iostat=ios) jmax
            Case ('tmax')
                Read(buffer, *, iostat=ios) time_max
            Case ('CFL')
                Read(buffer, *, iostat=ios) cfl
            Case ('gamma')
                Read(buffer, *, iostat=ios) gamma
            Case ('output_modulo')
                Read(buffer, *, iostat=ios) output_modulo
            Case ('case')
                Read(buffer, *, iostat=ios) case
            Case ('flux')
                Read(buffer, *, iostat=ios) numflux_name
            Case ('')
                ! Do nothing if it is an empty line
            Case Default
                Write(STDOUT,*) "Invalid label", label, " at line", line_number, "(skipping)"
            End Select
        End If
    End Do
    Close(111)
    
    Allocate(x(0:imax), y(0:jmax), xm(imax), ym(jmax))
    Allocate(Uvect(4,imax,jmax), Uvect_exacte(4,imax,jmax), vect_fluxF(4,0:imax, 0:jmax), vect_fluxG(4,0:imax, 0:jmax))
    Allocate(U_minus_x(4,imax,jmax),U_plus_x(4,imax,jmax))
    Allocate(U_minus_y(4,imax,jmax),U_plus_y(4,imax,jmax))
    Allocate(U_RK1(4,imax,jmax),U_RK2(4,imax,jmax))



    deltax = (xmax - xmin) / imax
    deltay = (ymax - ymin) / jmax
    x = (/ Real(PR) :: (xmin + i*deltax, i=0, imax) /)
    y = (/ Real(PR) :: (ymin + j*deltay, j=0, jmax) /)
    xm = (/ Real(PR) :: (xmin + .5_PR*deltax + i*deltax, i=0, imax-1) /)
    ym = (/ Real(PR) :: (ymin + .5_PR*deltay + j*deltay, j=0, jmax-1) /)

    !          ! DIRICHLET
!          !F(0,j) = -2._PR*D(xmin, ym(j)) * (T(1, j) - Touest(tn, ym(j))) / deltax
!          !F(imax,j) = 2._PR*D(xmax, ym(j)) * (T(imax,j) - Test(tn,ym(j))) / deltax
!          ! NEUMANN
!          !fluxF(:,imax,j) = (/ Real(PR) :: 0._PR, 0._PR, 0._PR, 0._PR /)
!         !fluxF(:,0,j) = (/ Real(PR) :: 0._PR, 0._PR, 0._PR, 0._PR /)
!       End Do
!       Do i=1, imax
!          Do j=1, jmax-1
!             fluxG(:,i,j) = Rusanov(Uvect(:,i,j), Uvect(:,i,j+1), gamma, fluxFuncG)
!             !fluxG(:,i,j) = HLL(Uvect(:,i,j), Uvect(:,i,j+1), gamma, fluxFuncG)
!          End Do
!          ! DIRICHLET
!          !G(i,0) = -2._PR*D(xm(i), ymin) * ( T(i,1) - Tnord(tn,xm(i)) ) / deltay
!          !G(i,jmax) = 2._PR*D(xm(i), ymax) * (T(i,jmax) - Tsud(tn,xm(i))) / deltay
!          ! NEUMANN
!          !fluxG(:,i,0) = (/ Real(PR) :: 0._PR, 0._PR, 0._PR, 0._PR /)
!          !fluxG(:,i,jmax) = (/ Real(PR) :: 0._PR, 0._PR, 0._PR, 0._PR /)
!       End Do
      
    
    
    !Initialisation de U
    
    Do i=1, imax
        Do j=1, jmax
            Uvect(:,i,j) = Uinit(xm(i), ym(j), gamma, case)
        End Do
    End Do


    !Call PerformTestsAndExit(gamma)

    ! ! Output initial state
    ! Call output_csv(Uvect, gamma, x, y, 0, 'sol')
    ! Call output_csv(Uvect, gamma, x, y, 0, 'exact')

    ! Time loop
    nb_iterations = 0
    time = 0._PR
    Do While (time < time_max)
        Call compute_CFL(Uvect, deltax, deltay, deltat, cfl)

        time = MIN( time + deltat, time_max )

    
        
        vect_fluxF = Compute_Flux('x',numflux_name,Uvect,case) 
        vect_fluxG = Compute_Flux('y',numflux_name,Uvect,case)

        
        !---------------------Stage 1 -------------------------
        
        Do i=1, imax
            Do j=1, jmax                
                
                U_RK1(:,i,j) = Uvect(:,i,j) &
                & - deltat/deltax * (vect_fluxF(:,i,j) - vect_fluxF(:,i-1,j)) &
                & - deltat/deltay * (vect_fluxG(:,i,j) - vect_fluxG(:,i,j-1))
                
            End Do
        End Do
        
        


        
        vect_fluxF = Compute_Flux('x',numflux_name,U_RK1,case) 
        vect_fluxG = Compute_Flux('y',numflux_name,U_RK1,case)
        
        

    !  !---------------------Stage 2 -------------------------

        
        Do i=1, imax
            Do j=1, jmax
                          
                U_RK2(:,i,j) = (3._PR/4._PR)*Uvect(:,i,j) + (1._PR/4._PR)*U_RK1(:,i,j) &
                & - (1._PR*deltat)/(4._PR*deltax) * (vect_fluxF(:,i,j) - vect_fluxF(:,i-1,j)) &
                & - (1._PR*deltat)/(4._PR*deltay) * (vect_fluxG(:,i,j) - vect_fluxG(:,i,j-1))
                
            End Do
        End Do


        vect_fluxF = Compute_Flux('x',numflux_name,U_RK2,case) 
        vect_fluxG = Compute_Flux('y',numflux_name,U_RK2,case)

  
        
        Do i=1, imax
            Do j=1, jmax        
                
                Uvect(:,i,j) = (1._PR/3._PR)*Uvect(:,i,j) + (2._PR/3._PR)*U_RK2(:,i,j) &
                & - (2._PR*deltat)/(3._PR*deltax) * (vect_fluxF(:,i,j) - vect_fluxF(:,i-1,j)) &
                & - (2._PR*deltat)/(3._PR*deltay) * (vect_fluxG(:,i,j) - vect_fluxG(:,i,j-1))
                
            End Do
        End Do





        If ( output_modulo > 0 .AND. Modulo(nb_iterations, output_modulo) == 0 ) Then
            Write(STDOUT, *) time, time_max
            Do i=1, imax
                Do j=1, jmax
                    Uvect_exacte(:,i,j) = Uexact(case, xm(i), ym(j), time, gamma)
                End Do
            End Do
            ! Call output_csv(Uvect, gamma, x, y, nb_iterations / output_modulo + 1, 'sol')
            ! Call output_csv(Uvect_exacte, gamma, x, y, nb_iterations / output_modulo + 1, 'exact')
        End If
        
        
        
        nb_iterations = nb_iterations + 1

 
        
    End Do

    Write(STDOUT, *) "Error:", error('L1', case, Uvect, time, gamma) 

    Deallocate(x, y, xm, ym)
    Deallocate(Uvect, Uvect_exacte, vect_fluxF, vect_fluxG)
    Deallocate(U_minus_x,U_plus_x)
    Deallocate(U_minus_y,U_plus_y)
    Deallocate(U_RK1,U_RK2)


Contains
    Subroutine compute_CFL(U, dx, dy, dt, cfl)
        
        Real(PR), Dimension(:,:,:), Intent(In) :: U
        Real(PR), Intent(In) :: dx, dy,  cfl
        Real(PR), Intent(Out) :: dt
    
        Real(PR) :: rho, velocity_u, velocity_v, e, q, p, a, bx, by, bx_max, by_max, l1, l3
        Integer  :: i,j

        bx_max = 0._PR
        by_max = 0._PR
        Do i=1, imax
            Do j=1 , jmax
                
                rho        = U(1,i,j)
                velocity_u = U(2,i,j) / rho
                velocity_v = U(3,i,j) / rho
                e          = U(4,i,j)

                q = .5_PR * ( velocity_u**2 + velocity_v**2 )
                p = (gamma - 1._PR)*(e - rho*q)
                a = SQRT(gamma*p/rho)

                l1 = ABS(velocity_u - a)
                l3 = ABS(velocity_u + a)
                bx = MAX(l1, l3)
                bx_max = MAX( bx, bx_max)
                
                l1 = ABS(velocity_v - a)
                l3 = ABS(velocity_v + a)
                by = MAX(l1, l3)
                by_max = MAX( by, by_max)
            
            End Do
        End Do

        dt = cfl * 0.5_PR * MIN(dx/bx_max, dy/by_max)
    End Subroutine compute_CFL

    Function error(norm, case, U, time, gamma)
        ! --- InOut
        Real(PR), Dimension(4,imax,jmax), Intent(In) :: U
        Real(PR), Intent(In) :: gamma, time
        Integer, Intent(In) :: case
        Character(len=*), Intent(In) :: norm
        Real(PR), Dimension(4) :: error
        ! --- Locals
        Real(PR), Dimension(4) :: exact_value
        Integer                :: i,j,k
        ! Real(PR)               :: error

        error = 0._PR
        ! sum_error = 0._PR
        Select Case (norm)
        Case ('L1') ! L1 norm
            Do i=1, imax
                Do j=1 , jmax
                    exact_value = Uexact(case, xm(i), ym(j), time, gamma)
                    error = error + ABS( (U(:,i,j) - exact_value)/exact_value )
                End Do
            End Do
            error = error / ( imax * jmax )


        Case ('L2') ! L2 norm
            Do i=1, imax
                Do j=1 , jmax-1
                    exact_value = Uexact(case, xm(i), ym(j), time, gamma)
                    error = error + ( U(:,i,j) - exact_value )**2
                End Do
            End Do
            error = SQRT( error / ( imax * jmax ) )

            Do k=1,4
                error = error + error(k)
            End Do

            error = (error/4)
        Case Default ! L_infinity norm
            Do i=1, imax
                Do j=1 , jmax
                    exact_value = Uexact(case, xm(i), ym(j), time, gamma)
                    error = MAX( error, ABS( U(:,i,j) - exact_value ) )
                End Do
            End Do

            Do k=1,4
                error = error + error(k)
            End Do

            error = (error/4)
        End Select
    End Function error

    Function Compute_Flux (axis,numflux_name,Uvect,case)

        Character(len=1), Intent(In)                :: axis
        Real(PR), Dimension(4,imax,jmax), Intent(In) :: Uvect
        Character(len=10), Intent(In)                :: numflux_name
        Integer, Intent(In)                          :: case

        Real(PR),Dimension(4,0:imax,0:jmax)          :: Compute_Flux 

       
       

        Real(PR),Dimension(4,0:imax,jmax)          :: U_g_x
        Real(PR),Dimension(4,1:imax+1,jmax)        :: U_d_x
        
        Real(PR),Dimension(4,imax,0:jmax)          :: U_g_y
        Real(PR),Dimension(4,imax,1:jmax+1)        :: U_d_y

        Real(PR),Dimension(4,jmax)               :: U_bounds_x_1,U_bounds_x_d1
        Real(PR),Dimension(4,jmax)               :: U_bounds_x_2,U_bounds_x_d2


        Real(PR),Dimension(4,imax)               :: U_bounds_y_1,U_bounds_y_d1
        Real(PR),Dimension(4,imax)               :: U_bounds_y_2,U_bounds_y_d2



        Integer                                      :: i,j
        
        If (axis=='x') Then
        
   
        
        Do j=1, jmax
            U_g_x(:,0,j)      = Uvect(:,imax,j)
            U_d_x(:,imax+1,j) = Uvect(:,1,j)
            Do i=1, imax-1
                Select Case (TRIM(ADJUSTL(numflux_name)))
                Case ('Rusanov')
                    Compute_Flux(:,i,j) = Rusanov('x', Uvect(:,i,j), Uvect(:,i+1,j), gamma)
                Case('WENO3')

                    U_g_x(:,i,j)    = Uvect(:,i,j)
                    U_d_x(:,i,j)    = Uvect(:,i,j)
                    
                    U_minus_x(:,i,j) = Reconstruction_L(U_g_x(:,i-1,j),Uvect(:,i,j),Uvect(:,i+1,j))
                    U_plus_x(:,i,j)  = Reconstruction_R(Uvect(:,i,j),Uvect(:,i+1,j),U_d_x(:,i+2,j))

                    Compute_Flux(:,i,j) = Rusanov('x',U_minus_x(:,i,j),U_plus_x(:,i,j),gamma)
                Case('Gudonov') 
                    Compute_Flux(:,i,j) = godunov_flux('x', Uvect(:,i,j), Uvect(:,i+1,j), gamma)

                Case Default ! Case ('HLL')
                    Compute_Flux(:,i,j) = HLL('x', Uvect(:,i,j), Uvect(:,i+1,j), gamma)
                End Select
            End Do
            ! Boundary
            Select Case (case)
            Case (2) ! Periodic
                ! Periodic
                Select Case (TRIM(ADJUSTL(numflux_name)))
                Case ('Rusanov')    
                    
                    Compute_Flux(:,0,j) = Rusanov('x', Uvect(:,imax,j), Uvect(:,1,j), gamma)
                    Compute_Flux(:,imax,j) = Compute_Flux(:,0,j)
                Case('WENO3')

                    U_bounds_x_1(:,j) = Reconstruction_L(Uvect(:,imax-1,j),Uvect(:,imax,j),Uvect(:,1,j))
                    U_bounds_x_2(:,j) = Reconstruction_R(Uvect(:,imax,j),Uvect(:,1,j),Uvect(:,2,j))
                    Compute_Flux(:,0,j) = Rusanov('x', U_bounds_x_1(:,j),U_bounds_x_2(:,j), gamma)

                    Compute_Flux(:,imax,j) = Compute_Flux(:,0,j)



                    ! U_bounds_x_d1(:,j)       = Reconstruction_L(Uvect(:,imax-1,j),Uvect(:,imax,j),Uvect(:,imax,j))
                    ! U_bounds_x_d2(:,j)       = Reconstruction_R(Uvect(:,imax,j),Uvect(:,imax,j),Uvect(:,imax,j))
                    ! Compute_Flux(:,imax,j) = Rusanov('x', U_bounds_x_d1(:,j), U_bounds_x_d2(:,j), gamma)

                    ! Compute_Flux(:,0,j) = Rusanov('x', Uvect(:,imax,j), Uvect(:,1,j), gamma)
                    ! Compute_Flux(:,imax,j) = Compute_Flux(:,0,j)




                Case('Gudonov')
                    Compute_Flux(:,0,j) = godunov_flux('x', Uvect(:,imax,j), Uvect(:,1,j), gamma)
                Case Default ! Case ('HLL')
                    Compute_Flux(:,0,j) = HLL('x', Uvect(:,imax,j), Uvect(:,1,j), gamma)
                End Select
                ! Compute_Flux(:,imax,j) = Compute_Flux(:,0,j)
            Case Default ! Absorbing
                Compute_Flux(:,0,j)    = fluxF( Uvect(:,1,j), gamma )
                Compute_Flux(:,imax,j) = fluxF( Uvect(:,imax,j), gamma )
            End Select
        End Do
    
    Else


        Do i=1, imax
            U_g_y(:,i,0)      = Uvect(:,i,jmax)
            U_d_y(:,i,jmax+1) = Uvect(:,i,1)
            
            Do j=1, jmax-1
                Select Case (TRIM(ADJUSTL(numflux_name)))
                Case ('Rusanov')
                    Compute_Flux(:,i,j) = Rusanov('y', Uvect(:,i,j), Uvect(:,i,j+1), gamma)
                Case ('WENO3')
                    U_g_y(:,i,j)    = Uvect(:,i,j)
                    U_d_y(:,i,j)    = Uvect(:,i,j)

                    U_minus_y(:,i,j) = Reconstruction_L(U_g_y(:,i,j-1),Uvect(:,i,j),Uvect(:,i,j+1))
                    U_plus_y(:,i,j)  = Reconstruction_R(Uvect(:,i,j),Uvect(:,i,j+1),U_d_y(:,i,j+2))

                    Compute_Flux(:,i,j) = Rusanov('y',U_minus_y(:,i,j),U_plus_y(:,i,j),gamma)
                Case('Gudonov')
                    Compute_Flux(:,i,j) = godunov_flux('y', Uvect(:,i,j), Uvect(:,i,j+1), gamma)

                Case Default ! Case ('HLL')
                    Compute_Flux(:,i,j) = HLL('y', Uvect(:,i,j), Uvect(:,i,j+1), gamma)
                End Select
            End Do
            ! Boundary
            Select Case (case)
            Case (2) ! Periodic
                Select Case (TRIM(ADJUSTL(numflux_name)))
                Case ('Rusanov')
                    Compute_Flux(:,i,0) = Rusanov('y', Uvect(:,i,jmax), Uvect(:,i,1), gamma)
                    Compute_Flux(:,i,jmax) = Compute_Flux(:,i,0)
                Case('WENO3')

                    U_bounds_y_1(:,i) = Reconstruction_L(Uvect(:,i,jmax-1),Uvect(:,i,jmax),Uvect(:,i,1))
                    U_bounds_y_2(:,i) = Reconstruction_R(Uvect(:,i,jmax),Uvect(:,i,1),Uvect(:,i,2))
                    Compute_Flux(:,i,0) = Rusanov('y', U_bounds_y_1(:,i),U_bounds_y_2(:,i), gamma)

                    Compute_Flux(:,i,jmax) = Compute_Flux(:,i,0)

                    ! U_bounds_y_d1(:,i) = Reconstruction_L(Uvect(:,i,jmax-1),Uvect(:,i,jmax),Uvect(:,i,jmax))
                    ! U_bounds_y_d2(:,i)= Reconstruction_R(Uvect(:,i,jmax),Uvect(:,i,jmax),Uvect(:,i,jmax))
                    ! Compute_Flux(:,i,jmax) = Rusanov('y', U_bounds_y_d1(:,i), U_bounds_y_d2(:,i), gamma)
                    
                    
                    !  Compute_Flux(:,i,0) = Rusanov('y', Uvect(:,i,jmax), Uvect(:,i,1), gamma)
                    !  Compute_Flux(:,i,jmax) = Compute_Flux(:,i,0)
                    
    
                Case('Gudonov')
                    Compute_Flux(:,i,0) = godunov_flux('y', Uvect(:,i,jmax), Uvect(:,i,1), gamma)
                Case Default ! Case ('HLL')
                    Compute_Flux(:,i,0) = HLL('y', Uvect(:,i,jmax), Uvect(:,i,1), gamma)
                End Select
                ! Compute_Flux(:,i,jmax) = Compute_Flux(:,i,0)
            Case Default ! Absorbing
                Compute_Flux(:,i,0)    = fluxG( Uvect(:,i,1), gamma )
                Compute_Flux(:,i,jmax) = fluxG( Uvect(:,i,jmax), gamma )
            End Select
        End Do


    End If


    End Function Compute_Flux


    ! Function Compute_Flux (axis,numflux_name,Uvect,case)

    !     Character(len=1), Intent(In)                :: axis
    !     Real(PR), Dimension(4,imax,jmax), Intent(In) :: Uvect
    !     Character(len=10), Intent(In)                :: numflux_name
    !     Integer, Intent(In)                          :: case

    !     Real(PR),Dimension(4,0:imax,0:jmax)          :: Compute_Flux 

    !     Integer                                      :: i,j
        
    !     If (axis=='x') Then
        
    !     Do j=1, jmax
    !         Do i=1, imax-1
    !             Select Case (TRIM(ADJUSTL(numflux_name)))
    !             Case ('Rusanov')
    !                 Compute_Flux(:,i,j) = Rusanov('x', Uvect(:,i,j), Uvect(:,i+1,j), gamma)
    !             Case('WENO3')

    !                 U_g_x(:,i,j)    = Uvect(:,i,j)
    !                 U_d_x(:,i,j)    = Uvect(:,i,j)
                    
    !                 U_minus_x(:,i,j) = Reconstruction_L(U_g_x(:,i-1,j),Uvect(:,i,j),Uvect(:,i+1,j))
    !                 U_plus_x(:,i,j)  = Reconstruction_R(Uvect(:,i,j),Uvect(:,i+1,j),U_d_x(:,i+2,j))

    !                 Compute_Flux(:,i,j) = Rusanov('x',U_minus_x(:,i,j),U_plus_x(:,i,j),gamma)
    !             Case('Gudonov') 
    !                 Compute_Flux(:,i,j) = godunov_flux('x', Uvect(:,i,j), Uvect(:,i+1,j), gamma)

    !             Case Default ! Case ('HLL')
    !                 Compute_Flux(:,i,j) = HLL('x', Uvect(:,i,j), Uvect(:,i+1,j), gamma)
    !             End Select
    !         End Do
    !         ! Boundary
    !         Select Case (case)
    !         Case (2) ! Periodic
    !             ! Periodic
    !             Select Case (TRIM(ADJUSTL(numflux_name)))
    !             Case ('Rusanov')    
                    
    !                 Compute_Flux(:,0,j) = Rusanov('x', Uvect(:,imax,j), Uvect(:,1,j), gamma)
    !             Case('WENO3')

    !                 ! U_bounds_i1(:,1,j) = Reconstruction_L(Uvect(:,imax,j),Uvect(:,1,j),Uvect(:,2,j))
    !                 ! U_bounds_i2(:,1,j) = Reconstruction_R(Uvect(:,1,j),Uvect(:,2,j),Uvect(:,3,j))
    !                 Compute_Flux(:,0,j) = HLL('x', Uvect(:,imax,j), Uvect(:,1,j), gamma)

    !                 ! U_bounds_k1(:,imax,j) = Reconstruction_L(Uvect(:,imax-1,j),Uvect(:,imax,j),Uvect(:,1,j))
    !                 ! U_bounds_k2(:,imax,j) = Reconstruction_R(Uvect(:,imax,j),Uvect(:,1,j),Uvect(:,2,j))
    !                 ! Compute_Flux(:,imax,j) = Rusanov('x', U_bounds_k1(:,imax,j), U_bounds_k2(:,imax,j), gamma)


    !             Case('Gudonov')
    !                 Compute_Flux(:,0,j) = godunov_flux('x', Uvect(:,imax,j), Uvect(:,1,j), gamma)
    !             Case Default ! Case ('HLL')
    !                 Compute_Flux(:,0,j) = HLL('x', Uvect(:,imax,j), Uvect(:,1,j), gamma)
    !             End Select
    !             Compute_Flux(:,imax,j) = Compute_Flux(:,0,j)
    !         Case Default ! Absorbing
    !             Compute_Flux(:,0,j)    = fluxF( Uvect(:,1,j), gamma )
    !             Compute_Flux(:,imax,j) = fluxF( Uvect(:,imax,j), gamma )
    !         End Select
    !     End Do
    
    ! Else

    !     Do i=1, imax
    !         Do j=1, jmax-1
    !             Select Case (TRIM(ADJUSTL(numflux_name)))
    !             Case ('Rusanov')
    !                 Compute_Flux(:,i,j) = Rusanov('y', Uvect(:,i,j), Uvect(:,i,j+1), gamma)
    !             Case ('WENO3')
    !                 U_g_y(:,i,j)    = Uvect(:,i,j)
    !                 U_d_y(:,i,j)    = Uvect(:,i,j)

    !                 U_minus_y(:,i,j) = Reconstruction_L(U_g_y(:,i,j-1),Uvect(:,i,j),Uvect(:,i,j+1))
    !                 U_plus_y(:,i,j)  = Reconstruction_R(Uvect(:,i,j),Uvect(:,i,j+1),U_d_y(:,i,j+2))

    !                 Compute_Flux(:,i,j) = Rusanov('y',U_minus_y(:,i,j),U_plus_y(:,i,j),gamma)
    !             Case('Gudonov')
    !                 Compute_Flux(:,i,j) = godunov_flux('y', Uvect(:,i,j), Uvect(:,i,j+1), gamma)

    !             Case Default ! Case ('HLL')
    !                 Compute_Flux(:,i,j) = HLL('y', Uvect(:,i,j), Uvect(:,i,j+1), gamma)
    !             End Select
    !         End Do
    !         ! Boundary
    !         Select Case (case)
    !         Case (2) ! Periodic
    !             Select Case (TRIM(ADJUSTL(numflux_name)))
    !             Case ('Rusanov')
    !                 Compute_Flux(:,i,0) = Rusanov('y', Uvect(:,i,jmax), Uvect(:,i,1), gamma)
    !             Case('WENO3')

    !                 ! U_bounds_j1(:,i,1) = Reconstruction_L(Uvect(:,i,jmax),Uvect(:,i,1),Uvect(:,i,2))
    !                 ! U_bounds_j2(:,i,1) = Reconstruction_R(Uvect(:,i,1),Uvect(:,i,2),Uvect(:,i,3))
    !                 Compute_Flux(:,i,0) = Rusanov('y', Uvect(:,i,jmax), Uvect(:,i,1), gamma)

    !                 ! U_bounds_l1(:,i,jmax) = Reconstruction_L(Uvect(:,i,jmax-1),Uvect(:,i,jmax),Uvect(:,i,1))
    !                 ! U_bounds_l2(:,i,jmax) = Reconstruction_R(Uvect(:,i,jmax),Uvect(:,i,1),Uvect(:,i,2))
    !                 ! Compute_Flux(:,i,jmax) = Rusanov('y', U_bounds_l1(:,i,jmax), U_bounds_l2(:,i,jmax), gamma)
                    
                    
    
    !             Case('Gudonov')
    !                 Compute_Flux(:,i,0) = godunov_flux('y', Uvect(:,i,jmax), Uvect(:,i,1), gamma)
    !             Case Default ! Case ('HLL')
    !                 Compute_Flux(:,i,0) = HLL('y', Uvect(:,i,jmax), Uvect(:,i,1), gamma)
    !             End Select
    !             Compute_Flux(:,i,jmax) = Compute_Flux(:,i,0)
    !         Case Default ! Absorbing
    !             Compute_Flux(:,i,0)    = fluxG( Uvect(:,i,1), gamma )
    !             Compute_Flux(:,i,jmax) = fluxG( Uvect(:,i,jmax), gamma )
    !         End Select
    !     End Do


    ! End If


    ! End Function Compute_Flux

    ! Subroutine Time_Scheme(Scheme_name)

    !     Character(len=*), Intent(In) :: Scheme_name


    ! End Subroutine Time_Scheme

    ! Subroutine ExplicitEuler(Ures, Uvect, fluxF, fluxG, deltax, deltay, deltat, imax, jmax)
    !     ! --- InOut
    !     Real(PR), Dimension(:,:,:), Intent(In) :: Uvect
    !     Real(PR), Intent(In) :: deltax, deltay, deltat
    !     Integer, Intent(In) :: imax, jmax
    !     Real(PR), Dimension(:,0:,0:), Intent(In) :: fluxF, fluxG
    !     Real(PR), Dimension(:,:,:), Intent(InOut) :: Ures
    !     ! --- Locals
    !     Integer :: i, j

    !     Do i=1, imax
    !         Do j=1, jmax
    !             Ures(:,i,j) = Uvect(:,i,j) &
    !                 & - deltat/deltax * (fluxF(:,i,j) - fluxF(:,i-1,j)) &
    !                 & - deltat/deltay * (fluxG(:,i,j) - fluxG(:,i,j-1))
    !         End Do
    !     End Do
    ! End Subroutine ExplicitEuler

End Program euler
