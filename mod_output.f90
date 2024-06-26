Module mod_output
    
    Use mod_parametres

    Implicit None

Contains
    Subroutine output(U, gamma, x, y, nb_filename, filename_prefix)
        

        ! Les variables d'entrée
        Integer, Intent(In)                    :: nb_filename
        Character(len=*)                       :: filename_prefix
        Real(PR), Dimension(0:), Intent(In)    :: x, y
        Real(PR), Dimension(:,:,:), Intent(In) :: U
        Real(PR), Intent(In)                   :: gamma

        ! Les variables locales
        Integer                                  :: imax, jmax, i, j
        Character(len=30)                        :: nomfichier
        Real(PR), Dimension(Size(x)-1,Size(y)-1) :: rho, vitesse_x, vitesse_y, energy, q, pressure

        imax = Size(x) - 1
        jmax = Size(y) - 1

        rho         = U(1,:,:)
        vitesse_x   = U(2,:,:) / rho
        vitesse_y   = U(3,:,:) / rho
        energy      = U(4,:,:)

        q = .5_PR * ( vitesse_x**2 + vitesse_y**2 )
        pressure = (gamma - 1._PR)*(energy - rho*q)

        Write(nomfichier,*) nb_filename
        Open(Unit=111, File=TRIM(ADJUSTL(filename_prefix))//'_'//TRIM(ADJUSTL(nomfichier))//'.vtk')

        Write(111,'(a)') '# vtk DataFile Version 2.0'
        Write(111,'(a)') 'Euler Equations'
        Write(111,'(a)') 'ASCII'
        Write(111,'(a)') 'DATASET RECTILINEAR_GRID'
        Write(111,'(1a10,3i4)') 'DIMENSIONS ',imax+1,jmax+1,1
        Write(111,fmt='(1a13,1i10,1a10)') 'X_COORDINATES',imax+1,' float'
        Do i=0,imax
            Write(111,*) x(i)
        End Do
        Write(111,fmt='(1a13,1i10,1a10)') 'Y_COORDINATES',jmax+1,' float'
        Do j=0,jmax
            Write(111,*) y(j)
        End Do
        Write(111,fmt='(1a13,1i10,1a10)') 'Z_COORDINATES',1,' float'
        Write(111,*) 0.0_PR

        Write(111,fmt='(1a9,1i10)') 'CELL_DATA', imax*jmax
        Write(111,'(a)') 'SCALARS density double'
        Write(111,'(a)') 'LOOKUP_TABLE default'
        Do j=1,jmax
            Do i=1,imax
                Write(111,*) rho(i,j)
            End Do
        End Do
        Write(111,'(a)') 'SCALARS pressure double'
        Write(111,'(a)') 'LOOKUP_TABLE default'
        Do j=1,jmax
            Do i=1,imax
                Write(111,*) pressure(i,j)
            End Do
        End Do
        Write(111,'(a)') 'VECTORS velocity double'
        Do j=1,jmax
            Do i=1,imax
                Write(111,*) vitesse_x(i,j), vitesse_y(i,j), 0._PR
            End do
        End do

        Close(111)

    End Subroutine output


    Subroutine output_csv(U, gamma, x, y, nb_filename, filename_prefix)
        
        
        !        ! Les variables d'entrée
        Integer, Intent(In)                    :: nb_filename
        Character(len=*)                       :: filename_prefix
        Real(PR), Dimension(0:), Intent(In)    :: x, y
        Real(PR), Dimension(:,:,:), Intent(In) :: U
        Real(PR), Intent(In)                   :: gamma

        ! Les variables locales
        Integer                                  :: imax, jmax, i, j
        Character(len=30)                        :: nomfichier
        Real(PR), Dimension(Size(x)-1,Size(y)-1) :: rho, vitesse_x, vitesse_y, energy, q, pressure

        imax = Size(x) - 1
        jmax = Size(y) - 1

        rho         = U(1,:,:)
        vitesse_x   = U(2,:,:) / rho
        vitesse_y   = U(3,:,:) / rho
        energy      = U(4,:,:)

        q = .5_PR * ( vitesse_x**2 + vitesse_y**2 )
        pressure = (gamma - 1._PR)*(energy - rho*q)
    
        ! Nom de fichier pour CSV ou texte simple
        Write(nomfichier,*) nb_filename
        Open(Unit=111, File=TRIM(ADJUSTL(filename_prefix))//'_'//TRIM(ADJUSTL(nomfichier))//'.csv')
    
        ! Écrire l'en-tête pour CSV
        Write(111,*) 'x, y, density, pressure, velocity_x, velocity_y'
    
        ! Écrire les données
        Do j=1,jmax
            Do i=1,imax
                Write(111, '(F8.3, ", ", F8.3, ", ", F8.3, ", ", F8.3, ", ", F8.3, ", ", F8.3)') &
                    x(i), y(j), rho(i,j), pressure(i,j), vitesse_x(i,j), vitesse_y(i,j)
            End Do
        End Do
    
        Close(111)
    
    End Subroutine output_csv
    

End Module mod_output
