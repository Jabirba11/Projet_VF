Module mod_schemas

   Use mod_parametres
   Implicit None

Contains
    
Function godunov_flux(direction, u, v, gamma) result(h)
    Character, Intent(In)              :: direction 
    real(PR), dimension(4), intent(in) :: u, v
    real(PR), intent(in)               :: gamma
    real(PR), dimension(4)             :: h
    real(PR), dimension(4)             :: f_u, f_v
    integer :: i

    If (direction == 'x') Then
        do i = 1, 4
            f_u = fluxF(u, gamma)
            f_v = fluxF(v, gamma)
            if (u(i) <= v(i)) then
                h(i) = min(f_u(i), f_v(i))
            else
                h(i) = max(f_u(i), f_v(i))
            endif
        end do
    Else
        do i = 1, 4
            f_u = fluxG(u, gamma)
            f_v = fluxG(v, gamma)
            if (u(i) <= v(i)) then
                h(i) = min(f_u(i), f_v(i))
            else
                h(i) = max(f_u(i), f_v(i))
            endif
        end do
    End If
End Function godunov_flux
   
    


Function Rusanov(direction, UL, UR, gamma)
        
        ! Variables d'entrée et de sortie
        Real(PR), Dimension(4), Intent(In) :: UL, UR
        Real(PR), Intent(In)               :: gamma
        Character, Intent(In)              :: direction 
        Real(PR), Dimension(4)             :: Rusanov
        
        ! Variables locales
        Real(PR) :: densL, vitesse_xL, vitesse_yL, energyL, pressionL, qL, aL, b, maxWaveL1, maxWaveL3
        Real(PR) :: densR, vitesse_xR, vitesse_yR, energyR, pressionR, qR, aR, maxWaveR1, maxWaveR3

        ! Traitement de l'état gauche
        densL      = UL(1)
        vitesse_xL = UL(2)/densL
        vitesse_yL = UL(3)/densL
        energyL    = UL(4)
        
        qL         = 0.5_PR * ( vitesse_xL**2 + vitesse_yL**2 ) 
        pressionL  = (gamma - 1._PR)*(energyL - densL*qL)
        aL         = SQRT(gamma * pressionL / densL)


        ! Traitement de l'état droit
        densR      = UR(1)
        vitesse_xR = UR(2)/densR
        vitesse_yR = UR(3)/densR
        energyR    = UR(4)
        
        qR         = 0.5_PR * ( vitesse_xR**2 + vitesse_yR**2 )     
        pressionR  = (gamma - 1._PR)*(energyR - densR*qR)
        aR         = SQRT(gamma * pressionR / densR)
 
  
        If (direction == 'y')Then
            maxWaveL1 = ABS(vitesse_yL - aL)
            maxWaveR1 = ABS(vitesse_yR - aR)
            maxWaveL3 = ABS(vitesse_yL + aL)
            maxWaveR3 = ABS(vitesse_yR + aR)
        Else
            maxWaveL1 = ABS(vitesse_xL - aL)
            maxWaveR1 = ABS(vitesse_xR - aR)
            maxWaveL3 = ABS(vitesse_xL + aL)
            maxWaveR3 = ABS(vitesse_xR + aR)
        End If

      
        b = MAX( MAX(maxWaveL1, maxWaveR1),MAX(maxWaveL3, maxWaveR3))

        Select Case (direction)
        Case ('y')
            Rusanov = fluxG(UL, gamma) + fluxG(UR, gamma)
        Case Default 
            Rusanov = fluxF(UL, gamma) + fluxF(UR, gamma)
        End Select
        Rusanov = 0.5_PR * Rusanov -  0.5_PR * b *(UR - UL) 
    
    End Function Rusanov




    Function HLL(direction, UL, UR, gamma)
        
        ! Variables d'entrée et de sortie
        Real(PR), Dimension(4), Intent(In) :: UL, UR
        Real(PR), Intent(In)               :: gamma
        Character, Intent(In)              :: direction 
        Real(PR), Dimension(4)             :: HLL
        ! Variables locales
        Real(PR) :: densL, vitesse_xL, vitesse_yL, energyL, pressionL, qL, aL, b_moins,maxWaveL1, maxWaveL3
        Real(PR) :: densR, vitesse_xR, vitesse_yR, energyR, pressionR, qR, aR,b_plus,maxWaveR1, maxWaveR3

        ! Traitement de l'état gauche
        densL      = UL(1)
        vitesse_xL = UL(2)/densL
        vitesse_yL = UL(3)/densL
        energyL    = UL(4)
        
        
        qL        = 0.5_PR * ( vitesse_xL**2 + vitesse_yL**2 )
        pressionL = (gamma - 1._PR)*(energyL - densL*qL)
        aL        = SQRT(gamma * pressionL / densL)

        
        ! Traitement de l'état droit
        densR      = UR(1)
        vitesse_xR = UR(2)/densR
        vitesse_yR = UR(3)/densR
        energyR    = UR(4)
        
        qR        = 0.5_PR * ( vitesse_xR**2 + vitesse_yR**2 )
        pressionR = (gamma - 1._PR)*(energyR - densR*qR)
        aR        = SQRT(gamma * pressionR / densR)

        
        If (direction=='y') Then
            maxWaveL1 = vitesse_yL - aL
            maxWaveR1 = vitesse_yR - aR
            maxWaveL3 = vitesse_yL + aL
            maxWaveR3 = vitesse_yR + aR
        
        Else
            maxWaveL1 = vitesse_xL - aL
            maxWaveR1 = vitesse_xR - aR
            maxWaveL3 = vitesse_xL + aL
            maxWaveR3 = vitesse_xR + aR
        End If
        


        b_moins = MIN( MIN( maxWaveL1, maxWaveL3 ), MIN( maxWaveR1, maxWaveR3 ), 0._PR )
        b_plus = MAX( MAX( maxWaveL1, maxWaveL3 ), MAX( maxWaveR1, maxWaveR3 ), 0._PR )

        Select Case (direction)
        Case ('y')
            HLL = b_plus * fluxG(UL, gamma) - b_moins * fluxG(UR, gamma)
            HLL = ( HLL   +   b_plus * b_moins * (UR - UL) ) / ( b_plus - b_moins )
        Case Default 
            HLL = b_plus * fluxF(UL, gamma) - b_moins * fluxF(UR, gamma)
            HLL = ( HLL   +   b_plus * b_moins * (UR - UL) ) / ( b_plus - b_moins )
        End Select
    End Function HLL


    Function Reconstruction_L(U_left,U_mid,U_right)result(U_minus)

        
        !Variables d'entrée et de sortie
        real(PR),dimension(4),Intent(In) :: U_left,U_mid,U_right
        real(PR),dimension(4)             :: U_minus

        !Variables locales
        integer  :: j
        real(PR) :: alpha0,alpha1,W0,W1
        real(PR) :: epsilon,d0,d1,beta0,beta1

        epsilon = 10e-6

        d0     = 1._PR/3._PR 
        d1     = 2._PR/3._PR
        

        Do j=1,4

            beta0 = (U_mid(j)-U_left(j))**2
            beta1 = (U_right(j)-U_mid(j))**2
            


            alpha0 = d0/((beta0+epsilon)**2)
            alpha1 = d1/((beta1+epsilon)**2)

            W0     = alpha0/(alpha0+alpha1)
            W1     = alpha1/(alpha0+alpha1)


            U_minus(j) = W0*(-0.5_PR*U_left(j) + (3._PR/2._PR)*U_mid(j) ) + W1*(0.5_PR*U_mid(j)+0.5_PR*U_right(j))
        
        End Do


        End Function Reconstruction_L

        Function Reconstruction_R(U,Ud,Udd)result(U_plus)

        
            !Variables d'entrée et de sortie
            real(PR),dimension(4),Intent(In) :: U,Ud,Udd
            real(PR),dimension(4)            :: U_plus
    
            !Variables locales
            integer  :: j
            real(PR) :: alpha0,alpha1,W0,W1
            real(PR) :: epsilon,d0,d1,beta0,beta1
    
            epsilon = 10e-6
    
            d0     = 2._PR/3._PR 
            d1     = 1._PR/3._PR
            
    
            Do j=1,4
    
                beta0 = (Ud(j)-U(j))**2
                beta1 = (Udd(j)-Ud(j))**2

                alpha0 = d0/((beta0+epsilon)**2)
                alpha1 = d1/((beta1+epsilon)**2)
        
                W0     = alpha0/(alpha0+alpha1)
                W1     = alpha1/(alpha0+alpha1)

                U_plus(j) = W0*(0.5_PR*Ud(j) + 0.5_PR*U(j) ) + W1*(-0.5_PR*Udd(j)+(3._PR/2._PR)*Ud(j))
        
                
            End Do



    End Function Reconstruction_R

    Function WENO5_Left(U_left_1,U_left,U_mid,U_right,U_right_1)

        
        !Variables d'entrée et de sortie
        real(PR),dimension(3),Intent(In) :: U_left,U_mid,U_right,U_left_1,U_right_1
        real(PR),dimension(3)            :: WENO5_Left
        

        !Variables locales
        integer  :: j
        real(PR) :: alpha0,alpha1,W0,W1,W2,alpha2
        real(PR) :: epsilon,d0,d1,d2,beta0,beta1,beta2

        epsilon = 10e-6

        d0     = 1._PR/10._PR 
        d1     = 3._PR/5._PR
        d2     = 3._PR/10._PR
        

        Do j=1,3

            beta0 = (13._PR/12._PR)*((U_left_1(j)-2._PR*U_left(j)+U_mid(j))**2) &
            & + ((U_left_1(j)-4._PR*U_left(j)+3._PR*U_mid(j))**2)/4._PR
            
            beta1 = (13._PR/12._PR)*((U_left(j)-2._PR*U_mid(j)+U_right(j))**2) &
            & + ((U_left(j)-U_right(j))**2)/4._PR 
            
            beta2 = (13._PR/12._PR)*((U_mid(j)-2._PR*U_right(j)+U_right_1(j))**2) &
            & + ((3._PR*U_mid(j)-4._PR*U_right(j)+U_right_1(j))**2)/4._PR

            alpha0 = d0/((beta0+epsilon)**2)
            alpha1 = d1/((beta1+epsilon)**2)
            alpha2 = d2/((beta2+epsilon)**2)

    
            W0     = alpha0/(alpha0+alpha1+alpha2)
            W1     = alpha1/(alpha0+alpha1+alpha2)
            W2     = alpha2/(alpha0+alpha1+alpha2)


            WENO5_Left(j) = W0*(U_left_1(j)/3._PR - (7._PR/6._PR)*U_left(j)+(11._PR/6._PR)*U_mid(j)) +&
            & W1*((-1._PR/6._PR)*U_left(j)+(5._PR/6._PR)*U_mid(j)+(1._PR/3._PR)*U_right(j)) +&
            & W2*((1._PR/3._PR)*U_mid(j)+(5._PR/6._PR)*U_right(j)-(1._PR/6._PR)*U_right_1(j))

            
        End Do


    End Function WENO5_Left

    Function WENO5_Right(U_left_1,U_left,U_mid,U_right,U_right_1)

    
        !Variables d'entrée et de sortie
        real(PR),dimension(3),Intent(In) :: U_left,U_mid,U_right,U_left_1,U_right_1
        real(PR),dimension(3)            :: WENO5_Right
        

        !Variables locales
        integer  :: j
        real(PR) :: alpha0,alpha1,W0,W1,W2,alpha2
        real(PR) :: epsilon,d0,d1,d2,beta0,beta1,beta2

        epsilon = 10e-6  

        d0     = 3._PR/10._PR 
        d1     = 3._PR/5._PR
        d2     = 1._PR/10._PR
        

        Do j=1,3

            beta0 = (13._PR/12._PR)*((U_left_1(j)-2._PR*U_left(j)+U_mid(j))**2) &
            & + ((U_left_1(j)-4._PR*U_left(j)+3._PR*U_mid(j))**2)/4._PR
            
            beta1 = (13._PR/12._PR)*((U_left(j)-2._PR*U_mid(j)+U_right(j))**2) &
            & + ((U_left(j)-U_right(j))**2)/4._PR 
            
            beta2 = (13._PR/12._PR)*((U_mid(j)-2._PR*U_right(j)+U_right_1(j))**2) &
            & + ((3._PR*U_mid(j)-4._PR*U_right(j)+U_right_1(j))**2)/4._PR

            alpha0 = d0/((beta0+epsilon)**2)
            alpha1 = d1/((beta1+epsilon)**2)
            alpha2 = d2/((beta2+epsilon)**2)

    
            W0     = alpha0/(alpha0+alpha1+alpha2)
            W1     = alpha1/(alpha0+alpha1+alpha2)
            W2     = alpha2/(alpha0+alpha1+alpha2)


            WENO5_Right(j) = W0*(U_mid(j)/3._PR + (5._PR/6._PR)*U_left(j)+(-1._PR/6._PR)*U_left_1(j)) +&
            & W1*((-1._PR/6._PR)*U_right(j)+(5._PR/6._PR)*U_mid(j)+(1._PR/3._PR)*U_left(j)) +&
            & W2*((1._PR/3._PR)*U_right_1(j)+(-7._PR/6._PR)*U_right(j)+(11._PR/6._PR)*U_mid(j))

            
        End Do





    End Function WENO5_Right





    ! Fluxes functions
    Function fluxF(Uvect, gamma)
       
       !Variables d'entrée et de sortie
        Real(PR), Dimension(4), Intent(In) :: Uvect
        Real(PR), Intent(In)               :: gamma
        Real(PR), Dimension(4)             :: fluxF
        
        !Variables locales
        Real(PR) :: r, ru, rv, e, u, v, p, q

        r  = Uvect(1)
        ru = Uvect(2)
        rv = Uvect(3)
        e  = Uvect(4)
        u  = ru/r
        v  = rv/r
        q  = 0.5_PR * ( u**2 + v**2 )
        p  = (gamma - 1._PR)*(e - r*q)

        fluxF(1) = ru
        fluxF(2) = ru*u + p
        fluxF(3) = rv*u
        fluxF(4) = (e + p)*u
    End Function fluxF

    Function fluxG(Uvect, gamma)
        
        ! Variables d'entrée et de sortie
        Real(PR), Dimension(4), Intent(In) :: Uvect
        Real(PR), Intent(In)               :: gamma
        Real(PR), Dimension(4)             :: fluxG
        
        ! Variables locales
        Real(PR) :: r, ru, rv, e, u, v, p, q

        r = Uvect(1)
        ru = Uvect(2)
        rv = Uvect(3)
        e = Uvect(4)
        u = ru/r
        v = rv/r
        q = 0.5_PR * ( u**2 + v**2 )
        p = (gamma - 1._PR)*(e - r*q)

        fluxG(1) = rv
        fluxG(2) = ru*v
        fluxG(3) = rv*v + p
        fluxG(4) = (e + p)*v
    End Function fluxG

End Module mod_schemas
