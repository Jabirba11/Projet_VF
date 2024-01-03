Module mod_schemas

   Use mod_parametres
   Implicit None

Contains
    
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
            
        End Do

        alpha0 = d0/((beta0+epsilon)**2)
        alpha1 = d1/((beta1+epsilon)**2)

        W0     = alpha0/(alpha0+alpha1)
        W1     = alpha1/(alpha0+alpha1)

        Do j=1,4

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
                
            End Do
    
            alpha0 = d0/((beta0+epsilon)**2)
            alpha1 = d1/((beta1+epsilon)**2)
    
            W0     = alpha0/(alpha0+alpha1)
            W1     = alpha1/(alpha0+alpha1)
    
            Do j=1,4
    
                U_plus(j) = W0*(0.5_PR*Ud(j) + 0.5_PR*U(j) ) + W1*(-0.5_PR*Udd(j)+(3._PR/2._PR)*Ud(j))
            
            End Do
        


    End Function Reconstruction_R


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
