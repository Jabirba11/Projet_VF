Module mod_fonctions

    Use mod_parametres
    Implicit None

Contains
   
    Function Uinit(x, y, gamma, case)
        
        Real(PR), Intent(In)   :: x, y, gamma
        Integer, Intent(In)    :: case
        Real(PR), Dimension(4) :: Uinit
        
        Real(PR) :: rho, u, v, q, p, e


        Select Case (case)
        
        Case (1) ! Tube a choc
            If (x < 0.5_PR) Then
                ! Gauche
                rho = 1._PR
                u = 0._PR
                v = 0._PR
                p = 1._PR
            Else
                ! Droit
                rho = 0.125_PR
                u = 0._PR
                v = 0._PR
                p = 0.1_PR
            End If
        
        Case (2) !  vortex isentropique
            Call IsentropicVortexCalculation(x, y, 0._PR, &
                & 45._PR, &
                & SQRT(2._PR/gamma), &
                & 1._PR, 1._PR, &
                & 1._PR, 5._PR, 1._PR, &
                & SQRT(EXP(1._PR)/gamma)*5._PR/(2._PR*PI), &
                & gamma, &
                & rho, u, v, p)
            
        Case(3)
            Call LiskaWendroff(x,y,rho,u,v,p,0._PR)
        
            Case Default
            rho = 1._PR
            u = 0._PR
            v = 0._PR
            p = 1._PR
        
        End Select

        q = 0.5_PR * ( u**2 + v**2 )
        e = p / (gamma - 1._PR) + rho*q

        Uinit = (/ Real(PR) :: rho, rho*u, rho*v, e /)
    
    End Function Uinit

    Function Uexact(case, x, y, t, gamma)
        
        Real(PR), Intent(In)   :: x, y, t, gamma
        Integer, Intent(In)    :: case
        Real(PR), Dimension(4) :: Uexact
        Real(PR) :: rho, u, v, q, p, e

        Select Case (case)
        Case (2) ! vortex isentropique
            Call IsentropicVortexCalculation(x, y, t, &
                & 45._PR, &
                & SQRT(2._PR/gamma), &
                & 1._PR, 1._PR, &
                & 1._PR, 5._PR, 1._PR, &
                & SQRT(EXP(1._PR)/gamma)*5._PR/(2._PR*PI), &
                & gamma, &
                & rho, u, v, p)
        Case(3) ! Liska Wendroff p9
            Call LiskaWendroff(x,y,rho,u,v,p,t)
        Case Default
            Write(STDERR, *) "Pas de solution trouvée dans ce cas (case ", case, ")"
            Call Exit(1)
        End Select

        q = .5_PR * ( u**2 + v**2 )
        e = p / (gamma - 1._PR) + rho*q

        Uexact = (/ Real(PR) :: rho, rho*u, rho*v, e /)
   
    End Function Uexact

    Subroutine IsentropicVortexCalculation(posX, posY, currentTime, angleDegrees, infinityMachNumber, &
        & infinityDensity, infinityPressure, radiusOfVortex, halfLengthOfDomain, standardDeviation, &
        & strengthOfPerturbation, specificHeatRatio, &
        & density, velocityX, velocityY, pressure)

        Real(PR), Intent(In) :: posX, posY, currentTime, angleDegrees, infinityMachNumber, &
            & infinityDensity, infinityPressure, radiusOfVortex, halfLengthOfDomain, &
            & standardDeviation, strengthOfPerturbation, specificHeatRatio
        Real(PR), Intent(Out) :: density, velocityX, velocityY, pressure

        Real(PR) :: radialFunction, vortexDisturbance, speedOfSoundAtInfinity
        Real(PR) :: freeStreamVelocityX, freeStreamVelocityY
        Real(PR) :: transformedX, transformedY
        Real(PR) :: angleRadians

        speedOfSoundAtInfinity = SQRT( specificHeatRatio * infinityPressure/infinityDensity )
        angleRadians = angleDegrees * PI / 180._PR
        freeStreamVelocityX = speedOfSoundAtInfinity * infinityMachNumber * COS( angleRadians )
        freeStreamVelocityY = speedOfSoundAtInfinity * infinityMachNumber * SIN( angleRadians )

        transformedX = MODULO( posX - freeStreamVelocityX*currentTime + halfLengthOfDomain, &
            & 2*halfLengthOfDomain ) - halfLengthOfDomain
        transformedY = MODULO( posY - freeStreamVelocityY*currentTime + halfLengthOfDomain, &
            & 2*halfLengthOfDomain ) - halfLengthOfDomain
        radialFunction = ( transformedX/radiusOfVortex )**2 + ( transformedY/radiusOfVortex )**2
        radialFunction = - radialFunction / ( 2._PR * standardDeviation**2 )
        vortexDisturbance = strengthOfPerturbation * EXP( radialFunction )

        density = infinityDensity * ( 1._PR  -  .5_PR*(specificHeatRatio - 1._PR)*vortexDisturbance**2 )** &
                (1._PR / (specificHeatRatio - 1._PR)) 
        velocityX = freeStreamVelocityX - posY/radiusOfVortex*vortexDisturbance
        velocityY = freeStreamVelocityY + posX/radiusOfVortex*vortexDisturbance
        pressure = infinityPressure / specificHeatRatio * ( 1._PR -.5_PR * ( specificHeatRatio - 1._PR ) * &
                vortexDisturbance**2 )**(specificHeatRatio / (specificHeatRatio - 1._PR))

    End Subroutine IsentropicVortexCalculation

    Subroutine LiskaWendroff(x,y,density,velocityX,velocityY,pressure,Time)

        Real(PR), Intent(In)  :: x,y,Time

        Real(PR), Intent(Out) :: density,velocityX,velocityY,pressure
        
        density    = 1._PR+0.2_PR*sin(PI*(x+y-Time))
        velocityX  = 1._PR 
        velocityY  = -0.5_PR
        pressure   = 1._PR
    
    End Subroutine LiskaWendroff



End Module mod_fonctions