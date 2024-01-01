Module mod_test

    Use mod_parametres
    Use mod_schemas
    Implicit None

Contains
    Subroutine PerformTestsAndExit(gamma)
        Real(PR), Intent(In) :: gamma
      
        Write(*,*) '=== Consistence des flux ==='
        Write(*,*) "-- Rusanov"
        Call testNumericalFlux(Rusanov, gamma, fluxF, fluxG)
        Write(*,*) " -> SUCCESS" 

        Write(*,*) "-- HLL"
        Call testNumericalFlux(HLL, gamma, fluxF, fluxG)
        Write(*,*) " -> SUCCESS"

        Call Exit()
    
    End Subroutine PerformTestsAndExit

    Function randomU(gamma)
        Real(PR), Parameter :: eps = EPSILON(Real(PR))
        Real(PR), Parameter :: hug = 1._PR/eps
        Real(PR), Dimension(4), Parameter :: lowerBounds = &
            & (/ Real(PR) :: eps, -hug, -hug, eps /)
        Real(PR), Dimension(4), Parameter :: upperBounds = &
            & (/ Real(PR) :: hug, hug, hug, hug /)

        Real(PR), Intent(In) :: gamma
        Real(PR), Dimension(4) :: randomU

        Real(PR), Dimension(4) :: randomUniform

        Real(PR) :: r, u, v, q, p, e

        Call RANDOM_NUMBER(randomUniform)
        randomU = upperBounds * randomUniform + lowerBounds * ( 1._PR - randomUniform )
        r = randomU(1)
        u = randomU(2)/r
        v = randomU(3)/r
        p = randomU(4)
        q = .5_PR * ( u**2 + v**2 )
        e = p / (gamma - 1._PR) + r*q
        
        randomU(4) = e
    End Function randomU

    Subroutine testNumericalFlux(numericalFlux, gamma, fluxF, fluxG)
        Real(PR), Parameter :: safety_factor = 2._PR
        Integer, Parameter :: number_of_tests = 1000
        ! ------------------- Intent In -----------------------
        Real(PR), Intent(In) :: gamma
        Interface
            Function numericalFlux(axis, UL, UR, gamma)
                Import PR
                Real(PR), Dimension(4), Intent(In) :: UL, UR
                Real(PR), Intent(In) :: gamma
                Character, Intent(In) :: axis
                Real(PR), Dimension(4) :: numericalFlux
            End Function numericalFlux
        End Interface
        Interface
            Function fluxF(Uvect, gamma)
                Import PR
                Real(PR), Dimension(4), Intent(In) :: Uvect
                Real(PR), Intent(In) :: gamma
                Real(PR), Dimension(4) :: fluxF
            End Function fluxF
        End Interface
        Interface
            Function fluxG(Uvect, gamma)
                Import PR
                Real(PR), Dimension(4), Intent(In) :: Uvect
                Real(PR), Intent(In) :: gamma
                Real(PR), Dimension(4) :: fluxG
            End Function fluxG
        End Interface
        ! ----------------------------------------------------------
        Integer :: i
        Real(PR), Dimension(4) :: U
        Logical :: failF, failG

        Do i=1, number_of_tests
            U = randomU(gamma)
            failF = SUM( ABS(numericalFlux('x', U, U, gamma) - fluxF(U, gamma)) ) > safety_factor*EPSILON(U(1))
            failG = SUM( ABS(numericalFlux('y', U, U, gamma) - fluxG(U, gamma)) ) > safety_factor*EPSILON(U(1))
            If (failF .OR. failG) Then
                Call Exit(1)
            End If
        End Do
    End Subroutine testNumericalFlux

End Module mod_test
