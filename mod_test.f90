Module mod_test

    Use mod_parametres
    Use mod_schemas
    Implicit None

Contains
    Subroutine PerformTest(gamma)
        Real(PR), Intent(In) :: gamma
      
        Write(*,*) '-- Consistence des flux --'
        Write(*,*) "- Rusanov"
        Call testFluxnum(Rusanov, gamma, fluxF, fluxG)
        Write(*,*) "- Reussi" 

        Write(*,*) "- HLL"
        Call testFluxnum(HLL, gamma, fluxF, fluxG)
        Write(*,*) " - Reussi"

        Call Exit()
    
    End Subroutine PerformTest

    Function random_U(gamma)
        Real(PR), Parameter :: eps = EPSILON(Real(PR))
        Real(PR), Parameter :: hug = 1._PR/eps
        Real(PR), Dimension(4), Parameter :: lowerBounds = &
            & (/ Real(PR) :: eps, -hug, -hug, eps /)
        Real(PR), Dimension(4), Parameter :: upperBounds = &
            & (/ Real(PR) :: hug, hug, hug, hug /)

        Real(PR), Intent(In) :: gamma
        Real(PR), Dimension(4) :: random_U

        Real(PR), Dimension(4) :: randomUniform

        Real(PR) :: r, u, v, q, p, e

        Call RANDOM_NUMBER(randomUniform)
        random_U = upperBounds * randomUniform + lowerBounds * ( 1._PR - randomUniform )
        r = random_U(1)
        u = random_U(2)/r
        v = random_U(3)/r
        p = random_U(4)
        q = .5_PR * ( u**2 + v**2 )
        e = p / (gamma - 1._PR) + r*q
        
        random_U(4) = e
    End Function random_U

    Subroutine testFluxnum(numericalFlux, gamma, fluxF, fluxG)
        Real(PR), Parameter :: tolerance = 1._PR
        Integer, Parameter  :: number_of_tests = 10000
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
            U = random_U(gamma)
            failF = SUM( ABS(numericalFlux('x', U, U, gamma) - fluxF(U, gamma)) ) > tolerance*EPSILON(U(1))
            failG = SUM( ABS(numericalFlux('y', U, U, gamma) - fluxG(U, gamma)) ) > tolerance*EPSILON(U(1))
            If (failF .OR. failG) Then
                Call Exit(1)
            End If
        End Do
    End Subroutine testFluxnum

End Module mod_test
