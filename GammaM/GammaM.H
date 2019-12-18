/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::GammaMLimiter

Description
    Class with limiter function which returns the limiter for the
    GammaM differencing scheme based on phict obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    GammaM.C

\*---------------------------------------------------------------------------*/

#ifndef GammaM_H
#define GammaM_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class GammaMLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class GammaMLimiter
:
    public LimiterFunc
{
    scalar k_;
    scalar l_;

public:

    GammaMLimiter(Istream& is)
    :
        k_(readScalar(is)),
        l_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorIn("GammaMLimiter(Istream& is)", is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        k_ = max(k_, SMALL);
    }

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        scalar phict = LimiterFunc::phict
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );

        if ((phict <= 0)||(phict >= 1))
        {
            if (phict >= 1) 
            {
               return 1;
            }
            if (phict <= 0)
            {
               return l_;
            }
        }

        if ((phict < 1)&&(phict >= k_))
        {
            return 1;
        }
        if ((phict < k_)&&(phict > 0))
        {
            return  min(max((phict/k_), 0), 1);
        }
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
