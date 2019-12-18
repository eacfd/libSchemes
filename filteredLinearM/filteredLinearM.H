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
    Foam::filteredLinearMLimiter

Description
    Class to generate weighting factors for the filteredLinearM
    differencing scheme.

    The aim is to remove high-frequency modes with "staggering"
    characteristics by comparing the face gradient relative to the
    background distribution represented by the neighbouring cell gradients
    with those gradients and introduce small amounts of upwind in order to
    damp these modes.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    filteredLinearM.C

\*---------------------------------------------------------------------------*/

#ifndef filteredLinearM_H
#define filteredLinearM_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class filteredLinearMWeight Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class filteredLinearMLimiter
:
    public LimiterFunc
{
    scalar k_;
    scalar l_;

public:

    filteredLinearMLimiter(Istream& is)
    :
        k_(readScalar(is)),
        l_(readScalar(is))
    {}

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
        scalar df = phiN - phiP;

        scalar dcP = d & gradcP;
        scalar dcN = d & gradcN;

        scalar limiter =
            k_
          - 0.5*min(mag(df - dcP), mag(df - dcN))
           /(max(mag(dcP), mag(dcN)) + SMALL);

        return max(min(limiter, 1), l_);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
