module mod_rdm2ll

    implicit none
    private

    public :: RDM2LL

contains

    subroutine RDM2LL(y,x,lat,lon)

        implicit none
        real x,y,lat,lon
        real referenceRdX,referenceRdY,referenceWgs84X,referenceWgs84Y
        real dX,dY
        real sumN,sumE

        !SOURCE:   https://www.roelvanlisdonk.nl/2012/11/21/simple-way-for-converting-rijksdriehoek-coordinates-to-lat-and-long-wgs84-in-c/
        !referentie coordinated RDM

        referenceRdX = 155000
        referenceRdY = 463000

        !The city "Amsterfoort" is used as reference "WGS84" coordinate.
        referenceWgs84X = 52.15517
        referenceWgs84Y = 5.387206

        dX = (x - referenceRdX) * 10**(-5.0)
        dY = (y - referenceRdY) * 10**(-5.0)
        sumN = (3235.65389 * dY) + (-32.58297 * (dX**2)) + (-0.2475 * (dY**2)) + (-0.84978 * (dX**2) * dY) + (-0.0655 * (dY**3)) + (-0.01709 * (dX**2) * (dY**2)) + (-0.00738 * dX) + (0.0053 * (dX**4)) + (-0.00039 * (dX**2) * (dY**3)) + (0.00033 * (dX**4) * dY) + (-0.00012 * dX * dY)

        sumE = (5260.52916 * dX) + (105.94684 * dX * dY) + (2.45656 * dX * (dY**2)) + (-0.81885 * (dX**3)) + (0.05594 * dX * (dY**3)) + (-0.05607 * (dX**3) * dY) + (0.01199 * dY) + (-0.00256 *(dX**3) * (dY**2)) + (0.00128 * dX * (dY**4)) + (0.00022 * (dY**2)) + (-0.00022 *(dX**2)) + (0.00026 * (dX**5))

        lat = referenceWgs84X + (sumN / 3600.)
        lon = referenceWgs84Y + (sumE / 3600.)

    end subroutine RDM2LL

end module mod_rdm2ll

