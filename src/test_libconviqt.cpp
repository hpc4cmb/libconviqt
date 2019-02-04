// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi.h"

#include <iostream>
#include <iomanip>
#include <exception>
#include <sstream>

#include "conviqt.hpp"

using namespace conviqt;

int main(int argc, char **argv) {

    MPI_Comm comm = MPI_COMM_WORLD;

    std::cout << "Checking that MPI is not yet initialized." << std::endl;

    int flag;
    MPI_Initialized(&flag);
    if (flag) throw std::runtime_error("ERROR: MPI was already initialized");

    std::cout << "Initializing MPI" << std::endl;

    if (MPI_Init(&argc, &argv)) throw std::runtime_error("ERROR: Failed to initialize MPI");

    int ntasks=0, rank=0;
    if (MPI_Comm_size(comm, &ntasks)) throw std::runtime_error("ERROR: Failed get MPI communicator size");
    if (MPI_Comm_rank(comm, &rank)) throw std::runtime_error("ERROR: Failed to get MPI rank");

    std::cout << std::setprecision(16);
    int verbosity=1;

    beam b;
    sky s;
    detector d("LFITEST");
    pointing pnt;

    d.set_epsilon(1.32495160e-04);

    for (int ipol=0; ipol<2; ++ipol) {

        std::cout << "Generating inputs, ipol = " << ipol << std::endl;

        long lmax=32; // default 5000
        long beamlmax=lmax;
        long beammmax=4; // default 14;
        bool pol=(ipol==0); // default false
        double fwhm=4.0;
        std::string beamfile("../data/mb_lfi_30_27_x_rescaled.alm");
        std::string skyfile("../data/slm.fits");

        b.read(beamlmax, beammmax, pol, beamfile, comm);
        s.read(lmax, pol, skyfile, fwhm, comm);
        s.remove_monopole();
        s.remove_dipole();

        // Populate the pointing array

        long ntheta = 3;
        long nphi = 3;
        long npsi = 3;
        long nsamp = ntheta*nphi*npsi;

        pnt.alloc(5 * nsamp);

        long row=0;
        for (long itheta=0; itheta < ntheta; ++itheta) {
            double theta = itheta * (pi / (double)ntheta);
            for (long iphi=0; iphi < nphi; ++iphi) {
                double phi = iphi * (twopi / (double)nphi);
                for (long ipsi=0; ipsi < npsi; ++ipsi) {
                    double psi = ipsi * (pi / (double)npsi);
                    pnt[row*5+0] = phi; // longitude
                    pnt[row*5+1] = theta; // latitude
                    pnt[row*5+2] = psi; // position angle
                    pnt[row*5+3] = 0; // TOD
                    pnt[row*5+4] = row; // time
                    ++row;
                }
            }
        }

        long order=3; // 5

        std::cout << "Instantiating convolver." << std::endl;

        convolver cnv(&s, &b, &d, pol, lmax, beammmax, order, verbosity, comm);

        std::cout << "Convolving." << std::endl;

        cnv.convolve(pnt);

        if (rank == 0) {

            std::cout << "Convolved TOD:" << std::endl;
            for (long row=0; row < nsamp; ++row) {
                std::cout << pnt[row*5+4] << " " << pnt[row*5+3] << std::endl;
            }

            if (pol) {
                if (fabs(pnt[ 0*5+3] -  0.8546349819096275) > 1e-6) {
                    std::ostringstream o;
                    o << "Row 0 should be 0.8546349819096275, not " << pnt[0*5+3];
                    throw std::runtime_error(o.str());
                }
                if (fabs(pnt[10*5+3] + 25.53734467183137  ) > 1e-4) {
                    std::ostringstream o;
                    o << "Row 10 should be -25.53734467183137, not " << pnt[10*5+3];
                    throw std::runtime_error(o.str());
                }
                if (fabs(pnt[15*5+3] + 76.04945574990082  ) > 1e-4) {
                    std::ostringstream o;
                    o << "Row 15 should be -76.04945574990082, not " << pnt[15*5+3];
                    throw std::runtime_error(o.str());
                }
            } else {
                if (fabs(pnt[ 0*5+3] -  0.8545846415739397) > 1e-6) {
                    std::ostringstream o;
                    o << "Row 0 should be 0.8545846415739397, not " << pnt[0*5+3];
                    throw std::runtime_error(o.str());
                }
                if (fabs(pnt[10*5+3] + 25.20150061107036  ) > 1e-4) {
                    std::ostringstream o;
                    o << "Row 10 should be -25.20150061107036, not " << pnt[10*5+3];
                    throw std::runtime_error(o.str());
                }
                if (fabs(pnt[15*5+3] + 76.14723911261254  ) > 1e-4) {
                    std::ostringstream o;
                    o << "Row 15 should be -76.14723911261254, not " << pnt[15*5+3];
                    throw std::runtime_error(o.str());
                }
            }

            std::cout << "Test passed." << std::endl;

        }

    } // ipol loop

    if (MPI_Finalize()) throw std::runtime_error("ERROR: Failed to finalize MPI");

    return 0;
}
