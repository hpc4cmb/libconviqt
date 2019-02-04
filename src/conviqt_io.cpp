#include "conviqt.hpp"

// This file will contain the I/O facilities: read beam, read sky, read detector

// An end user will have the choice between using these facilities or providing inputs directly

using namespace levels;

namespace conviqt {

int beam::read(long beamlmax, long beammmax, bool beampol,
               std::string infile_beam, MPI_Comm comm) {

    MPI_Manager mpiMgr(comm);
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (verbosity > 5) {
        std::cout << "Reading " << infile_beam << " on task = " << rank << std::endl;
    }

    pol = beampol;
    int lmax_file, mmax_file;
    if (rank == 0) {
        if (pol) {
            get_almsize_pol(infile_beam, lmax_file, mmax_file);
        } else {
            get_almsize(infile_beam, lmax_file, mmax_file);
        }
    }
    mpiMgr.bcast(lmax_file);
    mpiMgr.bcast(mmax_file);
    // Ensure that the beam in the file has a high enough expansion order.
    if (lmax_file < beamlmax) {
        throw std::runtime_error("Beam lmax exceeds available expansion order.");
    }
    if (mmax_file < beammmax) {
        throw std::runtime_error("Beam mmax exceeds available expansion order.");
    }
    // lmax
    if (beamlmax > 0) {
        lmax = beamlmax;
    } else {
        lmax = lmax_file;
    }
    // mmax
    if (beammmax > 0) {
        mmax = beamlmax;
    } else {
        mmax = std::min((int)lmax, mmax_file);
    }
    fname = infile_beam;

    blmT_.Set(lmax, mmax);
    if (pol) {
        blmG_.Set(lmax, mmax);
        blmC_.Set(lmax, mmax);
    }

    long blmsize = (mmax + 1) * (mmax + 2) + 2 * (mmax + 1) * (lmax - mmax);

    if (rank == 0) {
        if (pol) {
            read_Alm_from_fits(infile_beam, blmT_, blmG_, blmC_, lmax, mmax);
        } else {
            read_Alm_from_fits(infile_beam, blmT_, lmax, mmax);
        }
        if (verbosity > 5) {
            std::cout << "Done reading beam alms on task = " << rank << std::endl;
        }
    }

    mpiMgr.bcastRaw(&blmT_(0, 0).re, blmsize, 0);
    if (pol) {
        mpiMgr.bcastRaw(&blmG_(0, 0).re, blmsize, 0);
        mpiMgr.bcastRaw(&blmC_(0, 0).re, blmsize, 0);
    }

    return 0;
}

Alm< xcomplex<float> > & beam::blmT(void) {
    return blmT_;
}

Alm< xcomplex<float> > & beam::blmG(void) {
    if (!pol) {
        throw std::runtime_error("Requested polarized components of an unpolarized beam");
    }

    return blmG_;
}

Alm< xcomplex<float> > & beam::blmC(void) {
    if (!pol) {
        throw std::runtime_error("Requested polarized components of an unpolarized beam");
    }

    return blmC_;
}

int sky::read(long skylmax, bool skypol,
              std::string infile_sky,
              double fwhm_deconv_sky,
              MPI_Comm comm) {

    MPI_Manager mpiMgr(comm);
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (verbosity > 5) {
        std::cout << "Reading " << infile_sky << " on task = " << rank << std::endl;
    }

    pol = skypol;
    int lmax_file, mmax_file;
    if (rank == 0) {
        if (pol) {
            get_almsize_pol(infile_sky, lmax_file, mmax_file);
        } else {
            get_almsize(infile_sky, lmax_file, mmax_file);
        }
    }
    mpiMgr.bcast(lmax_file);
    mpiMgr.bcast(mmax_file);
    // If user provided skylmax, ensure that the file has high enough
    // expansion order
    if (lmax_file < skylmax) {
        throw std::runtime_error("Sky lmax exceeds available expansion order.");
    }
    if (skylmax > 0) {
        lmax = skylmax;
    } else {
        lmax = lmax_file;
    }
    // We require equal lmax and mmax
    if (mmax_file < lmax) {
        throw std::runtime_error("Sky mmax exceeds available expansion order.");
    }
    int mmax = lmax;
    fname = infile_sky;
    fwhm_deconv = fwhm_deconv_sky;

    slmT_.Set(lmax, lmax);
    if (pol) {
        slmG_.Set(lmax, mmax);
        slmC_.Set(lmax, mmax);
    }

    long slmsize = (lmax + 1) * (lmax + 2);

    if (rank == 0) {
        if (pol) {
            read_Alm_from_fits(infile_sky, slmT_, slmG_, slmC_, lmax, mmax);
        } else {
            read_Alm_from_fits (infile_sky, slmT_, lmax, mmax);
        }
        if (verbosity > 5) {
            std::cout << "Done reading sky alms on task = " << rank << std::endl;
        }
    }

    mpiMgr.bcastRaw(&slmT_(0, 0).re, slmsize, 0);
    if (pol) {
        mpiMgr.bcastRaw(&slmG_(0, 0).re, slmsize, 0);
        mpiMgr.bcastRaw(&slmC_(0, 0).re, slmsize, 0);
    }

    double fwhm = arcmin2rad * fwhm_deconv;
    if (pol) {
        smoothWithGauss(slmT_, slmG_, slmC_, -fwhm);
    } else {
        smoothWithGauss (slmT_, -fwhm);
    }

    return 0;
}

Alm< xcomplex<float> > & sky::slmT(void) {
    return slmT_;
}

Alm< xcomplex<float> > & sky::slmG(void) {
    if (!pol) {
        throw std::runtime_error("Requested polarized components of an unpolarized sky");
    }

    return slmG_;
}

Alm< xcomplex<float> > & sky::slmC(void) {
    if (!pol) {
        throw std::runtime_error("Requested polarized components of an unpolarized sky");
    }

    return slmC_;
}

} // namespace conviqt
