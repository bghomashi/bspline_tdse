
#include "tdse_propagators/cranknicolson.h"
#include "utility/index_manip.h"
#include "utility/spherical_harmonics.h"
#include "utility/logger.h"

void CrankNicolsonTDSE::FillInteractionX(Matrix& HI) {
    // if we made it here we assume there IS m->m+1 coupling
    // so a full -mmax to mmax matrix
    std::vector<complex> ddr(_N*_N);
    std::vector<complex> invR(_N*_N);
    
    for (int r = 0; r < _N; r++) {
        for (int c = 0; c < _N; c++) {
            ddr[r + c*_N] = _basis.Integrate(r+1, c+1, 0, 1);
            invR[r + c*_N] = _basis.Integrate(r+1, c+1, [](complex x) {
                return 1./x;
            });
        }
	}

    // for each m-block (block rows)
    for (int i = 0; i < _Ms.size(); i++) {              
        int m1 = _Ms[i];
        int m1Block = _mRows[i]/_N;                     // number of l-blocks to skip (rows)
        int m2;
        int m2Block;                                    // number of l-blocks to skip (cols)
        int blockRow, blockCol, l2;

        if (m1+1 <= _mmax) {                            // then there is a coupling
            m2 = m1+1;
            m2Block = RowFrom(m2, _Ms, _mRows)/_N;
            // for each l-block (block rows)
            for (int l1 = std::abs(m1); l1 <= _lmax; l1++) {
                // check one l-block up
                if (l1+1 <= _lmax && l1+1 >= std::abs(m2)) {
                    l2 = l1+1;
                    blockRow = m1Block + (l1-std::abs(m1));
                    blockCol = m2Block + (l2-std::abs(m2));

                    HI->FillBandedBlock(_order-1, _N, blockRow, blockCol, 
                    [=,&ddr,&invR](int row, int col) {
                        int i = row % _N, j = col % _N;                     // TODE: change i,j to anything else
                        
                        // m1=m2-1, l1=l2-1
                        //double a = -sqrt((l2+m2) * (l2+m2-1) / (2.*l2 + 1.) / (2.*l2 - 1.));
                        double a = -sqrt((l2+m2) * (l2+m2-1) / (2.*l2 + 1.) / (2.*l2 - 1.));

                        return (ddr[i + j*_N] + double(l2)*invR[i + j*_N])*YlmXYlm(l1,m1,l2,m2);
                    });
                }

                // check one l-block down
                if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                    l2 = l1-1;
                    blockRow = m1Block + (l1-std::abs(m1));
                    blockCol = m2Block + (l2-std::abs(m2));

                    HI->FillBandedBlock(_order-1, _N, blockRow, blockCol, 
                    [=,&ddr,&invR](int row, int col) {
                        int i = row % _N, j = col % _N;
                        
                        // m1=m2-1, l1=l2+1
                        double a = sqrt((l2-m2+1)*(l2-m2+2) / (2.*l2 + 1.) / (2.*l2 + 3.));

                        return (ddr[i + j*_N] - double(l2+1)*invR[i + j*_N])*YlmXYlm(l1,m1,l2,m2);;
                    });
                }
            }
        }
        if (m1-1 >= -_mmax) {                           // then there is a coupling
            m2 = m1-1;
            m2Block = RowFrom(m2, _Ms, _mRows)/_N;
            // for each l-block (block rows)
            for (int l1 = std::abs(m1); l1 <= _lmax; l1++) {
                // check one l-block up
                if (l1+1 <= _lmax && l1+1 >= std::abs(m2)) {
                    l2 = l1+1;
                    blockRow = m1Block + (l1-std::abs(m1));
                    blockCol = m2Block + (l2-std::abs(m2));

                    HI->FillBandedBlock(_order-1, _N, blockRow, blockCol, 
                    [=,&ddr,&invR](int row, int col) {
                        int i = row % _N, j = col % _N;
                        
                        // m1=m2+1, l1=l2-1
                        double a = sqrt((l2-m2)*(l2-m2-1) / (2.*l2 + 1.) / (2.*l2 - 1.));

                        return (ddr[i + j*_N] + double(l2)*invR[i + j*_N])*YlmXYlm(l1,m1,l2,m2);;
                    });
                }
                // check one l-block down
                if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                    l2 = l1-1;
                    blockRow = m1Block + (l1-std::abs(m1));
                    blockCol = m2Block + (l2-std::abs(m2));

                    HI->FillBandedBlock(_order-1, _N, blockRow, blockCol, 
                    [=,&ddr,&invR](int row, int col) {
                        int i = row % _N, j = col % _N;
                        
                        // m1=m2+1, l1=l2+1
                        double a = -sqrt((l2+m2+1)*(l2+m2+2) / (2.*l2 + 1.) / (2.*l2 + 3.));

                        return (ddr[i + j*_N] - double(l2+1)*invR[i + j*_N])*YlmXYlm(l1,m1,l2,m2);;
                    });
                }
            }
        }
    }

    HI->AssembleBegin();
    HI->AssembleEnd();

    // HI->Scale(-1.);

    // MatView(std::dynamic_pointer_cast<PetscMatrix>(HI)->_petsc_mat, 0);

    // if (HI->IsAntiSymmetric(1e-12))
    //     Log::info("Is Anti  Symmetric!");
    // else
    //     Log::info("Is not Anti Symmetric!");
    // exit(0);
}
