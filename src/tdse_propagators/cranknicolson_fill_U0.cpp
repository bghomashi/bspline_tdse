
#include "tdse_propagators/cranknicolson.h"
#include "utility/index_manip.h"

void CrankNicolsonTDSE::FillU0(Matrix& U0) {
    auto fOne = [](int r, int c) -> complex { return 1.; };

    // fill the main diagonal blocks
    U0->FillBandedBlock(_order-1, _N, fOne);

    if (_pol[Z]) {              // if there is z-polarization
        // fill the band above/block the main diagonal
        U0->FillBandedBlock(_order-1, _N, 1, fOne);
        U0->FillBandedBlock(_order-1, _N, -1, fOne);;
    }
    if (_pol[X] || _pol[Y]) {              // if there is x/y-polarization
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

                        U0->FillBandedBlock(_order-1, _N, blockRow, blockCol, fOne);
                    }

                    // check one l-block down
                    if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                        l2 = l1-1;
                        blockRow = m1Block + (l1-std::abs(m1));
                        blockCol = m2Block + (l2-std::abs(m2));

                        U0->FillBandedBlock(_order-1, _N, blockRow, blockCol, fOne);
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

                        U0->FillBandedBlock(_order-1, _N, blockRow, blockCol, fOne);
                    }
                    // check one l-block down
                    if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                        l2 = l1-1;
                        blockRow = m1Block + (l1-std::abs(m1));
                        blockCol = m2Block + (l2-std::abs(m2));

                        U0->FillBandedBlock(_order-1, _N, blockRow, blockCol, fOne);
                    }
                }
            }
        }

        U0->AssembleBegin();
        U0->AssembleEnd();
    }
} 