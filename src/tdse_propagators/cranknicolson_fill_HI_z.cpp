
#include "tdse_propagators/cranknicolson.h"
#include "utility/index_manip.h"
#include "utility/logger.h"
#include "math_libs/petsc/petsc_lib.h"

void CrankNicolsonTDSE::FillInteractionZ(Matrix& HI) {
    std::vector<complex> ddr(_N*_N);
    std::vector<complex> invR(_N*_N);
    
    // cache some common matrix elements
    for (int r = 0; r < _N; r++) {
        for (int c = 0; c < _N; c++) {
            ddr[r + c*_N] = _basis.Integrate(r+1, c+1, 0, 1);                   // <Bi|d/dr|Bj>
            invR[r + c*_N] = _basis.Integrate(r+1, c+1, [](complex x) {         // <Bi|1/r|Bj>
                return 1./x;
            });
        }
	}

    for (int i = 0; i < _Ms.size(); i++) {              
        int m1 = _Ms[i];
        int m1Block = _mRows[i]/_N;                     // number of l-blocks to skip (rows)
        int m2 = m1;
        int m2Block = m1Block;                          // number of l-blocks to skip (cols)
        int blockRow, blockCol, l2;


        for (int l1 = std::abs(m1); l1 <= _lmax; l1++) {
            // check one l-block up
            if (l1+1 <= _lmax && l1+1 >= std::abs(m2)) {
                l2 = l1+1;
                blockRow = m1Block + (l1-std::abs(m1));
                blockCol = m2Block + (l2-std::abs(m2));
                // std::cout << "m1 = " << m1;
                // std::cout << " m2 = " << m1;
                // std::cout << " l1 = " << l1;
                // std::cout << " l2 = " << l2;
                // std::cout << " a = " << sqrt((l2+m2) * (l2-m2) / (2.*l2 + 1.) / (2.*l2 - 1.)) << std::endl;;

                HI->FillBandedBlock(_order-1, _N, blockRow, blockCol, 
                [=,&ddr,&invR](int row, int col) {
                    int i = row % _N, j = col % _N;
                    
                    // m1=m2, l1=l2-1
                    double a = sqrt((l2+m2) * (l2-m2) / (2.*l2 + 1.) / (2.*l2 - 1.));

                    return (ddr[i + j*_N] + double(l2)*invR[i + j*_N])*a;
                });
            }
            // check one l-block down
            if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                l2 = l1-1;
                blockRow = m1Block + (l1-std::abs(m1));
                blockCol = m2Block + (l2-std::abs(m2));

                // std::cout << "m1 = " << m1;
                // std::cout << " m2 = " << m1;
                // std::cout << " l1 = " << l1;
                // std::cout << " l2 = " << l2;
                // std::cout << " a = " << sqrt((l2+m2+1.)*(l2 - m2 + 1.) / (2.*l2 + 1.) / (2.*l2 + 3.)) << std::endl;

                HI->FillBandedBlock(_order-1, _N, blockRow, blockCol, 
                [=,&ddr,&invR](int row, int col) {
                    int i = row % _N, j = col % _N;
                    // m1=m2, l1=l2+1
                    double a = sqrt((l2+m2+1.)*(l2 - m2 + 1.) / (2.*l2 + 1.) / (2.*l2 + 3.));

                    return (ddr[i + j*_N] - double(l2+1)*invR[i + j*_N])*a;
                });
            }
        }
    }
    // Fill L1 = l2 - 1 blocks
    // HI->FillBandedBlock(_order-1, _N, 1, [&](int row, int col) {
    //     int i, j, l1, l2, m1, m2;
    //     ILMFrom(row, i, l1, m1, _N, _Ms, _mRows);
    //     ILMFrom(col, j, l2, m2, _N, _Ms, _mRows);
        
    //     double a = sqrt((l2+m2) * (l2-m2) / (2.*l2 + 1.) / (2.*l2 - 1.));

    //     return (ddr[i + j*_N] + double(l2)*invR[i + j*_N])*a;
    // });
    // // Fill L1 = l2 + 1 blocks
    // HI->FillBandedBlock(_order-1, _N, -1, [&](int row, int col)  {
    //     int i, j, l1, l2, m1, m2;
    //     ILMFrom(row, i, l1, m1, _N, _Ms, _mRows);
    //     ILMFrom(col, j, l2, m2, _N, _Ms, _mRows);
        
    //     double a = sqrt((l2 + m2 + 1.)*(l2 - m2 + 1.) / (2.*l2 + 1.) / (2.*l2 + 3.));

    //     return (ddr[i + j*_N] - (l2+1.)*invR[i + j*_N])*a;
    // });
    HI->AssembleBegin();
    HI->AssembleEnd();


    // MatView(std::dynamic_pointer_cast<PetscMatrix>(HI)->_petsc_mat, 0);
    // exit(0);


}
