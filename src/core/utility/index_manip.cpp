#include "index_manip.h"
#include <limits>
#include <cmath>
#include <iostream>

// ----------------- index manipulations ----------------
int RowFrom(int i, int l, int m, int Nmax, const std::vector<int>& Ms, const std::vector<int>& mRows){
    int mRow = RowFrom(m, Ms, mRows);

    return mRow + (l-std::abs(m))*Nmax + i;
}
void ILMFrom(int row, int& i, int& l, int& m, int Nmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int mRow;
    i = row % Nmax;               
    m = MFrom(row, Ms, mRows);              // which m-block does this row sit in
    mRow = RowFrom(m, Ms, mRows);           // what is the first row of this m-block

    row -= mRow;      // subtract off all the rows before this block
    row /= Nmax;      // get ride of i and _N factor (this is the total l-blocks behind row)

    l = row + std::abs(m);
}
int MFrom(int row, const std::vector<int>&  Ms, const std::vector<int>& mRows) {
    for (int i = 0; i < Ms.size()-1; i++)
        if (row < mRows[i+1])
            return Ms[i];
    return Ms[0]; 
}
int RowFrom(int m, const std::vector<int>&  Ms, const std::vector<int>& mRows) {
    for (int i = 0; i < Ms.size(); i++)
        if (m == Ms[i]) 
            return mRows[i];    // this is the first row of m
    return 0; 
}