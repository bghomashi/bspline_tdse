
#pragma once 

#include <vector>

// ----------------- index manipulations ----------------
int RowFrom(int i, int l, int m, int Nmax, const std::vector<int>& Ms, const std::vector<int>& mRows);
int RowFrom(int m, const std::vector<int>&  Ms, const std::vector<int>& mRows);
void ILMFrom(int row, int& i, int& l, int& m, int Nmax, const std::vector<int>& Ms, const std::vector<int>& mRows);
int MFrom(int row, const std::vector<int>&  Ms, const std::vector<int>& mRows);