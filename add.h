#ifndef KURSACH_ADD_H
#define KURSACH_ADD_H

#include <vector>
#include "struct.h"

using namespace std;

void add_to_global(
    vector<double>& ggl,
    vector<double>& ggu,
    vector<double>& di,
    vector<int>& ig,
    vector<int>& jg,
    vector<vector<double>>& local_A,
    element element)
{
    for (int j = 0; j < local_A.size(); j++)
    {
        di[element.joints[j]] += local_A[j][j];
        int i_beg = ig[element.joints[j]];
        for (int k = 0; k < j; k++) {
            int i_end = ig[element.joints[j] + 1] - 1;
            while (jg[i_beg] != element.joints[k]) {
                int ind = (i_beg + i_end) / 2;
                if (jg[ind] < element.joints[k]) {
                    i_beg = ind + 1;
                }
                else {
                    i_end = ind;
                }
            }
            ggl[i_beg] += local_A[j][k];
            ggu[i_beg] += local_A[k][j];
            i_beg++;
        }
    }
}

#endif //KURSACH_ADD_H
