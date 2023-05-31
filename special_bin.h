
#ifndef GENETOH5_SPECIALBIN_H
#define GENETOH5_SPECIALBIN_H

#include <algorithm>
#include "gef.h"
#include "bgef_options.h"


class SpecialBin
{
public:
    SpecialBin();
    ~SpecialBin();
    //int calcE10(std::vector<float> &vec_e10_result);
    int calcE10(vector<pair<string, unsigned int>>& geneCnts, std::vector<float> &vec_e10_result);
    int createPNG_py(std::vector<int> &vecdnb);
private:
    double getInverseCDFValue(double p);
    double findppf(std::vector<float> &vec, float p);
private:
    BgefOptions *m_pcmd = nullptr;
};


#endif