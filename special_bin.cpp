#include "special_bin.h"
#include<stdlib.h>
#include <map>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#define NOMINMAX
#include <windows.h>
#else
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>
#endif
#include <cstdio>
#include <sstream>

SpecialBin::SpecialBin(/* args */)
{
    m_pcmd = BgefOptions::GetInstance();
}

SpecialBin::~SpecialBin()
{
}

double SpecialBin::getInverseCDFValue(double p) 
{
    double a1 = -39.69683028665376;
    double a2 = 220.9460984245205;
    double a3 = -275.9285104469687;
    double a4 = 138.3577518672690;
    double a5 =-30.66479806614716;
    double a6 = 2.506628277459239;

    double b1 = -54.47609879822406;
    double b2 = 161.5858368580409;
    double b3 = -155.6989798598866;
    double b4 = 66.80131188771972;
    double b5 = -13.28068155288572;

    double c1 = -0.007784894002430293;
    double c2 = -0.3223964580411365;
    double c3 = -2.400758277161838;
    double c4 = -2.549732539343734;
    double c5 = 4.374664141464968;
    double c6 = 2.938163982698783;

    double d1 = 0.007784695709041462;
    double d2 = 0.3224671290700398;
    double d3 = 2.445134137142996;
    double d4 = 3.754408661907416;

    double p_low =  0.02425;
    double p_high = 1 - p_low;
    long double  q, r, e, u;
    long double x = 0.0;

    if (0 < p && p < p_low) {
        q = sqrt(-2*log(p));
        x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
    }

    if (p_low <= p && p <= p_high) {
        q = p - 0.5;
        r = q*q;
        x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
    }

    if (p_high < p && p < 1) {
        q = sqrt(-2*log(1-p));
        x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
    }

    if(( 0 < p)&&(p < 1)){
        e = 0.5 * erfc(-x/sqrt(2)) - p;
        u = e * sqrt(2*M_PI) * exp(x*x/2);
        x = x - u/(1 + x*u/2);
    }
    return x;
}

double SpecialBin::findppf(std::vector<float> &vec, float p) //已知概率求分位点
{
    int sz = vec.size();
    double sum = 0.0;
    for(float f : vec)
    {
        sum += f;
    }
    double mean = sum/sz; //平均数
    double accum  = 0.0; 
    for(float f : vec)
    {
        accum += (f-mean)*(f-mean);  
    }
    double stdev = sqrt(accum/(sz-1)); //方差
    double ppf = getInverseCDFValue(p)*stdev + mean;
    return ppf;
}

int SpecialBin::calcE10(vector<pair<string, unsigned int>>& geneCnts, std::vector<float> &vec_e10_result)
{
    // map<string, int> genePos;
    // int i = 0;
    // for (auto& p : geneCnts)
    //     genePos[p.first] = i++;

    // std::stringstream ss;
    // ss.setf(std::ios::fixed);
    // ss.precision(2);
    // vec_e10_result.resize(geneCnts.size());
    // auto itor_rank = m_pcmd->vec_bin100_.begin();
    // for(;itor_rank != m_pcmd->vec_bin100_.end();itor_rank++)
    // {
    //     auto pos = genePos.find(itor_rank->geneid);
    //     if (pos == genePos.end()) continue;

    //     ss.str("");
    //     ss << itor_rank->e10;
    //     vec_e10_result[pos->second] = stof(ss.str());
    // }
    return 0;
}

// int SpecialBin::calcE10(std::vector<float> &vec_e10_result)
// {
//     std::vector<float> vec_e10;
//     std::vector<float> vec_c50;
//     auto itor_rank = m_pcmd->vec_bin100_.begin();
//     for(;itor_rank != m_pcmd->vec_bin100_.end();itor_rank++)
//     {
//         if(itor_rank->umicnt > 300)
//         {
//             memcpy(itor_rank->attribute, "non", 3);
//             vec_e10.emplace_back(itor_rank->e10);
//             vec_c50.emplace_back(itor_rank->c50);
//             //printf("%d  ", itor_rank->umicnt);
//         }
//         else
//         {
//             memcpy(itor_rank->attribute, "low", 3);
//         }
//     }

//     std::vector<GeneErank*> vec_pattern;
//     std::vector<GeneErank*> vec_non;
//     std::vector<GeneErank*> vec_low;
//     double e10_ppf = findppf(vec_e10, 0.9);
//     double c50_ppf = findppf(vec_c50, 0.1);
//     //printf("e10:%7.5f c50:%7.5f\n", e10_ppf, c50_ppf);
//     for(int i=0;i<m_pcmd->vec_bin100_.size();i++)
//     {
//         GeneErank & grank = m_pcmd->vec_bin100_[i];
//         if(grank.umicnt > 300 && grank.e10 > e10_ppf && grank.c50 < c50_ppf)
//         {
//             memcpy(grank.attribute, "pattern", 7);
//         }

//         if(memcmp(grank.attribute, "pattern", 7) == 0)
//         {
//             vec_pattern.emplace_back(&grank);
//         }
//         else if (memcmp(grank.attribute, "non", 3) == 0)
//         {
//             vec_non.emplace_back(&grank);
//         }
//         else
//         {
//             vec_low.emplace_back(&grank);
//         }
//     }

//     std::sort(vec_pattern.begin(), vec_pattern.end(), [](const GeneErank *a, const GeneErank *b){
//         return a->e10 > b->e10;
//         });
//     std::sort(vec_non.begin(), vec_non.end(), [](const GeneErank *a, const GeneErank *b){
//         return a->e10 > b->e10;
//         });
//     std::sort(vec_low.begin(), vec_low.end(), [](const GeneErank *a, const GeneErank *b){
//         return a->e10 > b->e10;
//         });
    

//     char buf[128]={0};
//     int len = 0;
//     string str_p, str_n, str_l;
//     vec_e10_result.reserve(m_pcmd->vec_bin100_.size());
//     for(GeneErank *ge : vec_pattern)
//     {
//         vec_e10_result.emplace_back(ge->e10);
//         // len = sprintf(buf, "%s\t%7.2f\t%7.2f\t pattern\n", ge->geneid, ge->e10, ge->c50);
//         // buf[len]='\0';
//         // str_p.append(buf);
//     }

//     for(GeneErank *ge : vec_non)
//     {
//         vec_e10_result.emplace_back(ge->e10);
//         // len = sprintf(buf, "%s\t%7.2f\t%7.2f\t non\n", ge->geneid, ge->e10, ge->c50);
//         // buf[len]='\0';
//         // str_n.append(buf);
//     }

//     for(GeneErank *ge : vec_low)
//     {
//         vec_e10_result.emplace_back(ge->e10);
//         // len = sprintf(buf, "%s\t%7.2f\t%7.2f\t low\n", ge->geneid, ge->e10, ge->c50);
//         // buf[len]='\0';
//         // str_l.append(buf);
//     }

//     // FILE *fout = fopen("erank.txt", "w");
//     // if(fout)
//     // {
//     //     fwrite(str_p.c_str(), 1, str_p.length(), fout);
//     //     fwrite(str_n.c_str(), 1, str_n.length(), fout);
//     //     fwrite(str_l.c_str(), 1, str_l.length(), fout);
//     //     fclose(fout);
//     // }
//     return 0;
// }

//char * getDirectory( char * buf, int count)
//{
    //int i;
    //int rslt = readlink("/proc/self/exe", buf, count - 1);
    //if (rslt < 0 || (rslt >= count - 1))
    //{
    //    return NULL;
    //}
    //buf[rslt] = '\0';
    //for (i = rslt; i >= 0; i--)
    //{
    //    if (buf[i] == '/')
    //    {
    //        buf[i + 1] = '\0';
    //        break;
    //    }
    //}
    //return buf;
//}

int SpecialBin::createPNG_py(std::vector<int> &vecdnb)
{
    //unsigned long x, y, cnt = 0;
    //std::string str_x, str_y, str_cnt;
    //char buf[32]={0};
    //int len = 0;
    //for(int i=0;i<vecdnb.size();i+=3)
    //{
    //    x = vecdnb[i];
    //    y = vecdnb[i+1];
    //    cnt = vecdnb[i+2];

    //    len = sprintf(buf, "%lu\t", x);
    //    buf[len] = '\0';
    //    str_x.append(buf);

    //    len = sprintf(buf, "%lu\t", y);
    //    buf[len] = '\0';
    //    str_y.append(buf);

    //    len = sprintf(buf, "%lu\t", cnt);
    //    buf[len] = '\0';
    //    str_cnt.append(buf);
    //}  
    //str_x.append("\n");
    //str_y.append("\n");
    //str_cnt.append("\n");
    //
    //std::string path(m_pcmd->output_file_);
    //path.append("/dnb_merge/");
    //mkdir(path.c_str(), 0777);
    //path.append("bin200.png");

    //std::string tmpfile(m_pcmd->output_file_);
    //tmpfile.append("/dnb_merge/.pngtmp");
    //FILE *fout = fopen(tmpfile.c_str(), "w");
    //if(fout)
    //{
    //    fwrite(str_x.c_str(), 1, str_x.length(), fout);
    //    fwrite(str_y.c_str(), 1, str_y.length(), fout);
    //    fwrite(str_cnt.c_str(), 1, str_cnt.length(), fout);
    //    fclose(fout);
    //}

    //char binpath[1024];
    //getDirectory(binpath, 1024);
    //std::string cmd("python3 ");
    //cmd += std::string(binpath) + "/png.py " + tmpfile + " " + path;
    //system(cmd.c_str());

    //remove(tmpfile.c_str());
    return 0;
}
