#ifndef GEFTOOLS_GEFTOGEM_H
#define GEFTOOLS_GEFTOGEM_H

#include "gef.h"
#include <unordered_map>
#include "opencv2/opencv.hpp"
#include "utils.h"
//using namespace cv;

struct GEFTOOLS_API cellmat
{
    int offsetx;
    int offsety;
    vector<cv::Point> vecP;
};


class GEFTOOLS_API geftogem
{
public:
    geftogem(const string &strout, const string &strsn, bool boutexon);
    ~geftogem();

    void bgeftogem(const string &strbgef, int binsize = 1);
    void cgeftogem(const string &strcgef, const string &strbgef);
    void bgeftocgem(const string &strmask, const string &strbgef);
private:
    void readBgef(const string &strinput);
    void readCgef(const string &strinput);
    void getBgefGene(hid_t file_id);
    void getBgefExp(hid_t file_id);
    void getdnb();
    void bgef2gem();
    void cgef2gem();
    void cgef2gem_exon();
    void readmask(const string &strmask);
private:
    bool m_bexon = false;
    bool m_boutexon = true;
    int m_bin = 1;
    Gene *m_genePtr = nullptr;
    Expression *m_expPtr = nullptr;
    uint32_t m_genencnt;
    uint64_t m_geneexpcnt;
    uint32_t m_cellcnt;
    int m_offsetX, m_offsetY;
    vector<string> m_vecgenename;
    int m_min_x, m_min_y, m_max_x, m_max_y;
    uint32_t m_resolution;
    unordered_map<uint64_t, vector<Dnbs>> m_hash_vecdnb;
    unordered_map<uint64_t, vector<Dnbs_exon>> m_hash_vecdnb_exon;
    string m_strout; //输出gem路径
    string m_strsn; //gem sn号
    //Mat m_fill_points;
    unordered_map<uint32_t, cellmat> m_hash_cellpoint;
};



#endif