/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:15:25
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-19 19:11:09
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFCELLGEM_H_
#define GEFTOOLS_CGEFCELLGEM_H_

#include "cgef_writer.h"
#include "thread_pool.h"
#include "opencv2/opencv.hpp"
#include "cgefUtil.h"
//using namespace cv;

struct GEFTOOLS_API celldata
{
    celldata(int cid, int lid):c_idx(cid),l_idx(lid){}
    int c_idx; //轮廓idx
    int l_idx; //连通域idx
};

struct GEFTOOLS_API cellt
{
    cellt(uint16_t m, uint16_t e, uint32_t c):
    mid(m),exon(e),cid(c){}
    uint16_t mid;
    uint16_t exon;
    uint32_t cid;
};


class GEFTOOLS_API cgefCellgem
{
public:
    cgefCellgem();
    ~cgefCellgem();
    //void readmask(const string &strmask);
    void writeFile(CgefWriter *cwptr, const string &strmask, const string &strinput);
    void writeAttr();
    // void writeCell();
    // void writeGene();
    // void getCelldata();
    // void getCelldata_celltype();
    // void writeCell_celltype();
    void addCellborder(int cx, int cy, vector<short> &vec_border, uint32_t idx);

    //void readcellgem_new();

    // void clabeltocid();
    // void writeGene_raw();
    // void writeCell_raw();

    void readBgef(const string &strinput);
    void writeGene_bgef();


    void readBgef_new(const string &strinput);
    void readmask_new(const string &strmask);
    void getCell();
    void writeCell_new();
    void writeGene_new();
    void gemPreAnalysis(const string &strmask, const string &strinput);
    // void readgem_4mask();
    // void readgem_5mask();
    // void readgem_5();
    // void readgem_6mask();
    // void readgem_6();
    // void readgem_6type();


    void cgem2cgef(CgefWriter *cwptr, const string &strin);
    void getCelldata_cgem();
    void writeCell_cgem();
    void writeGene_cgem();

public:
    unsigned int m_block_size[4] = {0};
    cv::Mat m_stats, m_outimg, m_centroids;
    unordered_map<uint64_t, vector<cellExp_Exon>> m_hash_vecdnb;
    GefQueue<cellUnit> *m_cellqueuePtr = nullptr;
private:
    bool m_bexon = false;
    unsigned int m_maskcellnum = 0; //从mask文件获取的细胞个数
    unsigned int m_blocknum;
    uint32_t m_labelcnt = 0;

    int m_min_x{INT_MAX}, m_max_x{0}, m_min_y{INT_MAX}, m_max_y{0}, m_rows{0}, m_cols{0};
    
    vector<vector<cv::Point>> m_contours;
    vector<vector<celldata>> m_vec_veccell;
    vector<vector<uint32_t>> m_vec_veccid;

    CgefWriter *m_cgefwPtr = nullptr;
    ThreadPool *m_thpoolPtr = nullptr;
    unordered_map<uint32_t, uint32_t> m_hash_clabel2cid; //建立从label到cellid的映射
    unordered_map<string, int> m_hash_gname2gid;//gname到geneid的映射
    unordered_map<string, int> m_hash_celltype;//celltype到typeid的映射
    vector<string> m_vec_celltype;

    vector<uint32_t> m_vec_blkidx;
    vector<uint32_t> m_vec_cellLabel;
    //map<uint32_t, bgef_cell *> m_mapcellexp;
    vector<bgef_cell *> m_vec_cellexp;

    uint32_t m_genecnt = 0;
    uint64_t m_geneExpcnt = 0;
    Gene *m_genePtr = nullptr;
    Expression *m_expPtr = nullptr;
    
    vector<vector<cellUnit*>> m_vec_vec_cellunit;
    unordered_map<uint32_t, geneUnit*> m_hash_geneunit; //key:gid 
    string m_stromics{"Transcriptomics"};

    uint32_t m_borcnt = 0;
    unordered_map<uint32_t, vector<cellt>> m_map_gene;//记录基因占据的cell信息
};

#endif