/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:15
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-21 10:55:23
 * @Description: file content
 */
#ifndef GEFTOOLS_CELLADJUST_H
#define GEFTOOLS_CELLADJUST_H

#include "gef.h"
#include "cgef_writer.h"
#include <unordered_map>
#include <unordered_set>
#include "opencv2/opencv.hpp"
//using namespace cv;


struct GEFTOOLS_API cellgem_label
{
    cellgem_label(uint32_t geneid, int x, int y, uint32_t midcnt, uint32_t cellid):
    geneid(geneid),x(x),y(y),midcnt(midcnt),cellid(cellid)
    {
    }
    uint32_t geneid;
    int x;
    int y;
    uint32_t midcnt;
    uint32_t cellid;
};

struct GEFTOOLS_API geneData
{
    geneData(uint16_t e, uint16_t m, uint32_t c):exon(e),mid(m),cid(c){}
    uint16_t exon;
    uint16_t mid;
    uint32_t cid;
};

class GEFTOOLS_API cellAdjust
{
public:
    cellAdjust();
    ~cellAdjust();
    void readBgef(const string &strinput);
    void readCgef(const string &strinput);
    uint32_t getCellLabelgem(vector<string> &genename, vector<cellgem_label> &vecCellgem);
    void writeCellAdjust(const string &outpath, const string &outline_path, Cell *cellptr, int cellcnt,
                         DnbExpression *dnbptr, int dnbcnt);
    bool addborder(unsigned int cid, vector<cv::Point> &vecPoint, vector<cv::Point> &border, vector<short> &vec_border);
    bool ParseBorderFile(const string &strInput);
    bool AddBorderFromFile(unsigned int cid, vector<cv::Point> &border, vector<short> &vecBorder);
    void writeCell(Cell *cellptr, unsigned int cellcnt, DnbExpression *dnbptr, unsigned int dnbcnt);
    void writeGene();
    void cgeftogem(const string &strbgef, const string &strcgef, const string &strout);
    void cgeftogem_exon(const string &strbgef, const string &strcgef, const string &strout);
    bool bexon(){return m_bexon;}
    void createRegionGef(const string &strout);
    void getRegionGenedata(vector<vector<int>> &m_vecpos);

    void readRawCgef(const string &strcgef);
    void getRegionCelldata(vector<vector<int>> &m_vecpos);
    void writeToCgef(const string &outpath);
    void writeCellToCgef();
    void writeGeneToCgef();
    void clear();

    void getSapRegion(const string &strinput, int bin, int thcnt, vector<vector<int>> &vecpos, vector<sapBgefData> &vecdata);
    void getRegionCelldataSap(vector<vector<int>> &m_vecpos);
    void getSapCellbinRegion(sapCgefData &vecdata);

private:
    bool m_bexon = false;
    uint32_t m_genencnt;
    uint64_t m_geneexpcnt;
    uint32_t m_cellcnt;
    int m_offsetX, m_offsetY;
    vector<string> m_vecgenename;
    int m_min_x, m_min_y, m_max_x, m_max_y;
    uint32_t m_resolution;
    unordered_map<uint64_t, vector<Dnbs_exon>> m_hash_vecdnb_exon;
    unordered_map<uint32_t, map<uint32_t, uint16_t>> m_hash_filter_cells;
    //map<uint32_t, Rect> m_hash_cellrect;
    cv::Mat m_fill_points;
    unsigned int m_block_size[4];
    CellData *m_cell_arrayptr = nullptr;
    CgefWriter *m_cgefwPtr = nullptr;
    map<uint32_t, vector<GeneExpData>> m_map_gene;
    BgefOptions *m_bgefopts = nullptr;

    char m_szomics[32]={0};
    int m_maxx = 0, m_maxy = 0;
    hid_t m_bgeffile_id = 0;
    short *m_borderdataPtr = nullptr;
    vector<cellgem_label> m_vecCellgem;
    unordered_set<uint64_t> m_setcell;

    int m_effective_rect[4];
    uint16_t m_ctypecnt = 0;
    S32 *m_ctypePtr = nullptr;
    CellExpData *m_cellexpPtr = nullptr;
    // For compatibility with older versions
    olderCellExpData *m_olderCellExpPtr = nullptr;
    bool isOldCellExpVersion = false;

    GeneData *m_genePtr = nullptr;
    uint16_t *m_cellexonPtr = nullptr;
    uint16_t *m_cellexonexpPtr = nullptr;
    map<uint32_t, vector<geneData>> m_map_genedata;
    BinStat *m_parry = nullptr;
    map<unsigned int, vector<cv::Point>> borderDatas;
    bool extend_method_ = false;
    int lasso_total_area_;
};



#endif