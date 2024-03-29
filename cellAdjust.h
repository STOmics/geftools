/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:15
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-21 10:55:23
 * @Description: file content
 */
#ifndef GEFTOOLS_CELLADJUST_H
#define GEFTOOLS_CELLADJUST_H

#include <thread>
#include <unordered_map>
#include <unordered_set>

#include "cgef_writer.h"
#include "gef.h"
#include "opencv2/opencv.hpp"
// using namespace cv;

struct GEFTOOLS_API cellgem_label {
    cellgem_label(uint32_t geneid, int x, int y, uint32_t midcnt, uint32_t cellid) :
        geneid(geneid), x(x), y(y), midcnt(midcnt), cellid(cellid) {}
    uint32_t geneid;
    int x;
    int y;
    uint32_t midcnt;
    uint32_t cellid;
};

struct GEFTOOLS_API geneData {
    geneData(uint16_t e, uint16_t m, uint32_t c) : exon(e), mid(m), cid(c) {}
    uint16_t exon;
    uint16_t mid;
    uint32_t cid;
};

struct Pos {
    double x;
    double y;
};

class GEFTOOLS_API cellAdjust {
  public:
    cellAdjust();
    ~cellAdjust();
    cellAdjust(cellAdjust const &) = delete;
    cellAdjust &operator=(cellAdjust const &) = delete;

    void readBgef(const string &strinput);
    void readCgef(const string &strinput);
    int CompleteBgefFile();

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
    bool bexon() { return m_bexon; }

    void createRegionGef(const string &strout);
    void getRegionGenedata(vector<vector<int>> &m_vecpos);

    void readRawCgef(const string &strcgef);
    void getRegionCelldata(vector<vector<int>> &m_vecpos);
    void writeToCgef(const string &outpath);
    void writeCellToCgef();
    void writeGeneToCgef();
    void clear();

    void getSapRegion(const string &strinput, int bin, int thcnt, vector<vector<int>> &vecpos,
                      vector<sapBgefData> &vecdata, float &region_area);
    void getSapRegionIndex(const string &strinput, int bin, int thcnt, vector<vector<int>> &vecpos,
                           vector<vector<int>> &vecdata);
    void getRegionCelldataSap(vector<vector<int>> &m_vecpos);
    void getSapCellbinRegion(sapCgefData &vecdata);

    void getMultiLabelInfoFromBgef(const string &strinput, vector<vector<int>> &vecpos, vector<LabelGeneData> &vecdata,
                                   uint32_t &total_mid, int bin, int thcnt);
    void getMultiLabelInfoFromCgef(const string &strcgef, vector<vector<int>> &vecpos,
                                   vector<LabelCellData> &vecdata, vector<LabelCellDataSum> &total_data);
    void GetPositionIndexByClusterId(const char *input_file, std::vector<int> cls_id,
                                     std::vector<std::vector<int>> &clusterpos_list);
    int GenerateFilterBgefFileByMidCount(const std::string input_file, const std::string output_file, int bin_size,
                                         std::vector<MidCntFilter> filter_genes, bool only_filter);
    int GenerateFilterBgefDuration();

    void DoGenerate(int bin_size, std::vector<MidCntFilter> filter_genes, bool only_filter);

    void FilterGeneInfo(int bin_size, std::vector<MidCntFilter> filter_genes,
                        std::map<std::string, std::set<uint64_t>> &filter_data, bool only_filter);

    int GenerateBgefByLasso(const string strinput, const string stroutput, vector<vector<int>> vecpos);
    void DoLassoGenerate(const string strinput, const string stroutput, vector<vector<int>> vecpos);
    int GenerateLassoBgefDuration();

    int createRegionBgefByCord(const string &strinput, const string &strout, vector<vector<int>> &m_vecpos,
                               int bin_size);
    int createRegionCgefByCord(const string &strinput, const string &strout, vector<vector<int>> &m_vecpos);
    void setLassoBinsize(std::vector<int> bin_size);

  private:
    BgefOptions *m_bgefopts = nullptr;
    
    hid_t m_bgeffile_id = 0;
    bool m_bexon = false;
    uint32_t m_genencnt = 0;
    uint64_t m_geneexpcnt = 0;
    uint32_t m_cellcnt = 0;
    int m_offsetX = 0, m_offsetY = 0;
    vector<string> m_vecgenename;
    vector<string> m_vecgeneidname;
    int m_min_x = INT_MAX, m_min_y = INT_MAX, m_max_x = 0, m_max_y = 0;  // from /geneExp/bin1/
    uint32_t m_resolution = 0;
    unordered_map<uint64_t, vector<Dnbs_exon>> m_hash_vecdnb_exon;
    unordered_map<uint32_t, map<uint32_t, uint16_t>> m_hash_filter_cells;
    // map<uint32_t, Rect> m_hash_cellrect;
    cv::Mat m_fill_points;
    unsigned int m_block_size[4];
    CellData *m_cell_arrayptr = nullptr;
    CgefWriter *m_cgefwPtr = nullptr;
    map<uint32_t, vector<GeneExpData>> m_map_gene;

    char m_szomics[32] = {0};
    int m_maxx = 0, m_maxy = 0;
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
    int lasso_total_area_ = 0;

    cv::Mat multilabel_img;
    int cellbin_minx = INT_MAX, cellbin_miny = INT_MAX, cellbin_maxx = 0, cellbin_maxy = 0;

    std::thread generate_bgef_thread_;
    int process_rate_ = 0;

    std::thread lasso_bgef_thread_;
    int lasso_bgef_rate_ = 0;

    int cgef_version_ = 0;
    std::vector<int> lasso_binsize {};
};

#endif