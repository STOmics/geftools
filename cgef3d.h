#ifndef GEFTOOLS_CGEF3D_H_
#define GEFTOOLS_CGEF3D_H_

#include "readFloatTask.h"

struct GEFTOOLS_API cell_3d
{
    cell_3d(uint16_t ctype, uint16_t area, uint16_t gcnt, uint16_t dcnt,
                             uint32_t id, float x, float y, float umi):
    ctypeid(ctype),area(area), genecnt(gcnt), dnbcnt(dcnt), id(id),x(x),y(y),sumumi(umi){}
    uint16_t ctypeid;
    uint16_t area;
    uint16_t genecnt;
    uint16_t dnbcnt;
    uint32_t id;
    float x;
    float y;
    float sumumi;
};

struct GEFTOOLS_API cellexp_3d
{
    cellexp_3d(uint16_t gid, float umi):gid(gid),umi(umi){}
    uint16_t gid;
    float umi;
};


struct GEFTOOLS_API gene_3d
{
    gene_3d(uint32_t off, uint32_t cellcnt, float umi, float maxumi, const char *ptr):
    offset(off),cellcnt(cellcnt),sumumi(umi),maxumi(maxumi){
        strcpy(gene,ptr);
    }
    uint32_t offset;
    uint32_t cellcnt;
    float sumumi;
    float maxumi;
    char gene[32]={0};
};

struct GEFTOOLS_API geneexp_3d
{
    geneexp_3d(uint32_t id, float umi):cid(id),umi(umi){}
    uint32_t cid;
    float umi;
};

class GEFTOOLS_API cgef3d
{
public:
    cgef3d(/* args */);
    ~cgef3d();
    void writeCgef(const string &strin, const string &strtxt, const string &strmask, const string &strout);
private:
    int gemAnalysis(const string &strinput);
    void readmask(const string &strmask);
    void readtxt(const string &strtxt);
    void readgem_5();
    //void readgem_6();
    void storeCell();
    void storeGene();
    void addCellborder(vector<float> &vec_border, vector<cv::Point2f> &vborder);
    void storeAttr(hid_t fileid);
private:
    hid_t m_gid_3d;
    ThreadPool *m_thpoolPtr = nullptr;
    unordered_map<uint32_t, std::vector<cellexp_3d>> m_hash_cell2gene;//cellid包含的gene
    unordered_map<uint32_t, uint16_t> m_hash_cell2ctype;//cellid到typeid的映射
};


#endif