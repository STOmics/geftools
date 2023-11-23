/** @file gef.h
    @brief Declare GEFs related structs.

    Created by huangzhibo on 2021/12/23.
*/

#ifndef GEFTOOLS_GEF_H
#define GEFTOOLS_GEF_H

#include <cstring>
#include <string>
#include <vector>

#include "hdf5.h"
#include "utils.h"

using namespace std;

/**
 * @brief Coordinate union
 */
union GEFTOOLS_API Coordinate {
    unsigned int pos[2];  ///< dnb coordinates x, y
    unsigned long long int pos_id;
};

struct GEFTOOLS_API DnbExpression {
    int x;                 ///< dnb coordinates x
    int y;                 ///< dnb coordinates x
    unsigned short count;  ///< expression count (MIDcount)
    unsigned int gene_id;
};

struct GEFTOOLS_API Dnbs {
    Dnbs(uint32_t gid, uint16_t cnt) : geneid(gid), midcnt(cnt) {};
    uint32_t geneid;
    uint16_t midcnt;
};

struct GEFTOOLS_API Dnbs_exon {
    Dnbs_exon(uint32_t gid, uint16_t cnt, uint16_t exon) : geneid(gid), midcnt(cnt), exon(exon) {};
    uint32_t geneid;
    uint16_t midcnt;
    uint16_t exon;
};

struct GEFTOOLS_API DnbAttr {
    int min_x;
    int len_x;
    int min_y;
    int len_y;
    unsigned int max_mid;
    unsigned int max_gene;
    unsigned int max_exon;
    unsigned long number;
    int max_x;
    int max_y;
};

// wholeExp Matrix
struct GEFTOOLS_API BinStat {
    unsigned int mid_count;
    unsigned short gene_count;
};

struct GEFTOOLS_API BinStatUS {
    unsigned short mid_count;
    unsigned short gene_count;
};

struct GEFTOOLS_API DnbMatrix {
    DnbAttr dnb_attr;
    BinStatUS *pmatrix_us;
    BinStat *pmatrix;
    unsigned short *pexon16;
    unsigned int *pexon32;
};

struct GEFTOOLS_API CoordinateInfo {
    CoordinateInfo(int x, int y, unsigned int count) : x(x), y(y), count(count) {};
    int x;
    int y;
    unsigned int count;
};

/**
 * @brief Expression struct
 */
struct GEFTOOLS_API Expression {
    Expression(int x, int y, unsigned int count, unsigned int exon = 0) : x(x), y(y), count(count), exon(exon) {};
    int x;               ///< dnb coordinates x
    int y;               ///< dnb coordinates x
    unsigned int count;  ///< expression count (MIDcount)
    unsigned int exon;   // exon cnt
};

/**
 * @brief ExpressionAttr struct, record the attributes of the dataset named expressione
 */
struct GEFTOOLS_API ExpressionAttr {
    int min_x;                ///< Min X of dnb coordinate
    int min_y;                ///< Min Y of dnb coordinate
    int max_x;                ///< Max X of dnb coordinate
    int max_y;                ///< Max Y of dnb coordinate
    unsigned int max_exp;     ///< Max expression count
    unsigned int resolution;  ///< The resolution of stereo chip
};

struct GEFTOOLS_API Gene {
    Gene(const char *g, unsigned int o, unsigned c) {
        int i = 0;
        while (g[i] != '\0') {
            gene[i] = g[i];
            ++i;
        }
        offset = o;
        count = c;
    }
    char gene[64] = {0};
    unsigned int offset;
    unsigned int count;
};

struct GEFTOOLS_API GeneStat {
    GeneStat(string &g, unsigned int m, float e) {
        memcpy(gene, g.c_str(), g.size());
        mid_count = m;
        E10 = e;
    }
    GeneStat(const char *g, unsigned int m, float e) {
        int len = strlen(g);
        memcpy(gene, g, len);
        mid_count = m;
        E10 = e;
    }
    char gene[64] = {0};
    unsigned int mid_count;
    float E10;
};

// simple gene data
struct GEFTOOLS_API GeneS {
    GeneS(const char *ptr) : geneid(ptr), vecptr(nullptr) {};
    const char *geneid;
    std::vector<Expression> *vecptr;
};

struct GEFTOOLS_API GeneInfo {
    GeneInfo(const char *ptr) : geneid(ptr), umicnt(0), e10(0.0), maxexp(0), maxexon(0), vecptr(nullptr) {};
    const char *geneid;
    unsigned long umicnt;
    float e10;
    // float c50;
    unsigned int maxexp;
    unsigned int maxexon;
    std::vector<Expression> *vecptr;
};

// struct GeneErank
// {
//     GeneErank(const char *ptr):geneid(ptr){};
//     const char *geneid;
//     unsigned long umicnt;
//     float e10;
//     //float c50;
//     char attribute[10];
// };

struct GEFTOOLS_API Cell {
    unsigned int cellid;
    unsigned int offset;
    unsigned short count;
};

/**
 * @brief Describe the Cell dataset in the cell bin GEF file
 */
struct GEFTOOLS_API CellData {
    unsigned int id;
    int x;                        ///< Coordinate X of center point in this cell
    int y;                        ///< Coordinate Y of center point in this cell
    unsigned int offset;          ///< Offset of current cell in cellExp, 0-based
    unsigned short gene_count;    ///< The number of gene in this cell
    unsigned short exp_count;     ///< The total expression count of all genes in this cell
    unsigned short dnb_count;     ///< Dnb number in this cell
    unsigned short area;          ///< The polygon area of this cell
    unsigned short cell_type_id;  ///< Cell type ID to index the CellTypeList
    unsigned short cluster_id;    ///< Cluster ID, should start from 1
    // unsigned short incnt;
};

struct GEFTOOLS_API CellAttr {
    float average_gene_count;
    float average_exp_count;
    float average_dnb_count;
    float average_area;
    float median_gene_count;
    float median_exp_count;
    float median_dnb_count;
    float median_area;
    int min_x;
    int min_y;
    unsigned short min_gene_count;
    unsigned short min_exp_count;
    unsigned short min_dnb_count;
    unsigned short min_area;
    int max_x;
    int max_y;
    unsigned short max_gene_count;
    unsigned short max_exp_count;
    unsigned short max_dnb_count;
    unsigned short max_area;
};

struct GEFTOOLS_API GeneData {
    GeneData() = default;
    GeneData(const char *g, unsigned int o, unsigned int c, unsigned int e, unsigned m) {
        int i = 0;
        while (g[i] != '\0') {
            gene_name[i] = g[i];
            ++i;
        }
        offset = o;
        cell_count = c;
        exp_count = e;
        max_mid_count = m;
    }
    char gene_name[64] = {0};
    unsigned int offset {0};  ///< Offset of current gene in geneExp, 0-based
    unsigned int cell_count {0};
    unsigned int exp_count {0};
    unsigned short max_mid_count {0};  ///< max MID count of current gene
};

struct GEFTOOLS_API CellExpData {
    //    explicit CellExpData(unsigned int gene_exp){
    //        geneID = gene_exp >> 16;
    //        count = gene_exp  & 0xFFFF;
    //    }
    //    CellExpData(unsigned short g, unsigned short c){
    //        geneID = g;
    //        count = c;
    //    }
    CellExpData() = default;
    CellExpData(unsigned int id, unsigned short cnt) : gene_id(id), count(cnt) {}
    unsigned int gene_id {0};
    unsigned short count {0};
};

// For compatibility with older versions
struct GEFTOOLS_API olderCellExpData {
    olderCellExpData() = default;
    olderCellExpData(unsigned short id, unsigned short cnt) : gene_id(id), count(cnt) {}
    unsigned short gene_id {0};
    unsigned short count {0};
};

struct GEFTOOLS_API GeneExpData {
    GeneExpData(unsigned int id, unsigned short cnt) : cell_id(id), count(cnt) {}
    unsigned int cell_id;
    unsigned short count;
};

/**
 * @brief Attributes of the cell bin GEF.
 */
struct GEFTOOLS_API CellBinAttr {
    unsigned int version;     ///< Cell Bin GEF version
    unsigned int resolution;  ///< Pitch (nm) between neighbor spots
    int offsetX;              ///< Minimum value of x-axis coordinate with offset
    int offsetY;              ///< Minimum value of y-axis coordinate with offset
    string omics;
};

struct GEFTOOLS_API sapBgefData {
    sapBgefData(uint16_t g, uint32_t m, int x, int y) : genecnt(g), midcnt(m), x(x), y(y) {}
    uint16_t genecnt;
    uint32_t midcnt;
    int x;
    int y;
};

struct GEFTOOLS_API LabelGeneData {
    LabelGeneData(const char *g, unsigned int c) {
        int i = 0;
        while (g[i] != '\0') {
            gene_name[i] = g[i];
            ++i;
        }
        mid_cnt = c;
    }
    char gene_name[64] = {0};
    unsigned int mid_cnt;
};

struct GEFTOOLS_API LabelCellDataSum {
    LabelCellDataSum() = default;
    LabelCellDataSum(uint16_t a, uint16_t b, float c, uint32_t d) : cluster_id(a), mid_cnt(b), area(c), cell_id(d) {}
    uint16_t cluster_id {0};
    uint16_t mid_cnt {0};
    float area {0.0};
    uint32_t cell_id {0};
};

struct GEFTOOLS_API LabelCellData {
    LabelCellData() = default;
    LabelCellData(uint16_t a, uint16_t b, float c, uint32_t d, int e, int f) : cluster_id(a), mid_cnt(b), area(c), cell_id(d), x(e), y(f) {}
    uint16_t cluster_id {0};
    uint16_t mid_cnt {0};
    float area {0.0};
    uint32_t cell_id {0};
    int x {0};
    int y {0};
};

struct GEFTOOLS_API sapCgefData {
    unsigned int cell_count;
    float total_area;
    float average_gene_count;
    float average_exp_count;
    float average_dnb_count;
    float average_area;
    float median_gene_count;
    float median_exp_count;
    float median_dnb_count;
    float median_area;
};

struct GEFTOOLS_API ClusterPosition {
    ClusterPosition(int x, int y) : x(x), y(y) {}
    int x;
    int y;
};

struct MidCntFilter {
    std::string gene_name;
    unsigned int min_mid;
    unsigned int max_mid;
};

struct levelgenednb {
    levelgenednb() {}
    levelgenednb(float tx, float ty, uint32_t tm, float tc) : x(tx), y(ty), midcnt(tm), color(tc) {}
    float x;
    float y;
    uint32_t midcnt;
    float color;
};

hid_t getMemtypeOfGeneData();
hid_t getMemtypeOfGeneExpData();
hid_t getMemtypeOfCellData();
hid_t getMemtypeOfCellExpData();
hid_t getMemtypeOfOlderCellExpData();
bool isOlderCellExpDataVersion(hid_t fileId);
std::string getOmicsName(hid_t file_id);
bool ParseOmicsType(const std::string &bgef_file, std::string &omics);
std::string getOmicsType(const std::string &file_path, const std::string &input_omics);

#endif  // GEFTOOLS_GEF_H
