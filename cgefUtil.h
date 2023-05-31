/*
 * @Author: zhaozijian
 * @Date: 2022-04-01 10:15:19
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-20 09:02:29
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFUTIL_H_
#define GEFTOOLS_CGEFUTIL_H_

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include "opencv2/opencv.hpp"
#include "utils.h"

struct cgef_gdata
{
    cgef_gdata(){};
    cgef_gdata(unsigned short u, unsigned short e):umi(u),exon(e){};
    unsigned short umi;
    unsigned short exon;
};

class cgef_cell
{
public:
    cgef_cell(int label):m_celllabel(label){}
    cgef_cell(int label, char *ptr,int len):
    m_celllabel(label)
    {
        memcpy(m_type, ptr, len);
        m_type[len]='\0';
    }
    ~cgef_cell(){}
    bool add(std::string &gene, unsigned short umi)
    {
        dnbcnt++;
        expcnt+=umi;

        if(m_map_cellexp.find(gene) != m_map_cellexp.end())
        {
            m_map_cellexp[gene].umi += umi;
        }
        else
        {
            cgef_gdata cg{umi, 0};
            m_map_cellexp.emplace(gene, std::move(cg));
        }
        return true;
    }

    bool add(std::string &gene, unsigned short umi, int x, int y)
    {
        m_vecPoint.emplace_back(x,y);
        dnbcnt++;
        expcnt+=umi;

        if(m_map_cellexp.find(gene) != m_map_cellexp.end())
        {
            m_map_cellexp[gene].umi += umi;
        }
        else
        {
            cgef_gdata cg{umi, 0};
            m_map_cellexp.emplace(gene, std::move(cg));
        }
        return true;
    }

    bool add(std::string &gene, unsigned short umi, int x, int y, unsigned short exon)
    {
        m_vecPoint.emplace_back(x,y);
        dnbcnt++;
        expcnt+=umi;
        m_exon+=exon;

        if(m_map_cellexp.find(gene) != m_map_cellexp.end())
        {
            m_map_cellexp[gene].umi += umi;
            m_map_cellexp[gene].exon += exon;
        }
        else
        {
            cgef_gdata cg{umi, exon};
            m_map_cellexp.emplace(gene, std::move(cg));
        }
        return true;
    }

    bool merge(cgef_cell &cell)
    {
        dnbcnt += cell.dnbcnt;
        expcnt += cell.expcnt;
        m_exon += cell.m_exon;

        auto itor = cell.m_map_cellexp.begin();
        for(;itor!= cell.m_map_cellexp.end();itor++)
        {
            if(m_map_cellexp.find(itor->first) != m_map_cellexp.end())
            {
                m_map_cellexp[itor->first].umi += itor->second.umi;
                m_map_cellexp[itor->first].exon += itor->second.exon;
            }
            else
            {
                cgef_gdata cg{itor->second.umi, itor->second.exon};
                m_map_cellexp.emplace(itor->first, std::move(cg));
            }
        }

        m_vecPoint.insert(m_vecPoint.end(), cell.m_vecPoint.begin(), cell.m_vecPoint.end());
        return true;
    }

    // bool getCenter_median(unsigned int *block_size, int offx, int offy) //中位数
    // {
    //     int sz = m_vecPoint.size();
    //     if(sz<3) return false;

    //     vector<int> vecx; vecx.reserve(sz);
    //     vector<int> vecy; vecy.reserve(sz);
    //     for(Point &pt : m_vecPoint)
    //     {
    //         vecx.emplace_back(pt.x);
    //         vecy.emplace_back(pt.y);
    //     }

    //     std::sort(vecx.begin(), vecx.end(), [](int a, int b){return a<b;});
    //     std::sort(vecy.begin(), vecy.end(), [](int a, int b){return a<b;});

    //     int pos = ceil( (sz + 1)*1.0 /2 );
    //     m_center.x = ceil (vecx[pos-2]*0.5+vecx[pos-1]*0.5);
    //     m_center.y = ceil (vecy[pos-2]*0.5+vecy[pos-1]*0.5);
        
    //     m_blkid = (m_center.x-offx)/block_size[0] + ((m_center.y-offy)/block_size[1])*block_size[2];
    //     assert(m_blkid < block_size[2] * block_size[3]);
    //     return true;
    // }

    bool getCenter_border(unsigned int *block_size, int offx, int offy)//凸多边形
    {
        if(m_vecPoint.size() <3) return false;
        vector<cv::Point> tmp, hull;
        convexHull(m_vecPoint, hull, true);
        m_vecPoint.swap(tmp);
        int sz = hull.size();
        if(sz <= 2)
        {
            return false;
        }
        else
        {
            if(sz > BORDERCNT)
            {
                double epsilon = 0.01 * arcLength(hull, true);
                approxPolyDP(hull, m_border, epsilon, true); 
            }
            else
            {
                m_border.swap(hull);
            }

            cv::Moments mu = cv::moments(m_border, true);
            if(mu.m00 == 0) return false;    
            m_center = cv::Point(static_cast<int>(mu.m10/mu.m00), static_cast<int>(mu.m01/mu.m00));
            m_area = mu.m00;
            
        }

        m_blkid = (m_center.x-offx)/block_size[0] + ((m_center.y-offy)/block_size[1])*block_size[2];
        assert(m_blkid < block_size[2] * block_size[3]);
        return true;
    }

public:
    cv::Point m_center;
    vector<cv::Point> m_vecPoint;
    vector<cv::Point> m_border;
    uint32_t m_blkid = 0;
    int m_celllabel = 0;
    unsigned short expcnt = 0; //细胞umi之和
    unsigned short dnbcnt = 0; //细胞所有坐标点个数 >= gene个数
    unsigned short m_exon = 0; //细胞exon之和
    char m_type[32]={0};
    uint16_t m_area = 0;
    std::map<std::string, cgef_gdata> m_map_cellexp;//包含的基因情况
};


class cgef_gene
{
public:
    cgef_gene(){}
    ~cgef_gene(){}
    // bool add(int label, unsigned short umi)
    // {
    //     expcnt += umi;
    //     if(m_map_geneexp.find(label) != m_map_geneexp.end())
    //     {
    //         m_map_geneexp[label].umi += umi;
    //     }
    //     else
    //     {
    //         m_map_geneexp.emplace(label, umi, 0);
    //     }
    //     return true;
    // }

    // bool add(int label, unsigned short umi, unsigned short exon)
    // {
    //     expcnt += umi;
    //     exonsum += exon;
    //     if(m_map_geneexp.find(label) != m_map_geneexp.end())
    //     {
    //         m_map_geneexp[label].umi += umi;
    //         m_map_geneexp[label].exon += exon;
    //     }
    //     else
    //     {
    //         m_map_geneexp.emplace(label, umi, exon);
    //     }
    //     return true;
    // }

    // bool merge(cgef_gene &gene)
    // {
    //     expcnt += gene.expcnt;
    //     exonsum += gene.exon;
    //     auto itor = gene.m_map_geneexp.begin();
    //     for(;itor != gene.m_map_geneexp.end();itor++)
    //     {
    //         if(m_map_geneexp.find(itor->first) != m_map_geneexp.end())
    //         {
    //             m_map_geneexp[itor->first].umi += itor->second.umi;
    //             m_map_geneexp[itor->first].exon += itor->second.exon;
    //         }
    //         else
    //         {
    //             m_map_geneexp.emplace(itor->first, itor->second);
    //         }
    //     }
    //     return true;
    // }
public:
    unsigned int expcnt = 0; //基因umi之和
    unsigned int exonsum = 0;
    std::map<int, cgef_gdata> m_map_geneexp;//出现的细胞情况
};

class bgef_gene
{
public:
    bgef_gene(){};
    ~bgef_gene(){};
    void add(uint32_t x, uint32_t y, uint32_t midcnt)
    {
        m_vecExp.emplace_back(x,y,midcnt);
    }
    void merge(bgef_gene &other)
    {
        m_vecExp.insert(m_vecExp.end(), other.m_vecExp.begin(), other.m_vecExp.end());
    }
public:
    vector<Expression> m_vecExp;
};

class bgef_cell
{
public:
    bgef_cell(uint32_t cid, int x, int y, uint16_t area, uint32_t label):
    m_cid(cid), m_cx(x), m_cy(y), m_area(area),m_clabel(label)
    {};
    ~bgef_cell(){};

    void add(uint16_t gid, uint16_t expcnt, uint16_t dnbcnt, uint16_t exoncnt)
    {
        m_vecCexp.emplace_back(gid, expcnt);
        m_vecCExon.emplace_back(exoncnt);
        m_expcnt += expcnt;
        m_dnbcnt += dnbcnt;
        m_exoncnt += exoncnt;
        m_maxexon = std::max(m_maxexon, exoncnt);
    }
public:
    uint16_t m_maxexon = 0;
    uint16_t m_exoncnt = 0;
    uint16_t m_expcnt = 0;
    uint16_t m_dnbcnt = 0;
    uint16_t m_area = 0;
    int m_cx,m_cy;
    uint32_t m_cid;
    uint32_t m_clabel = 0;
    vector<CellExpData> m_vecCexp;
    vector<uint16_t> m_vecCExon;
};

////////////////////////////////////////////////////////////////

struct gExp
{
    // gExp(uint16_t mid, uint16_t exon):
    // midcnt(mid),exoncnt(exon){};
    uint16_t midcnt;
    uint16_t exoncnt;
};

struct cellExp_Exon
{
    cellExp_Exon(uint32_t gid, uint16_t mid, uint16_t exon):
    geneid(gid),midcnt(mid),exoncnt(exon){};
    uint32_t geneid;
    uint16_t midcnt;
    uint16_t exoncnt;
};

class cellUnit
{
public:
    cellUnit(int x, int y, uint16_t area, uint32_t label, unsigned int *block_size):
    m_cx(x), m_cy(y), m_area(area), m_label(label)
    {
        m_blkid = x/block_size[0] + (y/block_size[1])*block_size[2];
    }
    ~cellUnit(){};
    void add(vector<cellExp_Exon> &vecdnb)
    {
        for(cellExp_Exon &dnb : vecdnb)
        {
            if(m_map_gExp.find(dnb.geneid) == m_map_gExp.end())
            {
                gExp tg{0,0};
                m_map_gExp.emplace(dnb.geneid, tg);
            }
            m_map_gExp[dnb.geneid].midcnt += dnb.midcnt;
            m_map_gExp[dnb.geneid].exoncnt += dnb.exoncnt;
            
            m_dnbcnt++;
            m_expcnt+=dnb.midcnt;
            m_exoncnt+=dnb.exoncnt;
        }
    }
public:
    uint16_t m_expcnt = 0;
    uint16_t m_dnbcnt = 0;
    uint16_t m_area = 0;
    uint16_t m_exoncnt = 0;
    int m_cx = 0;
    int m_cy = 0;
    uint32_t m_label = 0;
    uint32_t m_blkid = 0;
    std::map<uint32_t, gExp> m_map_gExp; //gid gExp
    vector<short> m_vecborder;
};

struct cexp
{
    cexp(uint32_t cid, uint16_t mid, uint16_t exon):
    cellid(cid), midcnt(mid),exoncnt(exon){}
    uint16_t midcnt;
    uint16_t exoncnt;
    uint32_t cellid;
};

class geneUnit
{
public:
    geneUnit(){};
    ~geneUnit(){};
    void add(uint32_t cid, uint16_t mid, uint16_t exon)
    {
        m_vec_cexp.emplace_back(cid, mid, exon);
        m_expcnt += mid;
        m_exoncnt += exon;
        m_maxmid = std::max(m_maxmid, mid);
    }
public:
    vector<cexp> m_vec_cexp;
    uint16_t m_expcnt = 0;
    uint16_t m_exoncnt = 0;
    uint16_t m_maxmid = 0;
};


#endif