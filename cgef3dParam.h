#ifndef GEFTOOLS_CGEF3DPARAM_H_
#define GEFTOOLS_CGEF3DPARAM_H_

#include <zlib.h>
#include <map>
#include <unordered_map>
#include <set>
#include "opencv2/opencv.hpp"
#include <float.h>

class GEFTOOLS_API cgef3d_cell
{
public:
    cgef3d_cell(){}
    ~cgef3d_cell(){}
    void add(float tx, float ty, float umi)
    {
        m_sumumi += umi;
        m_dnbcnt++;
        m_vecPoint.emplace_back(tx,ty);
    }
    void merge(cgef3d_cell &other)
    {
        m_vecPoint.insert(m_vecPoint.end(), other.m_vecPoint.begin(), other.m_vecPoint.end());
        m_dnbcnt += other.m_dnbcnt;
        m_sumumi += other.m_sumumi;
    }
    void setCellInfo(float tx, float ty, uint16_t area, vector<cv::Point> &vec)
    {
        binfo = true;
        m_x = tx;
        m_y = ty;
        m_area = area;
        m_border.insert(m_border.end(), vec.begin(), vec.end());
    }
    bool getCellInfo()
    {
        if(binfo) return true;
        int vsz = m_vecPoint.size();
        //if(vsz<3) return false;
        convexHull(m_vecPoint, m_border, true);
        int sz = m_border.size();
        if(sz <= 2)
        {
            if(m_vecPoint[0] == m_vecPoint[vsz-1]) 
            {
                m_x = m_vecPoint[0].x;
                m_y = m_vecPoint[0].y;
            }
            else
            {
                vector<float> vecx; vecx.reserve(sz);
                vector<float> vecy; vecy.reserve(sz);
                for(cv::Point2f &pt : m_vecPoint)
                {
                    vecx.emplace_back(pt.x);
                    vecy.emplace_back(pt.y);
                }

                std::sort(vecx.begin(), vecx.end(), [](float a, float b){return a<b;});
                std::sort(vecy.begin(), vecy.end(), [](float a, float b){return a<b;});

                int pos = ceil( (sz + 1)*1.0 /2 );
                m_x = ceil (vecx[pos-2]*0.5+vecx[pos-1]*0.5);
                m_y = ceil (vecy[pos-2]*0.5+vecy[pos-1]*0.5);
            }
        }
        else
        {
            // if(sz > BORDERCNT)
            // {
            //     double epsilon = 0.01 * arcLength(hull, true);
            //     approxPolyDP(hull, m_border, epsilon, true); 
            // }
            // else
            // {
            //     m_border.swap(hull);
            // }

            cv::Moments mu = cv::moments(m_border, true);
            if(mu.m00 == 0) return false;
            m_x = mu.m10/mu.m00;
            m_y = mu.m01/mu.m00;
            m_area = mu.m00;
        }
        return true;
    }
public:
    bool binfo = false;
    uint16_t m_dnbcnt = 0;
    uint16_t m_area = 0;
    float m_sumumi = 0.0;
    float m_x = 0.0;
    float m_y = 0.0;
    vector<cv::Point2f> m_vecPoint;
    vector<cv::Point2f> m_border;
};

class GEFTOOLS_API cgef3d_gene
{
public:
    cgef3d_gene(){};
    ~cgef3d_gene(){};
    void add(uint32_t cid, float umi)
    {
        if(m_map_cell.find(cid) == m_map_cell.end())
        {
            m_map_cell.emplace(cid, umi);
        }
        else
        {
            m_map_cell[cid] += umi;
        }
        m_sumumi += umi;
    }

    void merge(cgef3d_gene &other)
    {
        std::map<uint32_t, float> & other_map = other.m_map_cell;
        auto itor = other_map.begin();
        for(;itor != other_map.end();itor++)
        {
            if(m_map_cell.find(itor->first) != m_map_cell.end())
            {
                m_map_cell.emplace(itor->first, itor->second);
            }
            else
            {
                m_map_cell[itor->first] += itor->second;
            }
        }
        m_sumumi += other.m_sumumi;
    }
public:
    float m_sumumi = 0.0;
    std::map<uint32_t, float> m_map_cell;
};

class GEFTOOLS_API cgef3d_ctype
{
public:
    cgef3d_ctype(){}
    ~cgef3d_ctype(){}
    void add(uint32_t cid)
    {
        m_vec_cid.push_back(cid);
    }
    void merge(cgef3d_ctype &other)
    {
        // std::set<uint32_t> &tmpset = other.m_set_cid;
        // auto itor = tmpset.begin();
        // for(;itor != tmpset.end();itor++)
        // {
        //     m_set_cid.insert(*itor);
        // }
        m_vec_cid.insert(m_vec_cid.end(), other.m_vec_cid.begin(), other.m_vec_cid.end());
    }

public:
    std::vector<uint32_t> m_vec_cid;
};

class GEFTOOLS_API cgef3dParam
{

public:
    static cgef3dParam *GetInstance()
    {
        static cgef3dParam instance;
        return &instance;
    }

private:
    cgef3dParam(/* args */){};
    ~cgef3dParam(){};
public:
    gzFile m_infile;
    int m_threadcnt = 8;
    unordered_map<uint32_t, cgef3d_cell*> m_map_cell;
    unordered_map<string, cgef3d_gene*> m_map_gene;
    //unordered_map<string, cgef3d_ctype*> m_map_ctype;
};

#endif