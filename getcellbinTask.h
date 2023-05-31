/*
 * @Author: zhaozijian
 * @Date: 2022-05-05 17:30:29
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-10 15:32:54
 * @Description: file content
 */
#ifndef GEFTOOL_GETCELLBINTASK_H
#define GEFTOOL_GETCELLBINTASK_H

#include "thread_pool.h"
#include "gef.h"
#include <map>
#include <vector>
#include "cgefCellgem.h"
#include "gene_queue.h"
#include "opencv2/opencv.hpp"
//using namespace cv;


class getcellbinTask:public ITask
{
public:
    getcellbinTask(cgefCellgem *ptr, uint32_t label, cv::Rect &rect, vector<cv::Point> &vecpoint):
    m_cgefPtr(ptr),m_label(label),m_rect(rect),m_vecpoint(vecpoint)
    {
    }
    ~getcellbinTask()
    {

    }
    void doTask()
    {
        int cx = m_cgefPtr->m_centroids.at<double>(m_label,0);
        int cy = m_cgefPtr->m_centroids.at<double>(m_label,1);
        int area = m_cgefPtr->m_stats.at<int>(m_label, cv::CC_STAT_AREA);

        cellUnit *cptr = new cellUnit(cx, cy, area, m_label, m_cgefPtr->m_block_size);
        uint64_t l_id = 0;
        for(int j=m_rect.y;j<m_rect.y+m_rect.height;j++)
        {
            for(int i=m_rect.x;i<m_rect.x+m_rect.width;i++)
            {
                if(m_cgefPtr->m_outimg.at<unsigned int>(j,i) == m_label)
                {
                    l_id = i;
                    l_id = (l_id<<32) | j;
                    auto itor = m_cgefPtr->m_hash_vecdnb.find(l_id);
                    if(itor != m_cgefPtr->m_hash_vecdnb.end())
                    {
                        cptr->add(itor->second);
                    }
                }
            }
        }

        if(cptr->m_dnbcnt)
        {
            getborder(cptr);
        }

        m_cgefPtr->m_cellqueuePtr->addqueue(cptr);
    }

    void approx(vector<cv::Point> &in, vector<cv::Point> &out)
    {
        m_approxCount++;
        double epsilon;
        if (m_approxCount > 5) {
            epsilon = (1.0 * (m_approxCount / 5)) * arcLength(in, true);
        } else {
            epsilon = 0.01 * arcLength(in, true);
        }
        
        approxPolyDP(in, out, epsilon, true);
        if(out.size() > BORDERCNT)
        {
            in.swap(out);
            out.clear();
            approx(in, out);
        }
    }

    void getborder(cellUnit *cptr)
    {
        cptr->m_vecborder.reserve(BORDERCNT*2);
        int i=0;
        // int sz = m_vecpoint.size();
        // if(sz > BORDERCNT)
        // {
            vector<cv::Point> tmpborder;
            approx(m_vecpoint, tmpborder);
            m_approxCount = 0;
            int sz = tmpborder.size();
            assert(sz <= BORDERCNT);
            for(;i<sz;i++)
            {
                cptr->m_vecborder.emplace_back(tmpborder[i].x - cptr->m_cx);
                cptr->m_vecborder.emplace_back(tmpborder[i].y - cptr->m_cy);
            }
        // }
        // else
        // {
        //     for(;i<sz;i++)
        //     {
        //         cptr->m_vecborder.emplace_back(m_vecpoint[i].x - cptr->m_cx);
        //         cptr->m_vecborder.emplace_back(m_vecpoint[i].y - cptr->m_cy);
        //     }
        // }
        
        for(;i<BORDERCNT;i++) //不足补SHRT_MAX
        {
            cptr->m_vecborder.emplace_back(SHRT_MAX);
            cptr->m_vecborder.emplace_back(SHRT_MAX);
        }
    }

    // void getborder(cellUnit *cptr)
    // {
    //     cptr->m_vecborder.reserve(BORDERCNT*2);
    //     int i=0;
    //     int sz = m_vecpoint.size();

    //     vector<Point> tmpborder;
    //     approx(m_vecpoint, tmpborder);
    //     sz = tmpborder.size();
    //     for(;i<sz;i++)
    //     {
    //         cptr->m_vecborder.emplace_back(tmpborder[i].x - cptr->m_cx);
    //         cptr->m_vecborder.emplace_back(tmpborder[i].y - cptr->m_cy);
    //     }

    //     for(;i<BORDERCNT;i++) //不足补SHRT_MAX
    //     {
    //         cptr->m_vecborder.emplace_back(SHRT_MAX);
    //         cptr->m_vecborder.emplace_back(SHRT_MAX);
    //     }
    // }

private:
    uint32_t m_label = 0;
    uint32_t m_approxCount = 0;
    cv::Rect m_rect;
    vector<cv::Point> &m_vecpoint;
    cgefCellgem *m_cgefPtr = nullptr;
};

#endif