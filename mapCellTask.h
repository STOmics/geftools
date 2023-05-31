/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 16:56:28
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-02 16:13:11
 * @Description: file content
 */
#ifndef GEFTOOLS_MAPCELL_H_
#define GEFTOOLS_MAPCELL_H_

#include "thread_pool.h"
#include "gef.h"
#include "cgefParam.h"
#include <set>
#include "opencv2/opencv.hpp"
#include "polygon.h"

class mapCellTask:public ITask
{
public:
    mapCellTask(int id, const vector<GefTools::Polygon> &vec_poly):
    m_id(id), m_vec_poly(vec_poly)
    {
        m_cparamPtr = cgefParam::GetInstance();
        int len = vec_poly.size() / m_cparamPtr->m_threadcnt + 1;
        start = id*len;
        end = start+len;
        if(id == m_cparamPtr->m_threadcnt-1)
        {
            end = vec_poly.size();
        }
        printf("%d %d %d\n", id, start, end);
    }
    ~mapCellTask(){};
    void doTask()
    {
        std::set<int> set_clabel;
        double ret;
        int x, y;
        bool bfind = false;
        for(int i=start;i<end;i++)
        {
            auto itor = m_cparamPtr->m_map_cell.begin();
            for(;itor != m_cparamPtr->m_map_cell.end();itor++)
            {
                if(set_clabel.find(itor->first) != set_clabel.end())
                {
                    continue;
                }
                bfind = false;
                cgef_cell *cellptr = itor->second;
                x = cellptr->m_xtotal / cellptr->dnbcnt;
                y = cellptr->m_ytotal / cellptr->dnbcnt;

                ret = cv::pointPolygonTest(m_vec_poly[i].getBorder(), cv::Point2f(x, y), false);
                if(ret >= 0)
                {
                    bfind = true;
                    set_clabel.insert(itor->first);
                    m_cparamPtr->m_vec_clabel[i] = itor->first;
                    break;
                }
            }

            // if(!bfind)
            // {
            //     printf("err %d \n", i);
            // }
        }
    }
private:
    int m_id;
    int start,end;
    cgefParam *m_cparamPtr;
    const vector<GefTools::Polygon> &m_vec_poly;
};


#endif