#ifndef GEFTOOLS_GETSAPDATATASK_H
#define GEFTOOLS_GETSAPDATATASK_H

#include "gef.h"
#include "thread_pool.h"
#include "opencv2/opencv.hpp"
//using namespace cv;

class getsapdataTask:public ITask
{
public:
    getsapdataTask(int i, int tcnt, cv::Mat &mat, BinStat *ptr, vector<sapBgefData> &vecdata):
    m_mat(mat),m_ptr(ptr),m_vecdata(vecdata)
    {
        int cnt = m_mat.rows/tcnt + 1;
        m_start = i*cnt;
        if(i == tcnt -1)
        {
            m_end = m_mat.rows;
        }
        else
        {
            m_end = m_start+cnt;
        }
    }
    ~getsapdataTask(){}
    void doTask()
    {
        vector<sapBgefData> vettmp;
        int id = 0;
        for(uint32_t j=0;j<m_mat.cols;j++)
        {
            for(uint32_t i=m_start;i<m_end;i++)
            {
                id = j*m_mat.rows+i;
                if(m_mat.at<uchar>(i,j))
                {
                    if(m_ptr[id].gene_count)
                    {
                        vettmp.emplace_back(m_ptr[id].gene_count, m_ptr[id].mid_count, j, i);
                    }
                }
            }
        }

        {
            std::lock_guard<std::mutex> tlock(m_mtx);
            m_vecdata.insert(m_vecdata.end(), vettmp.begin(), vettmp.end());
        }
    }
private:
    int m_start, m_end;
    cv::Mat &m_mat;
    BinStat *m_ptr;
    vector<sapBgefData> &m_vecdata;
    static std::mutex m_mtx;
};

#endif