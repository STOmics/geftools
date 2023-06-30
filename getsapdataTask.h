#ifndef GEFTOOLS_GETSAPDATATASK_H
#define GEFTOOLS_GETSAPDATATASK_H

#include "gef.h"
#include "opencv2/opencv.hpp"
#include "thread_pool.h"
// using namespace cv;

class getsapdataTask : public ITask {
  public:
    getsapdataTask(int i, int tcnt, cv::Mat &mat, BinStat *ptr, vector<sapBgefData> &vecdata) :
        m_mat(mat), m_ptr(ptr), m_vecdata(vecdata) {
        int cnt = m_mat.rows / tcnt + 1;
        m_start = i * cnt;
        if (i == tcnt - 1) {
            m_end = m_mat.rows;
        } else {
            m_end = m_start + cnt;
        }
    }
    ~getsapdataTask() {}

    void doTask() {
        vector<sapBgefData> vettmp;
        int id = 0;
        for (uint32_t j = 0; j < m_mat.cols; j++) {
            for (uint32_t i = m_start; i < m_end; i++) {
                id = j * m_mat.rows + i;
                if (m_mat.at<uchar>(i, j)) {
                    if (m_ptr[id].gene_count) {
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

class getLabelInfoTask : public ITask {
  public:
    getLabelInfoTask(int i, int tcnt, cv::Mat &mat, Gene *ptr, Expression *e_ptr, std::vector<LabelGeneData> &vecdata,
                     uint32_t &total_mid, uint32_t &gene_cnt) :
        m_mat(mat), m_ptr(ptr), e_ptr_(e_ptr), m_vecdata(vecdata), total_mid_(total_mid), gene_cnt_(gene_cnt) {
        int cnt = gene_cnt / tcnt + 1;
        m_start = i * cnt;

        if (i == tcnt - 1) {
            m_end = gene_cnt;
        } else {
            m_end = m_start + cnt;
        }
    }
    ~getLabelInfoTask() {}

    void doTask() {
        std::vector<LabelGeneData> vettmp;
        int tmp_ttmidcnt;

        for (uint32_t i = m_start; i < m_end; i++) {
            Expression *ptr = e_ptr_ + m_ptr[i].offset;
            LabelGeneData tmp_ld("", 0);
            for (int j = 0; j < m_ptr[i].count; j++) {
                if (m_mat.at<uchar>(ptr[j].y, ptr[j].x)) {
                    strcpy(tmp_ld.gene_name, m_ptr[i].gene);
                    tmp_ld.mid_cnt += ptr[j].count;
                    tmp_ttmidcnt += ptr[j].count;
                }
            }
            if (tmp_ld.mid_cnt > 0) {
                vettmp.emplace_back(tmp_ld);
            }
        }

        {
            std::lock_guard<std::mutex> tlock(m_mtx_);
            if (!vettmp.empty()) {
                m_vecdata.insert(m_vecdata.end(), vettmp.begin(), vettmp.end());
                total_mid_ += tmp_ttmidcnt;
            }
        }
    }

  private:
    int m_start, m_end;
    cv::Mat &m_mat;
    Gene *m_ptr;
    Expression *e_ptr_;
    uint32_t &total_mid_;
    uint32_t &gene_cnt_;
    std::vector<LabelGeneData> &m_vecdata;
    static std::mutex m_mtx_;
};

#endif