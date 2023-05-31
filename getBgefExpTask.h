/*
 * @Author: zhaozijian
 * @Date: 2022-05-13 09:36:11
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-13 10:43:41
 * @Description: file content
 */
#ifndef GEFTOOLS_GETBGEFEXPTASK_H
#define GEFTOOLS_GETBGEFEXPTASK_H


#include "thread_pool.h"
#include "gef.h"

class getBgefExpTask:public ITask
{
public:
    getBgefExpTask(uint32_t count, Expression* expData, unsigned int *pcount, unsigned long long *pcellid):
    m_cnt(count), m_pexpData(expData), m_pcount(pcount), m_pcellid(pcellid)
    {

    }
    ~getBgefExpTask(){}
    void doTask()
    {
        unsigned long long uniq_cell_id;
        for(uint64_t i=0;i<m_cnt;i++)
        {
            uniq_cell_id = m_pexpData[i].x;
            uniq_cell_id = (uniq_cell_id << 32) | m_pexpData[i].y;
            m_pcount[i] = m_pexpData[i].count;
            m_pcellid[i] = uniq_cell_id;
        }
    }
private:
    uint64_t m_cnt;
    Expression* m_pexpData = nullptr;
    unsigned int *m_pcount = nullptr;
    unsigned long long *m_pcellid = nullptr;
};


#endif