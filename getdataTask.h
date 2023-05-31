/*
 * @Author: zhaozijian
 * @Date: 2022-02-21 18:48:17
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-02-22 10:12:59
 * @Description: file content
 */
#ifndef GEFTOOL_GETDATATASK_H
#define GEFTOOL_GETDATATASK_H

#include "thread_pool.h"
#include "gef.h"

class getdataTask:public ITask
{
public:
    getdataTask(unsigned short gene_id, Gene * gene, Expression * expression, 
                        unordered_map<string, vector<Expression>> &hashexp):
    m_gene_id(gene_id),m_genePtr(gene),m_expPtr(expression),m_hashExp(hashexp)
    {

    }
    ~getdataTask(){}
    void setRange(unsigned int min_x,unsigned int min_y, unsigned int max_x, unsigned int max_y)
    {
        m_min_x = min_x;
        m_min_y = min_y;
        m_max_x = max_x;
        m_max_y = max_y;
    }
    void doTask()
    {
        vector<Expression> exps;
        exps.reserve(m_genePtr[m_gene_id].count);
        unsigned int end = m_genePtr[m_gene_id].offset + m_genePtr[m_gene_id].count;
        for(unsigned int i = m_genePtr[m_gene_id].offset; i < end; i++){
            Expression &exp = m_expPtr[i];
            if(exp.x < m_min_x || exp.x > m_max_x || exp.y < m_min_y || exp.y > m_max_y){
                continue;
            }

            exps.emplace_back(exp);
        }

        {
            std::lock_guard<std::mutex> tlock(m_mtx);
            string str(m_genePtr[m_gene_id].gene);
            m_hashExp.emplace(str, std::move(exps));
        }
    }
private:
    unsigned short m_gene_id;
    unsigned int m_min_x, m_min_y, m_max_x, m_max_y;
    Gene *m_genePtr = nullptr;
    Expression *m_expPtr = nullptr;
    unordered_map<string, vector<Expression>> &m_hashExp;
    static std::mutex m_mtx;
};



#endif