
#ifndef GENETOH5_GENEINFOQUEUE_H
#define GENETOH5_GENEINFOQUEUE_H

#include <thread>
#include <mutex>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <condition_variable>
#include <iostream>
#include "gef.h"

using namespace std;

class GeneInfoQueue
{
  public:
    GeneInfoQueue(){};
    ~GeneInfoQueue(){};
    void addqueue(GeneS* ptr)
    {
        std::lock_guard<std::mutex> tlock(m_mtx_queue);
        m_vec_geneinfo.emplace_back(ptr);
        m_cv_queue.notify_all();
    }
    GeneS* getGeneInfo(unsigned int idx)
    {
        GeneS* ptr = nullptr;
        std::unique_lock<std::mutex> tlock(m_mtx_queue);

        m_cv_queue.wait(tlock, [idx,this] {return !m_vec_geneinfo.empty() && idx <= (m_vec_geneinfo.size()-1);});
        
        ptr = m_vec_geneinfo[idx];
        return ptr;
    }
    int getsize()
    {
        return m_vec_geneinfo.size();
    }

    void init(unsigned long size)
    {
        m_vec_geneinfo.reserve(size);
    }

    void clear(int bin)
    {
        if(bin == 1)
        {
            for(GeneS* ptr : m_vec_geneinfo)
            {
                delete ptr;
            }
        }
        else
        {
            for(GeneS* ptr : m_vec_geneinfo)
            {
                delete ptr->vecptr;
                delete ptr;
            }
        }
        m_vec_geneinfo.clear();
    }
  private:
    std::mutex m_mtx_queue; //队列锁
    std::condition_variable m_cv_queue;
    std::vector<GeneS*> m_vec_geneinfo;
};

#endif