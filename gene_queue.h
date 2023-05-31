
#ifndef GENETOH5_GENEQUEUE_H
#define GENETOH5_GENEQUEUE_H

#include <thread>
#include <mutex>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <queue>
#include <iostream>
#include <condition_variable>
#include "gef.h"

template<typename T>
class GefQueue
{
public:
    GefQueue(){};
    ~GefQueue(){};
    void addqueue(T *ptr)
    {
        std::lock_guard<std::mutex> tlock(m_mtx_queue);
        m_qgeneptr.emplace(ptr);
        m_cv_queue.notify_one();
    }

    T* getPtr()
    {
        T *ptr = nullptr;
        std::unique_lock<std::mutex> tlock(m_mtx_queue);

        m_cv_queue.wait(tlock, [this] {return !m_qgeneptr.empty();});
        ptr = m_qgeneptr.front();
        m_qgeneptr.pop();
        return ptr;
    }
private:
    std::mutex m_mtx_queue; //队列锁
    std::condition_variable m_cv_queue;
    std::queue<T *> m_qgeneptr;
};


#endif