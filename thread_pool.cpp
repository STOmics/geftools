#include "thread_pool.h"

ThreadPool::ThreadPool(int num)
{
    addThread(num);
}

ThreadPool::~ThreadPool()
{
    m_run = false;
    m_task_cv.notify_all();  // 唤醒所有线程执行
    for (thread& thread : m_vecThread)
    {
        // thread.detach(); // 让线程“自生自灭”
        if (thread.joinable())
            thread.join();  // 等待任务结束， 前提：线程一定会执行完
    }
}

int ThreadPool::addTask(ITask* task)
{
    lock_guard<mutex> lock(m_mtx);
    m_tasks.emplace(task);
    int cnt = m_tasks.size();
    m_task_cv.notify_one();
    return cnt;
}

int ThreadPool::addThread(int size)
{
    for(;m_vecThread.size() < 128 && size > 0; --size)
    {
        m_vecThread.emplace_back([this] {  //工作线程函数
            while (m_run)
            {
                ITask *ptask;  // 获取一个待执行的 task
                {
                    // unique_lock 相比 lock_guard 的好处是：可以随时 unlock() 和 lock()
                    unique_lock<mutex> lock{m_mtx};
                    m_task_cv.wait(lock, [this] { return !m_run || !m_tasks.empty(); });  // wait 直到有 task
                    if (!m_run && m_tasks.empty())
                        return;
                    ptask = m_tasks.front();  // 按先进先出从队列取一个 task
                    m_tasks.pop();
                    m_idlCnt--;
                }
                
                ptask->doTask();  //执行任务
                delete ptask;
                m_idlCnt++;
            }
        });
        m_idlCnt++;
    }
    return 0;
}

//空闲线程数量
int ThreadPool::idlCount()
{
    return m_idlCnt;
}
//线程数量
int ThreadPool::threadCount()
{
    return m_vecThread.size();
}

void ThreadPool::waitTaskDone()
{
    while (true)
    {
        if(m_idlCnt == m_vecThread.size() && m_tasks.empty())
        {
            break;
        }
#ifdef _WIN32
        Sleep(1);
#else
        sleep(1);
#endif
    }
    
}
