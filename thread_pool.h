#ifndef GENETOH5_ThreadPool_H
#define GENETOH5_ThreadPool_H

#include <atomic>
#include <queue>
#include <vector>
#include <condition_variable>
#include <thread>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <unistd.h>
#endif


using namespace std;

class ITask
{
public:
    ITask(){};
    virtual ~ITask(){};
    virtual void doTask() = 0;
};


class ThreadPool
{
public:
    ThreadPool(int num);
    ~ThreadPool();
    int addTask(ITask* task);
    int addThread(int sz);
    int idlCount();
    int threadCount();
    void waitTaskDone();
    void setdnbThreadnum(int num)
    {
        m_dnbnum = num;
    }
private:
    int m_dnbnum = 1;
    vector<thread>   m_vecThread;            //线程池
    queue<ITask*>      m_tasks;           //任务队列
    mutex              m_mtx;            //同步
    condition_variable m_task_cv;         //条件阻塞
    atomic<bool>     m_run{ true };     //线程池是否执行
    atomic<int>      m_idlCnt{ 0 };  //空闲线程数量
};

#endif //GENETOH5_ThreadPool_H