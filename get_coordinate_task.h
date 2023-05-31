/*
 * @Author: wanruiwen
 * @Date: 2023-03-3 14:53:03
 * @LastEditors: wanruiwen
 * @LastEditTime: 2023-03-3 14:53:03
 * @Description: get coordinate from gem file
 */

#ifndef GENETOH5_GETCOORDINATETASK_H
#define GENETOH5_GETCOORDINATETASK_H

#include <unordered_map>

#include "gef.h"
#include "thread_pool.h"
#include "utils.h"

// using namespace std;

typedef struct {
  int readlen;  // want read length
  int reallen;  // real read length
} RLen_;

class GetCoordinateTask : public ITask {
 public:
  GetCoordinateTask(
      gzFile file, int file_column, std::vector<int> &image_range,
      std::vector<CoordinateInfo> &coordinate_info);
  ~GetCoordinateTask();
  void doTask();

 private:
  void readbuf(RLen_ &rlen);
  int cuttail(char *pbuf);
  int GetCoordinateInfo();
  int MergeCoordinateinfo();

 private:
  int m_buflen = 0;
  int min_x = INT_MAX, min_y = INT_MAX, max_x = 0, max_y = 0;
  char *m_pbuf = nullptr;

  gzFile m_file;
  int file_column_;
  std::vector<int> &m_range;
  std::vector<CoordinateInfo> &coordinate_info_;
  std::vector<CoordinateInfo> coordinate_;
  static string m_leftstr;
  static mutex m_readmtx;   // read file lock
  static mutex m_mergemtx;  // merge lock
};

#endif  // GENETOH5_READTASK_H