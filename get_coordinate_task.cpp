#include <zlib.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "get_coordinate_task.h"

string GetCoordinateTask::m_leftstr;
mutex GetCoordinateTask::m_readmtx;
mutex GetCoordinateTask::m_mergemtx;

GetCoordinateTask::GetCoordinateTask(
    gzFile file, int file_column, std::vector<int> &image_range,
    std::vector<CoordinateInfo> &coordinate_info)
    : m_file(file),
      file_column_(file_column), 
      m_range(image_range),
      coordinate_info_(coordinate_info) {
  m_pbuf = new char[READLEN];
}

GetCoordinateTask::~GetCoordinateTask() { delete[] m_pbuf; }

void GetCoordinateTask::doTask() {
  RLen_ rlen{0, 0};
  while (true) {
    readbuf(rlen);
    GetCoordinateInfo();
    if (rlen.reallen < rlen.readlen)
    {
      break;
    }
  }
  MergeCoordinateinfo();
}

void GetCoordinateTask::readbuf(RLen_ &rlen) {
  lock_guard<mutex> lock(m_readmtx);
  char *pbuf = m_pbuf;
  int leftlen = m_leftstr.length();
  memcpy(pbuf, m_leftstr.c_str(), leftlen);
  m_leftstr.clear();
  pbuf += leftlen;
  rlen.readlen = READLEN - leftlen;
  rlen.reallen = gzread(m_file, pbuf, rlen.readlen);

  if (rlen.reallen == -1)  // error
  {
    int z_errnum = 0;
    const char *errmsg = gzerror(m_file, &z_errnum);
    if (z_errnum == Z_ERRNO) errmsg = strerror(errno);
    printf("read error %s", errmsg);
    char errMsg2File[32] = {0};
    sprintf(errMsg2File, "read error %s", errmsg);
    reportErrorCode2File(errorCode::E_PARSEFILEERROR, errMsg2File);
    exit(1);
  }

  m_buflen = rlen.reallen;
  if (rlen.reallen == rlen.readlen) {
    cuttail(m_pbuf);
  } else {
    if (m_buflen != 0) m_buflen += leftlen;
    // printf("m_buflen: %d\n", m_buflen);
  }
}

int GetCoordinateTask::cuttail(char *pbuf) {
  int i = READLEN - 1;
  for (; i > 0; i--) {
    if (pbuf[i] == '\n') {
      break;
    }
  }

  m_buflen = i + 1;
  m_leftstr.append(&pbuf[m_buflen], READLEN - m_buflen);
  return 0;
}

int GetCoordinateTask::GetCoordinateInfo() {
  int i = 0, k = 0;
  char *ptr = m_pbuf;
  CoordinateInfo expression{0, 0, 0};
  for (; i < m_buflen; i++) {
    if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' ||
        m_pbuf[i] == '\n') {
      switch (k) {
        case 0:
          k++;
          ptr = &m_pbuf[i + 1];
          break;
        case 1:
          expression.x = atoi(ptr);
          min_x = std::min(expression.x, min_x);
          max_x = std::max(expression.x, max_x);
          k++;
          ptr = &m_pbuf[i + 1];
          break;
        case 2:
          expression.y = atoi(ptr);
          min_y = std::min(expression.y, min_y);
          max_y = std::max(expression.y, max_y);
          k++;
          ptr = &m_pbuf[i + 1];
          break;
        case 3:
          expression.count = atoi(ptr);
          ptr = &m_pbuf[i + 1];
          coordinate_.emplace_back(expression);
          if ((k + 1) == file_column_) {
            k = 0;
          } else {
            k++;
          }
          break;
        default:
          if ((k + 1) == file_column_) {
            k = 0;
            ptr = &m_pbuf[i + 1];
          } else {
            k++;
            ptr = &m_pbuf[i + 1];
          }
          break;
      }
    }
  }

  return coordinate_.size();
}

int GetCoordinateTask::MergeCoordinateinfo() {
  lock_guard<mutex> lock(m_mergemtx);

  m_range[0] = std::min(m_range[0], min_x);
  m_range[1] = std::max(m_range[1], max_x);
  m_range[2] = std::min(m_range[2], min_y);
  m_range[3] = std::max(m_range[3], max_y);

  coordinate_info_.insert(coordinate_info_.end(), coordinate_.begin(),
                          coordinate_.end());

  return 0;
}
