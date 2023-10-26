#pragma once

#include "gef.h"
#include "thread_pool.h"

struct leveldnb {
    leveldnb() {}
    leveldnb(float tx, float ty, uint32_t tg, uint32_t tm, float tc) :
        x(tx), y(ty), genecnt(tg), midcnt(tm), color(tc) {}
    float x;
    float y;
    uint32_t genecnt;
    uint32_t midcnt;
    float color;
};

struct level_fractal {
    uint32_t fractal_size;
    uint32_t fractal_num;
    uint32_t mid_fractal_coor;
    uint32_t start_fractal_coor;
    uint32_t end_fractal_coor;
};

struct dnbbuf {
    uint32_t sz;
    uint32_t cnt;
    void *pbuf;
};

class getleveldnbtask : public ITask {
  public:
    getleveldnbtask(bool blevel, bool btop, int bin, uint32_t start, uint32_t len, uint32_t cols, uint32_t offset_x,
                    uint32_t offset_y, uint32_t maxmid, level_fractal &lf, dnbbuf &dbuf, uint64_t datacols,
                    vector<unsigned long long> &vecindex, BinStatUS *pdataus = nullptr, BinStat *pdata = nullptr) :
        blevel_(blevel),
        btop_(btop),
        bin_(bin),
        start_(start),
        len_(len),
        cols_(cols),
        offset_x_(offset_x),
        offset_y_(offset_y),
        maxmid_(maxmid),
        lf_(lf),
        dbuf_(dbuf),
        datacols_(datacols),
        vecindex_(vecindex),
        pdataus_(pdataus),
        pdata_(pdata) {}
    ~getleveldnbtask() {}

    bool leveltop(uint32_t x, uint32_t y) {
        uint32_t mx = x % lf_.fractal_num;
        bool xret = mx == lf_.start_fractal_coor | mx == lf_.mid_fractal_coor | mx == lf_.end_fractal_coor;
        if (xret) {
            uint32_t my = y % lf_.fractal_num;
            bool yret = my == lf_.start_fractal_coor | my == lf_.mid_fractal_coor | my == lf_.end_fractal_coor;
            if (yret) {
                return true;
            }
        }
        return false;
    }

    bool levelnormal(uint32_t x, uint32_t y) {
        uint32_t mx = x % lf_.fractal_num;
        bool xret = mx == lf_.start_fractal_coor | mx == lf_.mid_fractal_coor | mx == lf_.end_fractal_coor;
        if (xret) {
            uint32_t my = y % lf_.fractal_num;
            bool yret = my == lf_.start_fractal_coor | my == lf_.mid_fractal_coor | my == lf_.end_fractal_coor;
            if (yret) {
                if (!((mx == lf_.mid_fractal_coor) && (my == lf_.mid_fractal_coor)))  // 去掉中心点
                {
                    return true;
                }
            }
        }
        return false;
    }

    void doTask() {
        uint32_t x, y;
        uint64_t index = 0;
        if (pdataus_) {
            pdataus_ = pdataus_ + start_;
            if (btop_) {
                for (uint32_t i = 0; i < len_; i++) {
                    if (pdataus_[i].gene_count) {
                        x = ((i + start_) / cols_ + offset_x_) * bin_;
                        y = ((i + start_) % cols_ + offset_y_) * bin_;
                        bool ret = blevel_ ? leveltop(x, y) : true;
                        if (ret) {
                            m_vecdnb.emplace_back(x * bin_, y * bin_, pdataus_[i].mid_count, pdataus_[i].gene_count,
                                                  pdataus_[i].mid_count * 1.0 / maxmid_);
                            index = x;
                            index *= datacols_;
                            index += y;
                            vecindex_.push_back(index);
                        }

                    }
                }
            } else {
                for (uint32_t i = 0; i < len_; i++) {
                    if (pdataus_[i].gene_count) {
                        x = ((i + start_) / cols_ + offset_x_) * bin_;
                        y = ((i + start_) % cols_ + offset_y_) * bin_;
                        bool ret = blevel_ ? levelnormal(x, y) : true;
                        if (ret) {
                            m_vecdnb.emplace_back(x * bin_, y * bin_, pdataus_[i].mid_count, pdataus_[i].gene_count,
                                                  pdataus_[i].mid_count * 1.0 / maxmid_);
                            index = x;
                            index *= datacols_;
                            index += y;
                            vecindex_.push_back(index);
                        }

                    }
                }
            }
        } else if (pdata_) {
            pdata_ = pdata_ + start_;
            if (btop_) {
                for (uint32_t i = 0; i < len_; i++) {
                    if (pdata_[i].gene_count) {
                        x = (i + start_) / cols_ + offset_x_;
                        y = (i + start_) % cols_ + offset_y_;
                        bool ret = blevel_ ? leveltop(x, y) : true;
                        if (ret) {
                            m_vecdnb.emplace_back(x * bin_, y * bin_, pdata_[i].mid_count, pdata_[i].gene_count,
                                                  pdata_[i].mid_count * 1.0 / maxmid_);
                            index = x;
                            index *= datacols_;
                            index += y;
                            vecindex_.push_back(index);
                        }

                    }
                }
            } else {
                for (uint32_t i = 0; i < len_; i++) {
                    if (pdata_[i].gene_count) {
                        x = (i + start_) / cols_ + offset_x_;
                        y = (i + start_) % cols_ + offset_y_;
                        bool ret = blevel_ ? levelnormal(x, y) : true;
                        if (ret) {
                            m_vecdnb.emplace_back(x * bin_, y * bin_, pdata_[i].mid_count, pdata_[i].gene_count,
                                                  pdata_[i].mid_count * 1.0 / maxmid_);
                            index = x;
                            index *= datacols_;
                            index += y;
                            vecindex_.push_back(index);
                        }

                    }
                }
            }
        }

        {
            std::lock_guard<std::mutex> tlock(m_mtx);
            uint32_t sz = m_vecdnb.size() * sizeof(leveldnb);
        #ifdef _WIN32
            memcpy((char*)dbuf_.pbuf + dbuf_.sz, m_vecdnb.data(), sz);
        #else
            memcpy(dbuf_.pbuf + dbuf_.sz, m_vecdnb.data(), sz);
        #endif
            dbuf_.sz += sz;
            dbuf_.cnt += m_vecdnb.size();
        }
    }

  private:
    bool blevel_ = true, btop_ = false;  // 是否要进行分层过滤
    int bin_;
    uint32_t start_ = 0, len_ = 0, cols_ = 0, offset_x_ = 0, offset_y_ = 0, maxmid_ = 0;
    level_fractal &lf_;
    dnbbuf &dbuf_;
    BinStatUS *pdataus_ = nullptr;
    BinStat *pdata_ = nullptr;
    vector<leveldnb> m_vecdnb;
    static std::mutex m_mtx;
    uint64_t datacols_;
    vector<unsigned long long> &vecindex_;
};
