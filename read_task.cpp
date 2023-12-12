#include "read_task.h"

#include <zlib.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>

#include "bgef_options.h"

string ReadTask::m_leftstr;
mutex ReadTask::m_readmtx;
mutex ReadTask::m_mergemtx;

ReadTask::ReadTask(bool bexon, gzFile file, std::vector<int> &range,
                   std::unordered_map<std::string, std::vector<Expression>> &map_gene_exp) :
    m_bexon(bexon), m_file(file), m_range(range), m_map_gene_exp(map_gene_exp) {
    m_bexon = bexon;
    m_pbuf = new char[READLEN];
}

ReadTask::~ReadTask() { delete[] m_pbuf; }

void ReadTask::doTask() {
    RLen rlen {0, 0};
    std::function<int(ReadTask *)> parse;

    if (BgefOptions::GetInstance()->has_gene_name) {
        if (m_bexon) {
            parse = &ReadTask::getGeneInfoWithGenename_exon;
        } else {
            parse = &ReadTask::getGeneInfoWithGenename;
        }
    } else {
        if (m_bexon) {
            parse = &ReadTask::getGeneInfo_exon;
        } else {
            parse = &ReadTask::getGeneInfo;
        }
    }

    while (true) {
        readbuf(rlen);
        parse(this);
        if (rlen.reallen < rlen.readlen)  // file end
        {
            break;
        }
    }

    mergeGeneinfo();
}

void ReadTask::readbuf(RLen &rlen) {
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
        log_error << errorCode::E_PARSEFILEERROR << "read error " << errmsg;
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

int ReadTask::cuttail(char *pbuf) {
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

int ReadTask::getGeneInfo() {
    int i = 0, k = 0;
    char *ptr = m_pbuf;
    std::string geneid;
    Expression expression {0, 0, 0};
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    geneid.clear();
                    geneid.append(ptr, &m_pbuf[i] - ptr);
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
                    // printf("k %d\n", k);
                    expression.count = atoi(ptr);
                    k = 0;
                    ptr = &m_pbuf[i + 1];
                    // printf("%s %d %d %d\n", geneid.c_str(), gene.x, gene.y, gene.cnt);
                    m_map_gege[geneid].push_back(expression);
                    break;
                default:
                    break;
            }
        }
    }

    return m_map_gege.size();
}

int ReadTask::mergeGeneinfo() {
    lock_guard<mutex> lock(m_mergemtx);

    m_range[0] = std::min(m_range[0], min_x);
    m_range[1] = std::max(m_range[1], max_x);
    m_range[2] = std::min(m_range[2], min_y);
    m_range[3] = std::max(m_range[3], max_y);

    auto itor = m_map_gege.begin();
    for (; itor != m_map_gege.end(); itor++) {
        std::vector<Expression> &vec = m_map_gene_exp[itor->first];
        vec.insert(vec.end(), itor->second.begin(), itor->second.end());
    }

    if (BgefOptions::GetInstance()->has_gene_name) {
        BgefOptions::GetInstance()->gene_name_map.insert(m_genename_map_per_t.begin(), m_genename_map_per_t.end());
    }

    return 0;
}

int ReadTask::getGeneInfo_exon() {
    int i = 0, k = 0;
    char *ptr = m_pbuf;
    std::string geneid;
    Expression expression {0, 0, 0, 0};
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    geneid.clear();
                    geneid.append(ptr, &m_pbuf[i] - ptr);
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
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 4:
                    expression.exon = atoi(ptr);
                    k = 0;
                    ptr = &m_pbuf[i + 1];
                    m_map_gege[geneid].push_back(expression);
                    break;
                default:
                    break;
            }
        }
    }

    return m_map_gege.size();
}

int ReadTask::getGeneInfoWithGenename_exon() {
    int i = 0, k = 0;
    char *ptr = m_pbuf;
    std::string geneid;
    std::string genename;
    Expression expression {0, 0, 0, 0};
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    geneid.clear();
                    geneid.append(ptr, &m_pbuf[i] - ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 1:
                    genename.clear();
                    genename.append(ptr, &m_pbuf[i] - ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 2:
                    expression.x = atoi(ptr);
                    min_x = std::min(expression.x, min_x);
                    max_x = std::max(expression.x, max_x);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 3:
                    expression.y = atoi(ptr);
                    min_y = std::min(expression.y, min_y);
                    max_y = std::max(expression.y, max_y);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 4:
                    expression.count = atoi(ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 5:
                    expression.exon = atoi(ptr);
                    k = 0;
                    ptr = &m_pbuf[i + 1];
                    m_map_gege[geneid].push_back(expression);
                    m_genename_map_per_t[geneid] = genename;
                    break;
                default:
                    break;
            }
        }
    }

    return m_map_gege.size();
}

int ReadTask::getGeneInfoWithGenename() {
    int i = 0, k = 0;
    char *ptr = m_pbuf;
    std::string geneid;
    std::string genename;
    Expression expression {0, 0, 0, 0};
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    geneid.clear();
                    geneid.append(ptr, &m_pbuf[i] - ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 1:
                    genename.clear();
                    genename.append(ptr, &m_pbuf[i] - ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 2:
                    expression.x = atoi(ptr);
                    min_x = std::min(expression.x, min_x);
                    max_x = std::max(expression.x, max_x);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 3:
                    expression.y = atoi(ptr);
                    min_y = std::min(expression.y, min_y);
                    max_y = std::max(expression.y, max_y);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 4:
                    expression.count = atoi(ptr);
                    k = 0;
                    ptr = &m_pbuf[i + 1];
                    m_map_gege[geneid].push_back(expression);
                    m_genename_map_per_t[geneid] = genename;
                    break;
                default:
                    break;
            }
        }
    }

    return m_map_gege.size();
}