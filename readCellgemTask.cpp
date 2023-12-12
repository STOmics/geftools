/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:18:37
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-18 19:11:47
 * @Description: file content
 */

#include "readCellgemTask.h"

#include "cgefParam.h"

string readCellgemTask::m_leftstr;
mutex readCellgemTask::m_readmtx;
mutex readCellgemTask::m_mergemtx;

readCellgemTask::readCellgemTask() { m_pbuf = new char[READLEN]; }

readCellgemTask::~readCellgemTask() { delete[] m_pbuf; }

void readCellgemTask::doTask() {
    bool brun = true;
    while (brun) {
        brun = readbuf();
        getInfo();
    }

    mergeinfo();
    printf("read task end\n");
}

bool readCellgemTask::readbuf() {
    lock_guard<mutex> lock(m_readmtx);
    char *pbuf = m_pbuf;
    int leftlen = m_leftstr.length();
    memcpy(pbuf, m_leftstr.c_str(), leftlen);
    m_leftstr.clear();
    pbuf += leftlen;
    int readlen = READLEN - leftlen;
    int reallen = gzread(cgefParam::GetInstance()->m_infile, pbuf, readlen);

    m_buflen = reallen;
    if (reallen == readlen) {
        cuttail(m_pbuf);
    } else {
        if (m_buflen != 0) m_buflen += leftlen;
        return false;
    }
    return true;
}

int readCellgemTask::cuttail(char *pbuf) {
    int i = READLEN - 1;
    for (; i > 0; i--) {
        if (pbuf[i] == '\n') {
            break;
        }
    }
    // pbuf[i] = '\0';
    m_buflen = i + 1;
    m_leftstr.append(&pbuf[m_buflen], READLEN - m_buflen);
    return 0;
}

int readCellgemTask::mergeinfo() {
    lock_guard<mutex> lock(m_mergemtx);

    auto itor = m_map_cell.begin();
    auto &tmap_cell = cgefParam::GetInstance()->m_map_cell;
    for (; itor != m_map_cell.end(); itor++) {
        if (tmap_cell.find(itor->first) != tmap_cell.end()) {
            tmap_cell[itor->first]->merge(*(itor->second));
            delete itor->second;
        } else {
            tmap_cell.emplace(itor->first, itor->second);
        }
    }

    // auto itor_g = m_map_gene.begin();
    // auto &tmp_gene = cgefParam::GetInstance()->m_map_gene;
    // for(;itor_g!=m_map_gene.end();itor_g++)
    // {
    //     if(tmp_gene.find(itor_g->first) != tmp_gene.end())
    //     {
    //         tmp_gene[itor_g->first]->merge(*(itor_g->second));
    //         delete itor_g->second;
    //     }
    //     else
    //     {
    //         tmp_gene.emplace(itor_g->first, itor_g->second);
    //     }
    // }

    auto itor_g_i = m_map_gene_id.begin();
    auto &tmp_g_i = cgefParam::GetInstance()->m_map_gene_id;
    for (; itor_g_i != m_map_gene_id.end(); itor_g_i++) {
        if (tmp_g_i.find(itor_g_i->first) == tmp_g_i.end()) {
            tmp_g_i.emplace(itor_g_i->first, 0);
        }
    }

    if (cgefParam::GetInstance()->has_gene_name) {
        cgefParam::GetInstance()->cgem_genename_map.insert(m_genename_map_per_t.begin(), m_genename_map_per_t.end());
    }

    cgefParam::GetInstance()->m_min_x = std::min(cgefParam::GetInstance()->m_min_x, m_min_x);
    cgefParam::GetInstance()->m_max_x = std::max(cgefParam::GetInstance()->m_max_x, m_max_x);
    cgefParam::GetInstance()->m_min_y = std::min(cgefParam::GetInstance()->m_min_y, m_min_y);
    cgefParam::GetInstance()->m_max_y = std::max(cgefParam::GetInstance()->m_max_y, m_max_y);

    return 0;
}

int readCellgemTask::getInfo() {
    int i = 0, k = 0;
    char *ptr = m_pbuf;

    int len = 0, x = 0, y = 0, umi = 0;
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 1:
                    x = atoi(ptr);
                    m_min_x = std::min(m_min_x, x);
                    m_max_x = std::max(m_max_x, x);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 2:
                    y = atoi(ptr);
                    m_min_y = std::min(m_min_y, y);
                    m_max_y = std::max(m_max_y, y);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 3:
                    k = 0;
                    ptr = &m_pbuf[i + 1];
                    break;
                default:
                    break;
            }
        }
    }
    return 0;
}

////////////////////////////////////////////////raw////////////////////////////////////////////////////
int readCellgemTask_raw::getInfo() {
    int i = 0, k = 0;
    char *ptr = m_pbuf;

    string gid;
    int len = 0, x = 0, y = 0, umi = 0;
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    len = &m_pbuf[i] - ptr;
                    gid.clear();
                    gid.append(ptr, len);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 1:
                    x = atoi(ptr);
                    m_min_x = std::min(m_min_x, x);
                    m_max_x = std::max(m_max_x, x);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 2:
                    y = atoi(ptr);
                    m_min_y = std::min(m_min_y, y);
                    m_max_y = std::max(m_max_y, y);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 3:
                    k = 0;
                    ptr = &m_pbuf[i + 1];
                    if (m_map_bgene.find(gid) == m_map_bgene.end()) {
                        bgef_gene *gptr = new bgef_gene();
                        m_map_bgene.emplace(gid, gptr);
                    }
                    m_map_bgene[gid]->add(x, y, umi);
                    break;
                default:
                    break;
            }
        }
    }
    return m_map_bgene.size();
}

int readCellgemTask_raw::mergeinfo() {
    lock_guard<mutex> lock(m_mergemtx);

    cgefParam::GetInstance()->m_min_x = std::min(cgefParam::GetInstance()->m_min_x, m_min_x);
    cgefParam::GetInstance()->m_min_y = std::min(cgefParam::GetInstance()->m_min_y, m_min_y);
    cgefParam::GetInstance()->m_max_x = std::max(cgefParam::GetInstance()->m_max_x, m_max_x);
    cgefParam::GetInstance()->m_max_y = std::max(cgefParam::GetInstance()->m_max_y, m_max_y);

    auto itor_g = m_map_bgene.begin();
    auto &tmp_gene = cgefParam::GetInstance()->m_map_bgene;
    for (; itor_g != m_map_bgene.end(); itor_g++) {
        if (tmp_gene.find(itor_g->first) != tmp_gene.end()) {
            tmp_gene[itor_g->first]->merge(*(itor_g->second));
            delete itor_g->second;
        } else {
            tmp_gene.emplace(itor_g->first, itor_g->second);
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int readCellgemTask_cell::getInfo() {
    std::function<int(readCellgemTask_cell *)> parse;

    if (cgefParam::GetInstance()->has_gene_name) {
        if (m_bexon) {
            parse = &readCellgemTask_cell::getdataWithGenename_exon;
        } else {
            parse = &readCellgemTask_cell::getdataWithGenename;
        }
    } else {
        if (m_bexon) {
            parse = &readCellgemTask_cell::getdata_exon;
        } else {
            parse = &readCellgemTask_cell::getdata;
        }
    }
    return parse(this);
}

int readCellgemTask_cell::getdata() {
    int i = 0, k = 0, celllabel = 0;
    char *ptr = m_pbuf;

    string gid;
    int len = 0, x = 0, y = 0, umi = 0, exon = 0;
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    len = &m_pbuf[i] - ptr;
                    gid.clear();
                    gid.append(ptr, len);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 1:
                    x = atoi(ptr);
                    m_min_x = std::min(m_min_x, x);
                    m_max_x = std::max(m_max_x, x);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 2:
                    y = atoi(ptr);
                    m_min_y = std::min(m_min_y, y);
                    m_max_y = std::max(m_max_y, y);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 3:
                    umi = atoi(ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 4:
                    k = 0;
                    celllabel = atoi(ptr);
                    ptr = &m_pbuf[i + 1];

                    if (celllabel > 0) {
                        if (m_map_cell.find(celllabel) == m_map_cell.end()) {
                            cgef_cell *cptr = new cgef_cell(celllabel);
                            m_map_cell.emplace(celllabel, cptr);
                        }
                        m_map_cell[celllabel]->add(gid, umi, x, y);

                        // if(m_map_gene.find(gid) == m_map_gene.end())
                        // {
                        //     cgef_gene *gptr = new cgef_gene();
                        //     m_map_gene.emplace(gid, gptr);
                        // }
                        // m_map_gene[gid]->add(celllabel, umi);
                        if (m_map_gene_id.find(gid) == m_map_gene_id.end()) {
                            m_map_gene_id.emplace(gid, 0);
                        }
                    }

                    break;
                default:
                    break;
            }
        }
    }

    return m_map_cell.size();
}

int readCellgemTask_cell::getdata_exon() {
    int i = 0, k = 0, celllabel = 0;
    char *ptr = m_pbuf;

    string gid;
    int len = 0, x = 0, y = 0, umi = 0, exon = 0;
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    len = &m_pbuf[i] - ptr;
                    gid.clear();
                    gid.append(ptr, len);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 1:
                    x = atoi(ptr);
                    m_min_x = std::min(m_min_x, x);
                    m_max_x = std::max(m_max_x, x);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 2:
                    y = atoi(ptr);
                    m_min_y = std::min(m_min_y, y);
                    m_max_y = std::max(m_max_y, y);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 3:
                    umi = atoi(ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 4:
                    exon = atoi(ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 5:
                    k = 0;
                    celllabel = atoi(ptr);
                    ptr = &m_pbuf[i + 1];

                    if (celllabel > 0) {
                        if (m_map_cell.find(celllabel) == m_map_cell.end()) {
                            cgef_cell *cptr = new cgef_cell(celllabel);
                            m_map_cell.emplace(celllabel, cptr);
                        }
                        m_map_cell[celllabel]->add(gid, umi, x, y, exon);

                        // if(m_map_gene.find(gid) == m_map_gene.end())
                        // {
                        //     cgef_gene *gptr = new cgef_gene();
                        //     m_map_gene.emplace(gid, gptr);
                        // }
                        // m_map_gene[gid]->add(celllabel, umi, exon);
                        if (m_map_gene_id.find(gid) == m_map_gene_id.end()) {
                            m_map_gene_id.emplace(gid, 0);
                        }
                    }

                    break;
                default:
                    break;
            }
        }
    }

    return m_map_cell.size();
}

int readCellgemTask_cell::getdataWithGenename() {
    int i = 0, k = 0, celllabel = 0;
    char *ptr = m_pbuf;

    string gid;
    string gname;
    int len = 0, x = 0, y = 0, umi = 0, exon = 0;
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    len = &m_pbuf[i] - ptr;
                    gid.clear();
                    gid.append(ptr, len);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 1:
                    len = &m_pbuf[i] - ptr;
                    gname.clear();
                    gname.append(ptr, len);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 2:
                    x = atoi(ptr);
                    m_min_x = std::min(m_min_x, x);
                    m_max_x = std::max(m_max_x, x);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 3:
                    y = atoi(ptr);
                    m_min_y = std::min(m_min_y, y);
                    m_max_y = std::max(m_max_y, y);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 4:
                    umi = atoi(ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 5:
                    k = 0;
                    celllabel = atoi(ptr);
                    ptr = &m_pbuf[i + 1];

                    if (celllabel > 0) {
                        if (m_map_cell.find(celllabel) == m_map_cell.end()) {
                            cgef_cell *cptr = new cgef_cell(celllabel);
                            m_map_cell.emplace(celllabel, cptr);
                        }
                        m_map_cell[celllabel]->add(gid, umi, x, y);

                        if (m_map_gene_id.find(gid) == m_map_gene_id.end()) {
                            m_map_gene_id.emplace(gid, 0);
                        }
                        if (m_genename_map_per_t.find(gid) == m_genename_map_per_t.end()) {
                            m_genename_map_per_t.emplace(gid, gname);
                        }
                    }

                    break;
                default:
                    break;
            }
        }
    }

    return m_map_cell.size();
}

int readCellgemTask_cell::getdataWithGenename_exon() {
    int i = 0, k = 0, celllabel = 0;
    char *ptr = m_pbuf;

    string gid;
    string gname;
    int len = 0, x = 0, y = 0, umi = 0, exon = 0;
    for (; i < m_buflen; i++) {
        if (m_pbuf[i] == ',' || m_pbuf[i] == ';' || m_pbuf[i] == '\t' || m_pbuf[i] == '\n') {
            switch (k) {
                case 0:
                    len = &m_pbuf[i] - ptr;
                    gid.clear();
                    gid.append(ptr, len);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 1:
                    len = &m_pbuf[i] - ptr;
                    gname.clear();
                    gname.append(ptr, len);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 2:
                    x = atoi(ptr);
                    m_min_x = std::min(m_min_x, x);
                    m_max_x = std::max(m_max_x, x);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 3:
                    y = atoi(ptr);
                    m_min_y = std::min(m_min_y, y);
                    m_max_y = std::max(m_max_y, y);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 4:
                    umi = atoi(ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 5:
                    exon = atoi(ptr);
                    k++;
                    ptr = &m_pbuf[i + 1];
                    break;
                case 6:
                    k = 0;
                    celllabel = atoi(ptr);
                    ptr = &m_pbuf[i + 1];

                    if (celllabel > 0) {
                        if (m_map_cell.find(celllabel) == m_map_cell.end()) {
                            cgef_cell *cptr = new cgef_cell(celllabel);
                            m_map_cell.emplace(celllabel, cptr);
                        }
                        m_map_cell[celllabel]->add(gid, umi, x, y, exon);

                        if (m_map_gene_id.find(gid) == m_map_gene_id.end()) {
                            m_map_gene_id.emplace(gid, 0);
                        }
                        if (m_genename_map_per_t.find(gid) == m_genename_map_per_t.end()) {
                            m_genename_map_per_t.emplace(gid, gname);
                        }
                    }

                    break;
                default:
                    break;
            }
        }
    }

    return m_map_cell.size();
}