/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:23
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-21 11:35:58
 * @Description: file content
 */
#include "cellAdjust.h"

#include <fstream>
#include <sstream>

#include "bgef_writer.h"
#include "bin_task.h"
#include "dnb_merge_task.h"
#include "getsapdataTask.h"
#include "timer.h"

std::mutex getsapdataTask::m_mtx;
std::mutex getLabelInfoTask::m_mtx_;
cellAdjust::cellAdjust() {}

cellAdjust::~cellAdjust() {
    if (m_cell_arrayptr) {
        free(m_cell_arrayptr);
    }
    if (m_borderdataPtr) {
        free(m_borderdataPtr);
    }
    if (m_bgeffile_id) {
        H5Fclose(m_bgeffile_id);
    }
    if (generate_bgef_thread_.joinable()) generate_bgef_thread_.join();
}

void cellAdjust::readBgef(const string &strinput) {
    timer st(__FUNCTION__);
    m_bgeffile_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t dims[1];
    hid_t gene_did = H5Dopen(m_bgeffile_id, "/geneExp/bin1/gene", H5P_DEFAULT);
    hid_t gene_sid = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(gene_sid, dims, nullptr);

    m_genencnt = dims[0];
    Gene *genePtr = (Gene *)malloc(dims[0] * sizeof(Gene));
    hid_t genememtype, strtype;
    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 64);

    genememtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(genememtype, "gene", HOFFSET(Gene, gene), strtype);
    H5Tinsert(genememtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(genememtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);
    H5Dread(gene_did, genememtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, genePtr);
    H5Tclose(genememtype);
    H5Sclose(gene_sid);
    H5Dclose(gene_did);

    hid_t exp_did = H5Dopen(m_bgeffile_id, "/geneExp/bin1/expression", H5P_DEFAULT);
    hid_t exp_sid = H5Dget_space(exp_did);
    H5Sget_simple_extent_dims(exp_sid, dims, nullptr);

    m_geneexpcnt = dims[0];

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    Expression *expPtr = (Expression *)calloc(dims[0], sizeof(Expression));
    H5Dread(exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, expPtr);

    if (H5Lexists(m_bgeffile_id, "/geneExp/bin1/exon", H5P_DEFAULT) > 0) {
        m_bexon = true;
        hsize_t edims[1];
        hid_t did = H5Dopen(m_bgeffile_id, "/geneExp/bin1/exon", H5P_DEFAULT);
        hid_t sid = H5Dget_space(did);
        H5Sget_simple_extent_dims(sid, edims, nullptr);
        assert(edims[0] == m_geneexpcnt);
        unsigned int *exonPtr = new unsigned int[edims[0]];
        H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, exonPtr);
        H5Sclose(sid);
        H5Dclose(did);
        for (uint64_t i = 0; i < m_geneexpcnt; i++) {
            expPtr[i].exon = exonPtr[i];
        }
        delete[] exonPtr;
    }

    hid_t attr = H5Aopen(exp_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_x);
    attr = H5Aopen(exp_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_y);
    attr = H5Aopen(exp_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_x);
    attr = H5Aopen(exp_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_y);
    attr = H5Aopen(exp_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_resolution);
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", m_min_x, m_min_y, m_max_x, m_max_y);
    H5Aclose(attr);
    H5Tclose(memtype);
    H5Sclose(exp_sid);
    H5Dclose(exp_did);
    H5Tclose(strtype);

    if (H5Aexists(m_bgeffile_id, "omics") > 0) {
        hid_t omics_strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(omics_strtype, 32);
        hid_t fattr = H5Aopen(m_bgeffile_id, "omics", H5P_DEFAULT);
        H5Aread(fattr, omics_strtype, m_szomics);
        H5Tclose(omics_strtype);
    } else {
        strcpy(m_szomics, "Transcriptomics");
    }

    uint64_t l_id = 0;
    for (int i = 0; i < m_genencnt; i++) {
        m_vecgenename.emplace_back(genePtr[i].gene);
        Expression *ptr = expPtr + genePtr[i].offset;
        for (int j = 0; j < genePtr[i].count; j++) {
            l_id = ptr[j].x;
            l_id = (l_id << 32) | ptr[j].y;

            if (m_hash_vecdnb_exon.find(l_id) == m_hash_vecdnb_exon.end()) {
                vector<Dnbs_exon> tvec;
                m_hash_vecdnb_exon.emplace(l_id, tvec);
            }
            m_hash_vecdnb_exon[l_id].emplace_back(i, ptr[j].count, ptr[j].exon);
        }
    }
    log_info << util::Format("gene:{0} geneexp:{1} hashcnt:{2}", m_genencnt, m_geneexpcnt, m_hash_vecdnb_exon.size());

    free(genePtr);
    free(expPtr);
}

void cellAdjust::readCgef(const string &strinput) {
    timer st(__FUNCTION__);
    hid_t file_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t cell_did = H5Dopen(file_id, "/cellBin/cell", H5P_DEFAULT);

    hsize_t dims[1];
    hid_t cell_sid = H5Dget_space(cell_did);
    H5Sget_simple_extent_dims(cell_sid, dims, nullptr);

    m_cellcnt = dims[0];
    hid_t memtype = getMemtypeOfCellData();
    m_cell_arrayptr = (CellData *)malloc(dims[0] * sizeof(CellData));
    H5Dread(cell_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cell_arrayptr);

    H5Tclose(memtype);
    H5Sclose(cell_sid);
    H5Dclose(cell_did);

    // hid_t d_id = H5Dopen(file_id, "/cellBin/blockIndex", H5P_DEFAULT);
    // hid_t s_id = H5Dget_space(d_id);
    // H5Sget_simple_extent_dims(s_id, dims, nullptr);
    // m_blkidxPtr = static_cast<unsigned int *>(calloc(dims[0], sizeof(unsigned
    // int))); H5Dread(d_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
    // m_blkidxPtr);

    // H5Sclose(s_id);
    // H5Dclose(d_id);

    hid_t d_id = H5Dopen(file_id, "/cellBin/blockSize", H5P_DEFAULT);
    if (d_id < 0) {
        log_error << errorCode::E_MISSINGFILEINFO << "can't find /cellBin/blockSize in file";
    } else {
        H5Dread(d_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_block_size);
        H5Dclose(d_id);
    }

    hsize_t bdims[3];
    hid_t border_id = H5Dopen(file_id, "/cellBin/cellBorder", H5P_DEFAULT);
    hid_t border_sid = H5Dget_space(border_id);
    H5Sget_simple_extent_dims(border_sid, bdims, nullptr);

    m_borderdataPtr = (short *)calloc(bdims[0] * bdims[1] * bdims[2], 2);
    H5Dread(border_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_borderdataPtr);

    hid_t cell_exp_did = H5Dopen(file_id, "/cellBin/cellExp", H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(cell_exp_did);
    H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
    int cellexpcnt = dims[0];
    if (isOlderCellExpDataVersion(file_id)) {
        isOldCellExpVersion = true;
        hid_t memtype = getMemtypeOfOlderCellExpData();
        m_olderCellExpPtr = (olderCellExpData *)malloc(dims[0] * sizeof(olderCellExpData));
        H5Dread(cell_exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_olderCellExpPtr);
        H5Tclose(memtype);
    } else {
        isOldCellExpVersion = false;
        hid_t memtype = getMemtypeOfCellExpData();
        m_cellexpPtr = (CellExpData *)malloc(dims[0] * sizeof(CellExpData));
        H5Dread(cell_exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cellexpPtr);
        H5Tclose(memtype);
    }
    H5Sclose(dataspace_id);
    H5Dclose(cell_exp_did);

    uint32_t fcnt = 0;
    uint64_t l_id = 0;
    vector<cv::Point> vecpoint;
    m_vecCellgem.reserve(m_geneexpcnt);

    int x, y;
    vector<cv::Point> vecborder;
    vector<cv::Point> rvecborder;
    // m_fill_points = Mat::zeros(m_max_y-m_min_y+1, m_max_x-m_min_x+1, CV_8UC1);
    short *ptmp = m_borderdataPtr;
    m_hash_filter_cells.clear();
    for (uint64_t i = 0; i < bdims[0]; i++) {
        vecborder.clear();
        for (uint64_t j = 0; j < bdims[1]; j++) {
            x = ptmp[2 * j];
            y = ptmp[2 * j + 1];
            if (x == SHRT_MAX && y == SHRT_MAX) {
                break;
            }
            x += m_cell_arrayptr[i].x;
            y += m_cell_arrayptr[i].y;
            vecborder.emplace_back(x, y);
        }

        if (!vecborder.empty()) {
            rvecborder.clear();
            cv::Rect rtmp = cv::boundingRect(vecborder);
            cv::Mat t = cv::Mat::zeros(rtmp.height, rtmp.width, CV_8UC1);
            for (cv::Point &pt : vecborder) {
                rvecborder.emplace_back(pt.x - rtmp.x, pt.y - rtmp.y);
            }
            cv::fillPoly(t, rvecborder, 1);
            cv::findNonZero(t, vecpoint);
            bool bfind = false;

            for (cv::Point &pt : vecpoint) {
                x = pt.x + rtmp.x;
                y = pt.y + rtmp.y;
                l_id = x;
                l_id = (l_id << 32) | y;
                auto dnb_itor = m_hash_vecdnb_exon.find(l_id);
                if (dnb_itor != m_hash_vecdnb_exon.end()) {
                    for (Dnbs_exon &dnbs : dnb_itor->second) {
                        m_vecCellgem.emplace_back(dnbs.geneid, x, y, dnbs.midcnt, i + 1);
                    }
                    m_hash_vecdnb_exon.erase(l_id);
                    bfind = true;
                }
            }

            if (bfind) {
                fcnt++;
            } else {
                log_info << i << " " << m_cell_arrayptr[i].dnb_count << " " << m_cell_arrayptr[i].area;
            }
        } else {
            if (isOldCellExpVersion) {
                map<uint32_t, uint16_t> gid_count;
                gid_count.emplace(m_olderCellExpPtr[i].gene_id, m_olderCellExpPtr[i].count);
                m_hash_filter_cells.emplace(i + 1, std::move(gid_count));
            } else {
                map<uint32_t, uint16_t> gid_count;
                gid_count.emplace(m_cellexpPtr[i].gene_id, m_cellexpPtr[i].count);
                m_hash_filter_cells.emplace(i + 1, std::move(gid_count));
            }
            log_info << "empty cid " << i;
        }
        ptmp += BORDERCNT * 2;
    }

    log_info << "cellcnt:" << m_cellcnt << " fcnt:" << fcnt;

    auto itor_s = m_hash_vecdnb_exon.begin();
    for (; itor_s != m_hash_vecdnb_exon.end(); itor_s++) {
        x = (itor_s->first) >> 32;
        y = (itor_s->first) & 0xFFFFFFFF;
        for (Dnbs_exon &dnbs : itor_s->second) {
            m_vecCellgem.emplace_back(dnbs.geneid, x, y, dnbs.midcnt, 0);
        }
    }

    int min_x, min_y, max_x, max_y;
    hid_t attr = H5Aopen(border_id, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &min_x);
    attr = H5Aopen(border_id, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &min_y);
    attr = H5Aopen(border_id, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_x);
    attr = H5Aopen(border_id, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_y);
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", min_x, min_y, max_x, max_y);

    attr = H5Aopen(file_id, "offsetX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetX);

    attr = H5Aopen(file_id, "offsetY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetY);
    printf("offsetx:%d offsety:%d\n", m_offsetX, m_offsetY);
    H5Aclose(attr);
    H5Sclose(border_sid);
    H5Dclose(border_id);
    H5Fclose(file_id);
}

uint32_t cellAdjust::getCellLabelgem(vector<string> &genename, vector<cellgem_label> &vecCellgem) {
    timer st(__FUNCTION__);
    genename.reserve(m_vecgenename.size());
    genename.insert(genename.end(), m_vecgenename.begin(), m_vecgenename.end());
    vecCellgem.swap(m_vecCellgem);
    return vecCellgem.size();
}

bool cellAdjust::addborder(unsigned int cid, vector<cv::Point> &vecPoint, vector<cv::Point> &border,
                           vector<short> &vec_border) {
    convexHull(vecPoint, border, true);
    if (border.size() < 3) {
        // printf("borderwarn cid=%d pot=%d bor=%d dnb=%d\n", cid, vecPoint.size(),
        // border.size(), m_cell_arrayptr[cid].dnb_count);
        return false;
    }

    int bsz = border.size();
    uint64_t i = 0;

    if (bsz > BORDERCNT) {
        vector<cv::Point> tmpborder;
        double epsilon = 0.01 * arcLength(border, true);
        approxPolyDP(border, tmpborder, epsilon, true);

        bsz = tmpborder.size();
        for (; i < bsz; i++) {
            vec_border.emplace_back(tmpborder[i].x - m_cell_arrayptr[cid].x);
            vec_border.emplace_back(tmpborder[i].y - m_cell_arrayptr[cid].y);
        }
    } else {
        for (; i < bsz; i++) {
            vec_border.emplace_back(border[i].x - m_cell_arrayptr[cid].x);
            vec_border.emplace_back(border[i].y - m_cell_arrayptr[cid].y);
        }
    }

    for (; i < BORDERCNT; i++)  // 不足补SHRT_MAX
    {
        vec_border.emplace_back(SHRT_MAX);
        vec_border.emplace_back(SHRT_MAX);
    }

    return true;
}

bool cellAdjust::ParseBorderFile(const string &strInput) {
    vector<string> cellAdjustDatas = readLines(strInput);
    if (cellAdjustDatas.empty()) {
        return false;
    }
    borderDatas.clear();

    for (uint64_t i = 0; i < cellAdjustDatas.size(); i++) {
        vector<string> cellAdjustRowDatas = split(cellAdjustDatas[i], '\t');

        if (!cellAdjustRowDatas.empty()) {
            unsigned int cellId = std::stoi(cellAdjustRowDatas[0]) - 1;
            vector<cv::Point> border;
            for (uint64_t i = 1; i < cellAdjustRowDatas.size(); i++) {
                vector<string> cellAdjustCoordinateDatas = split(cellAdjustRowDatas[i], ' ');
                border.emplace_back(std::stoi(cellAdjustCoordinateDatas[0]), std::stoi(cellAdjustCoordinateDatas[1]));
            }
            borderDatas.emplace(cellId, border);
        } else {
            return false;
        }
    }
    return true;
}

bool cellAdjust::AddBorderFromFile(unsigned int cid, vector<cv::Point> &border, vector<short> &vecBorder) {
    if (borderDatas.find(cid) != borderDatas.end()) {
        border = borderDatas[cid];
        uint64_t i = 0;
        for (; i < border.size(); i++) {
            vecBorder.emplace_back(border[i].x - m_cell_arrayptr[cid].x);
            vecBorder.emplace_back(border[i].y - m_cell_arrayptr[cid].y);
        }
        for (; i < BORDERCNT; i++) {
            vecBorder.emplace_back(SHRT_MAX);
            vecBorder.emplace_back(SHRT_MAX);
        }
        return true;
    } else {
        return false;
    }
}

void cellAdjust::writeCell(Cell *cellptr, unsigned int cellcnt, DnbExpression *dnbptr, unsigned int dnbcnt) {
    timer st(__FUNCTION__);
    uint16_t gene_count, exp_count, dnb_count, area, cell_type_id, maxExpmid = 0;
    map<uint32_t, uint16_t> map_gene_cnt;  // gid midcnt
    uint32_t offset = 0;
    vector<cv::Point> vecPoint;
    vector<cv::Point> border;
    log_info << util::Format("rawcellcnt:{0} newcellcnt:{1} dnbcnt:{2}", m_cellcnt, cellcnt, dnbcnt);

    unsigned int blocknum = m_block_size[2] * m_block_size[3];
    vector<vector<Cell>> vec_vec_cell;
    for (int i = 0; i < blocknum; i++) {
        vector<Cell> vectmp;
        vec_vec_cell.emplace_back(std::move(vectmp));
    }
    unsigned int blkid = 0, cid = 0;

    for (unsigned int ci = 0; ci < cellcnt; ci++) {
        cid = cellptr[ci].cellid - 1;
        CellData &cd = m_cell_arrayptr[cid];
        blkid = cd.x / m_block_size[0] + (cd.y / m_block_size[1]) * m_block_size[2];
        vec_vec_cell[blkid].emplace_back(cellptr[ci]);
    }
    auto itor_f = m_hash_filter_cells.begin();
    for (; itor_f != m_hash_filter_cells.end(); itor_f++) {
        CellData &cd = m_cell_arrayptr[itor_f->first];
        blkid = cd.x / m_block_size[0] + (cd.y / m_block_size[1]) * m_block_size[2];
        Cell old_cell_data = {itor_f->first + 1, 0, 0};
        vec_vec_cell[blkid].emplace_back(std::move(old_cell_data));
    }

    vector<uint32_t> vec_blkidx;
    vec_blkidx.reserve(blocknum + 1);
    vector<short> vec_border;
    vec_border.reserve(m_cellcnt * 2 * BORDERCNT);
    // vec_border.reserve(cellcnt*2*BORDERCNT);

    uint32_t newcid = 0, offcnt = 0, blkcnt = 0;
    int minx = INT_MAX, miny = INT_MAX, maxx = 0, maxy = 0;
    for (vector<Cell> &vecC : vec_vec_cell) {
        blkcnt = 0;
        for (Cell &ce : vecC) {
            cid = ce.cellid - 1;
            vecPoint.clear();
            border.clear();
            gene_count = 0;
            exp_count = 0;
            map_gene_cnt.clear();
            dnb_count = 0;

            auto cell_itor = m_hash_filter_cells.find(cid);
            if (cell_itor != m_hash_filter_cells.end()) {
                gene_count = m_cell_arrayptr[cell_itor->first].gene_count;
                exp_count = m_cell_arrayptr[cell_itor->first].exp_count;
                area = m_cell_arrayptr[cell_itor->first].area;
                dnb_count = m_cell_arrayptr[cell_itor->first].dnb_count;

                auto itor = cell_itor->second.begin();
                for (; itor != cell_itor->second.end(); itor++) {
                    m_cgefwPtr->cell_exp_list_.emplace_back(itor->first, itor->second);
                    maxExpmid = std::max(maxExpmid, itor->second);
                    if (m_map_gene.find(itor->first) == m_map_gene.end()) {
                        vector<GeneExpData> tvec;
                        m_map_gene.emplace(itor->first, tvec);
                    }
                    m_map_gene[itor->first].emplace_back(newcid, itor->second);
                }

                short *ptr = m_borderdataPtr + cid * 2 * BORDERCNT;
                vec_border.insert(vec_border.end(), ptr, ptr + 2 * BORDERCNT);
            } else {
                for (uint32_t i = 0; i < ce.count; i++) {
                    DnbExpression &dnb = dnbptr[ce.offset + i];
                    if (map_gene_cnt.find(dnb.gene_id) != map_gene_cnt.end()) {
                        map_gene_cnt[dnb.gene_id] += dnb.count;
                    } else {
                        gene_count++;
                        map_gene_cnt.emplace(dnb.gene_id, dnb.count);
                    }
                    exp_count += dnb.count;
                    vecPoint.emplace_back(dnb.x, dnb.y);
                }
                dnb_count = ce.count;
                bool ret = false;
                if (!extend_method_) {
                    ret = addborder(cid, vecPoint, border, vec_border);
                } else {
                    ret = AddBorderFromFile(cid, border, vec_border);
                }
                if (!ret) {
                    short *ptr = m_borderdataPtr + cid * 2 * BORDERCNT;
                    vec_border.insert(vec_border.end(), ptr, ptr + 2 * BORDERCNT);
                }

                cv::Moments mu = cv::moments(border, true);
                cv::Rect trect = boundingRect(border);

                minx = std::min(minx, trect.x);
                maxx = std::max(maxx, trect.x + trect.width);
                miny = std::min(miny, trect.y);
                maxy = std::max(maxy, trect.y + trect.height);

                area = mu.m00;
                auto itor = map_gene_cnt.begin();
                for (; itor != map_gene_cnt.end(); itor++) {
                    m_cgefwPtr->cell_exp_list_.emplace_back(itor->first, itor->second);
                    maxExpmid = std::max(maxExpmid, itor->second);
                    if (m_map_gene.find(itor->first) == m_map_gene.end()) {
                        vector<GeneExpData> tvec;
                        m_map_gene.emplace(itor->first, tvec);
                    }
                    m_map_gene[itor->first].emplace_back(newcid, itor->second);
                }
            }

            cell_type_id =
                m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand() % (m_cgefwPtr->random_cell_type_num_ + 1);
            CellData cell = {newcid++,
                             m_cell_arrayptr[cid].x,
                             m_cell_arrayptr[cid].y,
                             offset,
                             gene_count,
                             exp_count,
                             dnb_count,
                             area,
                             cell_type_id};
            offset += gene_count;

            m_cgefwPtr->cell_attr_.min_x = std::min(m_cgefwPtr->cell_attr_.min_x, cell.x);
            m_cgefwPtr->cell_attr_.max_x = std::max(m_cgefwPtr->cell_attr_.max_x, cell.x);

            m_cgefwPtr->cell_attr_.min_y = std::min(m_cgefwPtr->cell_attr_.min_y, cell.y);
            m_cgefwPtr->cell_attr_.max_y = std::max(m_cgefwPtr->cell_attr_.max_y, cell.y);

            m_cgefwPtr->cell_attr_.min_area = std::min(m_cgefwPtr->cell_attr_.min_area, area);
            m_cgefwPtr->cell_attr_.max_area = std::max(m_cgefwPtr->cell_attr_.max_area, area);

            m_cgefwPtr->cell_attr_.min_gene_count = std::min(m_cgefwPtr->cell_attr_.min_gene_count, gene_count);
            m_cgefwPtr->cell_attr_.max_gene_count = std::max(m_cgefwPtr->cell_attr_.max_gene_count, gene_count);

            m_cgefwPtr->cell_attr_.min_exp_count = std::min(m_cgefwPtr->cell_attr_.min_exp_count, exp_count);
            m_cgefwPtr->cell_attr_.max_exp_count = std::max(m_cgefwPtr->cell_attr_.max_exp_count, exp_count);

            m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, cell.dnb_count);
            m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, cell.dnb_count);

            m_cgefwPtr->expression_num_ += gene_count;
            m_cgefwPtr->exp_count_sum_ += exp_count;
            m_cgefwPtr->dnb_count_sum_ += cell.dnb_count;
            m_cgefwPtr->area_sum_ += area;
            m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
            blkcnt++;
        }

        vec_blkidx.emplace_back(offcnt);
        offcnt += blkcnt;
    }

    vec_blkidx.emplace_back(newcid);
    m_cgefwPtr->cell_num_ = newcid;
    m_cgefwPtr->max_mid_count_ = maxExpmid;

    int effective_rect[4] = {minx, miny, maxx, maxy};
    m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), m_cgefwPtr->cell_num_, effective_rect);
    m_cgefwPtr->storeCell(m_block_size[2] * m_block_size[3], vec_blkidx.data(), m_block_size);
    m_cgefwPtr->storeCellExp();
    m_cgefwPtr->storeCellTypeList();
}

void cellAdjust::writeGene() {
    timer st(__FUNCTION__);
    log_info << util::Format("genecnt:{0} hashcnt:{1} geneexpcnt:{2}", m_genencnt, m_map_gene.size(),
                             m_cgefwPtr->expression_num_);
    m_cgefwPtr->gene_num_ = m_genencnt;
    GeneData *gene_data_list = static_cast<GeneData *>(calloc(m_cgefwPtr->gene_num_, sizeof(GeneData)));

    unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
    unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
    unsigned short max_MID_count = 0;
    vector<GeneExpData> gene_exp_list;
    gene_exp_list.reserve(m_cgefwPtr->expression_num_);

    m_cgefwPtr->max_mid_count_ = 0;
    for (uint32_t i = 0; i < m_genencnt; i++) {
        exp_count = 0;
        max_MID_count = 0;
        auto itor = m_map_gene.find(i);
        string &strgene = m_vecgenename[i];
        if (itor != m_map_gene.end()) {
            for (GeneExpData &gexp : itor->second) {
                gene_exp_list.emplace_back(gexp);
                max_MID_count = std::max(max_MID_count, gexp.count);
                m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, gexp.count);
                exp_count += gexp.count;
            }

            cell_count = itor->second.size();
            gene_data_list[i].cell_count = cell_count;
            gene_data_list[i].exp_count = exp_count;
            memcpy(gene_data_list[i].gene_name, strgene.c_str(), strgene.length());
            gene_data_list[i].max_mid_count = max_MID_count;
            gene_data_list[i].offset = offset;
            offset += cell_count;
        } else {
            memcpy(gene_data_list[i].gene_name, strgene.c_str(), strgene.length());
            gene_data_list[i].cell_count = 0;
            gene_data_list[i].exp_count = 0;
            gene_data_list[i].max_mid_count = 0;
            gene_data_list[i].offset = 0;
        }

        m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, max_MID_count);
        min_exp_count = std::min(min_exp_count, exp_count);
        max_exp_count = std::max(max_exp_count, exp_count);
        min_cell_count = std::min(min_cell_count, cell_count);
        max_cell_count = std::max(max_cell_count, cell_count);
    }

    m_cgefwPtr->expression_num_ = gene_exp_list.size();
    m_cgefwPtr->storeGeneAndGeneExp(min_exp_count, max_exp_count, min_cell_count, max_cell_count, gene_data_list,
                                    gene_exp_list);
    free(gene_data_list);
}

void cellAdjust::writeCellAdjust(const string &outpath, const string &outline_path, Cell *cellptr, int cellcnt,
                                 DnbExpression *dnbptr, int dnbcnt) {
    if (outline_path.empty()) {
        printf("No cell outline file, will be handled by default");
    } else {
        if (!ParseBorderFile(outline_path)) {
            printf("Can not parse input cell border file");
            return;
        }
        extend_method_ = true;
    }
    m_cgefwPtr = new CgefWriter();
    m_cgefwPtr->setOutput(outpath);
    CellBinAttr cell_bin_attr = {/*.version = */ 2,
                                 /*.resolution = */ m_resolution,
                                 /*.offsetX = */ m_min_x,
                                 /*.offsetY = */ m_min_y};
    m_cgefwPtr->storeAttr(cell_bin_attr);
    writeCell(cellptr, cellcnt, dnbptr, dnbcnt);
    writeGene();
    delete m_cgefwPtr;
}

herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata) {
    auto group_names = reinterpret_cast<std::vector<std::string> *>(opdata);
    group_names->push_back(name);
    return 0;
}

void cellAdjust::createRegionGef(const string &out) {
    timer st(__FUNCTION__);
    hid_t gid = H5Gopen(m_bgeffile_id, "/geneExp", H5P_DEFAULT);
    if (gid < 0) {
        std::cout << "can not find input spatial bin gef file... " << std::endl;
        return;
    }
    std::vector<std::string> group_names;
    herr_t idx = H5Literate(gid, H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, &group_names);
    H5Gclose(gid);

    m_bgefopts->bin_sizes_.clear();
    for (string &str : group_names) {
        int bin = std::stoi(str.substr(3));
        m_bgefopts->bin_sizes_.push_back(bin);
    }

    if (m_bgefopts->map_gene_exp_.size() == 0) {
        log_error << "region has no information. please check. ";
        return;
    }

    m_bgefopts->m_genes_queue.init(m_bgefopts->map_gene_exp_.size());
    ThreadPool thpool(m_bgefopts->thread_ * 2);

    m_bgefopts->m_stromics.append(m_szomics);
    BgefWriter bgef_writer(out, false, m_bexon, m_bgefopts->m_stromics);
    bgef_writer.setResolution(m_resolution);

    int genecnt = 0;
    for (unsigned int bin : m_bgefopts->bin_sizes_) {
        auto &dnb_matrix = m_bgefopts->dnbmatrix_;
        auto &dnbAttr = m_bgefopts->dnbmatrix_.dnb_attr;

        dnbAttr.min_x = (m_min_x / bin) * bin;
        dnbAttr.len_x = (m_maxx) / bin + 1;
        dnbAttr.min_y = (m_min_y / bin) * bin;
        dnbAttr.len_y = (m_maxy) / bin + 1;

        dnbAttr.max_x = (m_maxx / bin) * bin;
        dnbAttr.max_y = (m_maxy / bin) * bin;

        dnbAttr.max_gene = 0;
        dnbAttr.max_mid = 0;
        dnbAttr.number = 0;
        unsigned long matrix_len = (unsigned long)(dnbAttr.len_x) * dnbAttr.len_y;

        log_info << "bin " << bin << " matrix: min_x=" << dnbAttr.min_x << " len_x=" << dnbAttr.len_x
                 << " min_y=" << dnbAttr.min_y << " len_y=" << dnbAttr.len_y << " matrix_len=" << matrix_len;

        if (bin == 1) {
            dnb_matrix.pmatrix_us = (BinStatUS *)calloc(matrix_len, sizeof(BinStatUS));
            assert(dnb_matrix.pmatrix_us);
            if (!dnb_matrix.pmatrix_us) {
                log_error << errorCode::E_ALLOCMEMORYFAILED << "can not alloc memory for wholeExp matrix. ";
                return;
            }
            if (m_bexon) {
                dnb_matrix.pexon16 = (unsigned short *)calloc(matrix_len, 2);
                assert(dnb_matrix.pexon16);
                if (!dnb_matrix.pexon16) {
                    log_error << errorCode::E_ALLOCMEMORYFAILED << "can not alloc memory for wholeExp matrix. ";
                    return;
                }
            }
        } else {
            dnb_matrix.pmatrix = (BinStat *)calloc(matrix_len, sizeof(BinStat));
            assert(dnb_matrix.pmatrix);
            if (!dnb_matrix.pmatrix) {
                log_error << errorCode::E_ALLOCMEMORYFAILED << "can not alloc memory for wholeExp matrix. ";
                return;
            }
            if (m_bexon) {
                dnb_matrix.pexon32 = (unsigned int *)calloc(matrix_len, 4);
                assert(dnb_matrix.pexon32);
                if (!dnb_matrix.pexon32) {
                    log_error << errorCode::E_ALLOCMEMORYFAILED << "can not alloc memory for wholeExp matrix. ";
                    return;
                }
            }
        }

        for (int i = 0; i < m_bgefopts->thread_; i++) {
            auto *task = new DnbMergeTask(m_bgefopts->map_gene_exp_.size(), i, bin);
            thpool.addTask(task);
        }

        auto itor = m_bgefopts->map_gene_exp_.begin();
        for (; itor != m_bgefopts->map_gene_exp_.end(); itor++) {
            auto *task = new BinTask(bin, itor->first.c_str());
            thpool.addTask(task);
        }

        unsigned int offset = 0;
        unsigned int maxexp = 0;
        unsigned int maxexon = 0;
        genecnt = 0;
        map<string, vector<Expression>> gene_info;
        while (true) {
            GeneInfo *pgeneinfo = m_bgefopts->m_geneinfo_queue.getPtr();
            gene_info.insert(map<string, vector<Expression>>::value_type(pgeneinfo->geneid, *pgeneinfo->vecptr));

            maxexp = std::max(maxexp, pgeneinfo->maxexp);
            maxexon = std::max(maxexon, pgeneinfo->maxexon);

            if (bin == 100) {
                m_bgefopts->m_vec_bin100.emplace_back(pgeneinfo->geneid, pgeneinfo->umicnt, pgeneinfo->e10);
            }

            delete pgeneinfo;
            genecnt++;
            if (genecnt == m_bgefopts->map_gene_exp_.size()) {
                break;
            }
        }

        for (auto itor : gene_info) {
            if (bin == 1) {
                m_bgefopts->expressions_.insert(m_bgefopts->expressions_.end(), itor.second.begin(), itor.second.end());
            } else {
                for (auto g : itor.second) {
                    g.x *= bin;
                    g.y *= bin;
                    m_bgefopts->expressions_.push_back(std::move(g));
                }
            }

            m_bgefopts->genes_.emplace_back(itor.first.c_str(), offset, static_cast<unsigned int>(itor.second.size()));
            offset += itor.second.size();
        }

        bgef_writer.storeGene(m_bgefopts->expressions_, m_bgefopts->genes_, dnb_matrix.dnb_attr, maxexp, bin);
        bgef_writer.storeGeneExon(m_bgefopts->expressions_, maxexon, bin);
        m_bgefopts->expressions_.clear();
        m_bgefopts->genes_.clear();

        thpool.waitTaskDone();
        m_bgefopts->m_genes_queue.clear(bin);

        // write dnb
        if (bin == 100) {
            vector<GeneStat> &geneStat = m_bgefopts->m_vec_bin100;
            std::sort(geneStat.begin(), geneStat.end(), [](const GeneStat &p1, const GeneStat &p2) {
                if (p1.mid_count > p2.mid_count) {
                    return true;
                } else if (p1.mid_count == p2.mid_count) {
                    int ret = strcmp(p1.gene, p2.gene);
                    return ret < 0;
                } else {
                    return false;
                }
            });
            bgef_writer.storeStat(geneStat);
        }

        vector<unsigned int> vec_mid;
        unsigned long number = 0;

        if (bin == 1) {
            for (unsigned long i = 0; i < matrix_len; i++) {
                if (dnb_matrix.pmatrix_us[i].gene_count) {
                    ++number;
                    vec_mid.push_back(dnb_matrix.pmatrix_us[i].mid_count);
                }
            }
        } else {
            for (unsigned long i = 0; i < matrix_len; i++) {
                if (dnb_matrix.pmatrix[i].gene_count) {
                    ++number;
                    vec_mid.push_back(dnb_matrix.pmatrix[i].mid_count);
                }
            }
        }

        int sz = vec_mid.size();
        sort(vec_mid.begin(), vec_mid.end(), [](const unsigned int a, const unsigned int b) { return a < b; });
        if (bin > 50) {
            dnbAttr.max_mid = vec_mid[sz - 1];
        } else {
            int limit = sz * 0.999;
            dnbAttr.max_mid = vec_mid[limit];
        }

        dnbAttr.number = number;
        bgef_writer.storeDnb(dnb_matrix, bin);
        bgef_writer.storeWholeExon(dnb_matrix, bin);

        if (bin == 1) {
            if (dnb_matrix.pmatrix_us != nullptr) {
                free(dnb_matrix.pmatrix_us);
                dnb_matrix.pmatrix_us = nullptr;
                if (m_bexon) {
                    free(dnb_matrix.pexon16);
                    dnb_matrix.pexon16 = nullptr;
                }
            }
        } else {
            if (dnb_matrix.pmatrix != nullptr) {
                free(dnb_matrix.pmatrix);
                dnb_matrix.pmatrix = nullptr;
                if (m_bexon) {
                    free(dnb_matrix.pexon32);
                    dnb_matrix.pexon32 = nullptr;
                }
            }
        }
    }
}

void cellAdjust::getRegionGenedata(vector<vector<int>> &m_vecpos) {
    timer st(__FUNCTION__);
    m_bgefopts = BgefOptions::GetInstance();
    m_bgefopts->map_gene_exp_.clear();
    int num = m_vecpos.size();
    uint64_t l_id = 0;
    vector<cv::Point> non_zerovecpoint;
    vector<cv::Point> relativepoint;
    int x, y;
    unsigned long totalSize = 0;

    for (uint32_t i = 0; i < num; i++) {
        relativepoint.clear();
        non_zerovecpoint.clear();

        int cnt = m_vecpos[i].size();
        int *pos = m_vecpos[i].data();
        int minx = INT_MAX, miny = INT_MAX, maxx = 0, maxy = 0;
        for (uint32_t j = 0; j < cnt; j += 2) {
            minx = std::min(minx, pos[j]);
            maxx = std::max(maxx, pos[j]);
            miny = std::min(miny, pos[j + 1]);
            maxy = std::max(maxy, pos[j + 1]);
        }
        m_maxx = std::max(m_maxx, maxx);
        m_maxy = std::max(m_maxy, maxy);

        for (uint32_t j = 0; j < cnt; j += 2) {
            relativepoint.emplace_back(pos[j] - minx, pos[j + 1] - miny);
        }

        int rows = maxy - miny + 1;
        int cols = maxx - minx + 1;
        cv::Mat fill_points = cv::Mat::zeros(rows, cols, CV_8UC1);
        fillPoly(fill_points, relativepoint, 1);
        findNonZero(fill_points, non_zerovecpoint);

        for (cv::Point &pt : non_zerovecpoint) {
            x = pt.x + minx;
            y = pt.y + miny;
            l_id = x;
            l_id = (l_id << 32) | y;
            auto dnb_itor = m_hash_vecdnb_exon.find(l_id);
            if (dnb_itor != m_hash_vecdnb_exon.end()) {
                for (Dnbs_exon &dnbs : dnb_itor->second) {
                    string str(m_vecgenename[dnbs.geneid]);
                    auto itor_t = m_bgefopts->map_gene_exp_.find(str);
                    if (itor_t == m_bgefopts->map_gene_exp_.end()) {
                        vector<Expression> tvec;
                        m_bgefopts->map_gene_exp_.emplace(str, std::move(tvec));
                    }
                    m_bgefopts->map_gene_exp_[str].emplace_back(x, y, dnbs.midcnt, dnbs.exon);
                }
                m_hash_vecdnb_exon.erase(l_id);
                totalSize += dnb_itor->second.size();
            }
        }
    }

    m_bgefopts->expressions_.clear();
    m_bgefopts->genes_.clear();
    m_bgefopts->expressions_.reserve(totalSize);
    m_bgefopts->genes_.reserve(m_bgefopts->map_gene_exp_.size());
}

//////////////////////////////cgef lasso////////////////////////////////////
void cellAdjust::getRegionCelldata(vector<vector<int>> &m_vecpos) {
    timer st(__FUNCTION__);
    int num = m_vecpos.size();
    uint64_t l_id = 0;
    vector<cv::Point> non_zerovecpoint;
    vector<cv::Point> relativepoint;
    int x, y;

    for (uint32_t i = 0; i < num; i++) {
        relativepoint.clear();
        non_zerovecpoint.clear();

        int cnt = m_vecpos[i].size();
        int *pos = m_vecpos[i].data();
        int minx = INT_MAX, miny = INT_MAX, maxx = 0, maxy = 0;
        for (uint32_t j = 0; j < cnt; j += 2) {
            minx = std::min(minx, pos[j]);
            maxx = std::max(maxx, pos[j]);
            miny = std::min(miny, pos[j + 1]);
            maxy = std::max(maxy, pos[j + 1]);
        }

        for (uint32_t j = 0; j < cnt; j += 2) {
            relativepoint.emplace_back(pos[j] - minx, pos[j + 1] - miny);
        }

        int rows = maxy - miny + 1;
        int cols = maxx - minx + 1;
        cv::Mat fill_points = cv::Mat::zeros(rows, cols, CV_8UC1);
        cv::fillPoly(fill_points, relativepoint, 1);
        cv::findNonZero(fill_points, non_zerovecpoint);
        for (cv::Point &pt : non_zerovecpoint) {
            x = pt.x + minx;
            y = pt.y + miny;
            l_id = x;
            l_id = (l_id << 32) | y;
            m_setcell.insert(l_id);
        }
    }
}

void cellAdjust::readRawCgef(const string &strcgef) {
    timer st(__FUNCTION__);
    hid_t file_id = H5Fopen(strcgef.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        log_info << "can not open gef file. please check file. ";
        return;
    }
    if (H5Aexists(file_id, "omics") > 0) {
        hid_t omics_strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(omics_strtype, 32);
        hid_t fattr = H5Aopen(file_id, "omics", H5P_DEFAULT);
        H5Aread(fattr, omics_strtype, m_szomics);
        H5Aclose(fattr);
        H5Tclose(omics_strtype);
    } else {
        strcpy(m_szomics, "Transcriptomics");
    }

    uint32_t cellexpcnt = 0;
    // cell
    hid_t cell_did = H5Dopen(file_id, "/cellBin/cell", H5P_DEFAULT);
    hsize_t dims[1];
    hid_t cell_sid = H5Dget_space(cell_did);
    H5Sget_simple_extent_dims(cell_sid, dims, nullptr);
    m_cellcnt = dims[0];
    hid_t memtype = getMemtypeOfCellData();
    m_cell_arrayptr = (CellData *)malloc(dims[0] * sizeof(CellData));
    H5Dread(cell_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cell_arrayptr);
    H5Tclose(memtype);
    H5Sclose(cell_sid);
    H5Dclose(cell_did);

    // cell bor
    hid_t cell_bor_did = H5Dopen(file_id, "/cellBin/cellBorder", H5P_DEFAULT);
    hsize_t dims_b[3];
    hid_t dataspace_id = H5Dget_space(cell_bor_did);
    H5Sget_simple_extent_dims(dataspace_id, dims_b, nullptr);
    m_borderdataPtr = (short *)calloc(dims_b[0] * dims_b[1] * dims_b[2], 2);
    H5Dread(cell_bor_did, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_borderdataPtr);

    hid_t d_id = H5Dopen(file_id, "/cellBin/blockSize", H5P_DEFAULT);
    H5Dread(d_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_block_size);
    H5Dclose(d_id);

    // int min_x, min_y, max_x, max_y;     // replace with global var 2023-6-13
    hid_t attr = H5Aopen(cell_bor_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &cellbin_minx);
    attr = H5Aopen(cell_bor_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &cellbin_miny);
    attr = H5Aopen(cell_bor_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &cellbin_maxx);
    attr = H5Aopen(cell_bor_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &cellbin_maxy);
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", cellbin_minx, cellbin_miny, cellbin_maxx, cellbin_maxy);
    m_effective_rect[0] = cellbin_minx;
    m_effective_rect[1] = cellbin_miny;
    m_effective_rect[2] = cellbin_maxx;
    m_effective_rect[3] = cellbin_maxy;

    H5Sclose(dataspace_id);
    H5Dclose(cell_bor_did);

    // cell type
    hid_t cell_type_did = H5Dopen(file_id, "/cellBin/cellTypeList", H5P_DEFAULT);
    dataspace_id = H5Dget_space(cell_type_did);
    H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
    hid_t str32_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type, 32);
    m_ctypecnt = dims[0];
    m_ctypePtr = new S32[dims[0]];
    H5Dread(cell_type_did, str32_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_ctypePtr);
    H5Tclose(str32_type);
    H5Sclose(dataspace_id);
    H5Dclose(cell_type_did);

    // cell exp
    hid_t cell_exp_did = H5Dopen(file_id, "/cellBin/cellExp", H5P_DEFAULT);
    dataspace_id = H5Dget_space(cell_exp_did);
    H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
    cellexpcnt = dims[0];
    if (isOlderCellExpDataVersion(file_id)) {
        isOldCellExpVersion = true;
        memtype = getMemtypeOfOlderCellExpData();
        m_olderCellExpPtr = (olderCellExpData *)malloc(dims[0] * sizeof(olderCellExpData));
        H5Dread(cell_exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_olderCellExpPtr);
    } else {
        isOldCellExpVersion = false;
        memtype = getMemtypeOfCellExpData();
        m_cellexpPtr = (CellExpData *)malloc(dims[0] * sizeof(CellExpData));
        H5Dread(cell_exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cellexpPtr);
    }
    H5Tclose(memtype);
    H5Sclose(dataspace_id);
    H5Dclose(cell_exp_did);

    // gene
    hid_t gene_did = H5Dopen(file_id, "/cellBin/gene", H5P_DEFAULT);
    dataspace_id = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
    m_genencnt = dims[0];
    memtype = getMemtypeOfGeneData();
    m_genePtr = (GeneData *)malloc(dims[0] * sizeof(GeneData));
    H5Dread(gene_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_genePtr);
    H5Tclose(memtype);
    H5Sclose(dataspace_id);
    H5Dclose(gene_did);

    // gene exp
    //  hid_t gene_exp_did = H5Dopen(file_id, "/cellBin/geneExp", H5P_DEFAULT);
    //  memtype = getMemtypeOfGeneExpData();
    //  dataspace_id = H5Dget_space(gene_exp_did);
    //  H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
    //  m_geneexpcnt = dims[0];
    //  GeneExpData* gene_exp_data = (GeneExpData*)malloc(dims[0] *
    //  sizeof(GeneExpData)); H5Dread(gene_exp_did, memtype, H5S_ALL, H5S_ALL,
    //  H5P_DEFAULT, gene_exp_data); H5Tclose(memtype); H5Sclose(dataspace_id);
    //  H5Dclose(gene_exp_did);

    // 存在exon信息
    if (H5Lexists(file_id, "/cellBin/cellExon", H5P_DEFAULT) > 0) {
        m_bexon = true;
        // cell exon
        hid_t cellexon_did = H5Dopen(file_id, "/cellBin/cellExon", H5P_DEFAULT);
        m_cellexonPtr = (uint16_t *)malloc(m_cellcnt * 2);
        H5Dread(cellexon_did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cellexonPtr);
        H5Dclose(cellexon_did);

        // cellExp exon
        hid_t cellExpexon_did = H5Dopen(file_id, "/cellBin/cellExpExon", H5P_DEFAULT);
        m_cellexonexpPtr = (uint16_t *)malloc(cellexpcnt * 2);
        H5Dread(cellExpexon_did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cellexonexpPtr);
        H5Dclose(cellExpexon_did);

        // gene exon
        //  hid_t geneexon_did = H5Dopen(file_id, "/cellBin/geneExon", H5P_DEFAULT);
        //  uint32_t *pgexon = new uint32_t[m_genencnt];
        //  H5Dread(geneexon_did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
        //  pgexon); H5Dclose(geneexon_did);

        // geneExp exon
        //  hid_t geneExpexon_did = H5Dopen(file_id, "/cellBin/geneExpExon",
        //  H5P_DEFAULT); uint32_t *pgexon_exp = new uint32_t[m_geneexpcnt];
        //  H5Dread(geneExpexon_did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL,
        //  H5P_DEFAULT, pgexon_exp); H5Dclose(geneExpexon_did);
    }

    attr = H5Aopen(file_id, "offsetX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetX);
    attr = H5Aopen(file_id, "offsetY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetY);
    attr = H5Aopen(file_id, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_resolution);
    printf("offsetx:%d offsety:%d\n", m_offsetX, m_offsetY);
    H5Aclose(attr);
    H5Fclose(file_id);
}

void cellAdjust::clear() {
    if (m_cell_arrayptr) {
        free(m_cell_arrayptr);
        m_cell_arrayptr = nullptr;
    }

    if (m_cellexpPtr) {
        free(m_cellexpPtr);
        m_cellexpPtr = nullptr;
    }

    if (m_olderCellExpPtr) {
        free(m_olderCellExpPtr);
        m_olderCellExpPtr = nullptr;
    }

    if (m_genePtr) {
        free(m_genePtr);
        m_genePtr = nullptr;
    }

    if (m_cellexonPtr) {
        free(m_cellexonPtr);
        m_cellexonPtr = nullptr;
    }

    if (m_cellexonexpPtr) {
        free(m_cellexonexpPtr);
        m_cellexonexpPtr = nullptr;
    }
}

void cellAdjust::writeToCgef(const string &outpath) {
    m_cgefwPtr = new CgefWriter();
    m_cgefwPtr->setOutput(outpath);
    CellBinAttr cell_bin_attr = {/*.version = */ 2,
                                 /*.resolution = */ m_resolution,
                                 /*.offsetX = */ m_offsetX,
                                 /* .offsetY = */ m_offsetY, m_szomics};
    m_cgefwPtr->storeAttr(cell_bin_attr);

    writeCellToCgef();
    writeGeneToCgef();
    clear();
    delete m_cgefwPtr;
}

void cellAdjust::writeCellToCgef() {
    timer st(__FUNCTION__);
    vector<uint32_t> vec_cid;
    uint64_t l_id = 0;
    for (uint32_t i = 0; i < m_cellcnt; i++) {
        l_id = m_cell_arrayptr[i].x;
        l_id = (l_id << 32) | m_cell_arrayptr[i].y;
        if (m_setcell.find(l_id) != m_setcell.end()) {
            vec_cid.push_back(m_cell_arrayptr[i].id);
        }
    }

    // vector<uint32_t> vec_cid_tt;
    // int x,y;
    // for(uint32_t i=0;i<m_cellcnt;i++)
    // {
    //     if(m_cell_arrayptr[i].x >= m_minx &&
    //     m_cell_arrayptr[i].x <= m_maxx &&
    //     m_cell_arrayptr[i].y >= m_miny &&
    //     m_cell_arrayptr[i].y <= m_maxy)
    //     {
    //         x = m_cell_arrayptr[i].x - m_minx;
    //         y = m_cell_arrayptr[i].y - m_miny;
    //         if(m_fill_points.at<uchar>(y,x))
    //         {
    //             vec_cid_tt.push_back(m_cell_arrayptr[i].id);
    //         }
    //     }
    // }

    log_info << util::Format("rawcellcnt:{0} newcellcnt:{1}", m_cellcnt, vec_cid.size());
    if (vec_cid.size() == 0) {
        log_error << "lasso region has no infomation. please check. ";
        return;
    }

    uint16_t cell_type_id;

    uint32_t offset = 0;
    vector<uint16_t> vec_cellexon;
    vector<uint16_t> vec_cellexon_exp;
    uint16_t minExon = USHRT_MAX, maxExon = 0, maxExpExon = 0, maxExpmid = 0;

    unsigned int blocknum = m_block_size[2] * m_block_size[3];
    vector<vector<uint32_t>> vec_vec_cell;
    for (uint32_t i = 0; i < blocknum; i++) {
        vector<uint32_t> vectmp;
        vec_vec_cell.emplace_back(std::move(vectmp));
    }
    unsigned int blkid = 0;
    unordered_map<uint32_t, uint32_t> map_gid;  // 记录新旧gid对应
    uint32_t newgid = 0;
    for (uint32_t cid : vec_cid) {
        CellData &cd = m_cell_arrayptr[cid];
        blkid = cd.x / m_block_size[0] + (cd.y / m_block_size[1]) * m_block_size[2];
        vec_vec_cell[blkid].emplace_back(cid);

        for (uint32_t i = 0; i < cd.gene_count; i++) {
            if (isOldCellExpVersion) {
                olderCellExpData *pce = m_olderCellExpPtr + cd.offset;
                if (map_gid.find(pce[i].gene_id) == map_gid.end()) {
                    map_gid.emplace(pce[i].gene_id, newgid++);
                }
            } else {
                CellExpData *pce = m_cellexpPtr + cd.offset;
                if (map_gid.find(pce[i].gene_id) == map_gid.end()) {
                    map_gid.emplace(pce[i].gene_id, newgid++);
                }
            }
        }
    }
    log_info << util::Format("rawgene:{0} newgene:{1}", m_genencnt, map_gid.size());

    vector<uint32_t> vec_blkidx;
    vec_blkidx.reserve(blocknum + 1);
    vector<short> vec_border;
    vec_border.reserve(vec_cid.size() * 2 * BORDERCNT);

    m_cgefwPtr->cell_type_list_.insert(m_cgefwPtr->cell_type_list_.end(), m_ctypePtr, m_ctypePtr + m_ctypecnt);

    uint32_t newcid = 0, offcnt = 0, blkcnt = 0, newctypeid = 0;
    for (vector<uint32_t> &vecC : vec_vec_cell) {
        blkcnt = 0;
        for (uint32_t cid : vecC) {
            short *ptr = m_borderdataPtr + cid * 2 * BORDERCNT;
            vec_border.insert(vec_border.end(), ptr, ptr + 2 * BORDERCNT);

            // m_cgefwPtr->cell_exp_list_.insert(m_cgefwPtr->cell_exp_list_.end(),
            // pce, pce + m_cell_arrayptr[cid].gene_count);

            uint16_t *pexonexp = nullptr;
            if (m_bexon) {
                vec_cellexon.emplace_back(m_cellexonPtr[cid]);
                pexonexp = m_cellexonexpPtr + m_cell_arrayptr[cid].offset;
            }
            if (isOldCellExpVersion) {
                olderCellExpData *pce = m_olderCellExpPtr + m_cell_arrayptr[cid].offset;
                for (uint32_t i = 0; i < m_cell_arrayptr[cid].gene_count; i++) {
                    newgid = map_gid[pce[i].gene_id];
                    m_cgefwPtr->cell_exp_list_.emplace_back(newgid, pce[i].count);
                    if (m_map_genedata.find(newgid) == m_map_genedata.end()) {
                        vector<geneData> tmp;
                        m_map_genedata.emplace(newgid, std::move(tmp));
                    }

                    maxExpmid = std::max(maxExpmid, pce[i].count);
                    if (m_bexon) {
                        vec_cellexon_exp.push_back(pexonexp[i]);
                        minExon = std::min(minExon, pexonexp[i]);
                        maxExon = std::max(maxExon, pexonexp[i]);
                        maxExpExon = std::max(maxExpExon, pexonexp[i]);
                        m_map_genedata[newgid].emplace_back(pexonexp[i], pce[i].count, newcid);
                    } else {
                        m_map_genedata[newgid].emplace_back(0, pce[i].count, newcid);
                    }
                }
            } else {
                CellExpData *pce = m_cellexpPtr + m_cell_arrayptr[cid].offset;
                for (uint32_t i = 0; i < m_cell_arrayptr[cid].gene_count; i++) {
                    newgid = map_gid[pce[i].gene_id];
                    m_cgefwPtr->cell_exp_list_.emplace_back(newgid, pce[i].count);
                    if (m_map_genedata.find(newgid) == m_map_genedata.end()) {
                        vector<geneData> tmp;
                        m_map_genedata.emplace(newgid, std::move(tmp));
                    }

                    maxExpmid = std::max(maxExpmid, pce[i].count);
                    if (m_bexon) {
                        vec_cellexon_exp.push_back(pexonexp[i]);
                        minExon = std::min(minExon, pexonexp[i]);
                        maxExon = std::max(maxExon, pexonexp[i]);
                        maxExpExon = std::max(maxExpExon, pexonexp[i]);
                        m_map_genedata[newgid].emplace_back(pexonexp[i], pce[i].count, newcid);
                    } else {
                        m_map_genedata[newgid].emplace_back(0, pce[i].count, newcid);
                    }
                }
            }

            CellData cell {newcid++,
                           m_cell_arrayptr[cid].x,
                           m_cell_arrayptr[cid].y,
                           offset,
                           m_cell_arrayptr[cid].gene_count,
                           m_cell_arrayptr[cid].exp_count,
                           m_cell_arrayptr[cid].dnb_count,
                           m_cell_arrayptr[cid].area,
                           m_cell_arrayptr[cid].cell_type_id};
            offset += cell.gene_count;

            m_cgefwPtr->cell_attr_.min_x = std::min(m_cgefwPtr->cell_attr_.min_x, cell.x);
            m_cgefwPtr->cell_attr_.max_x = std::max(m_cgefwPtr->cell_attr_.max_x, cell.x);

            m_cgefwPtr->cell_attr_.min_y = std::min(m_cgefwPtr->cell_attr_.min_y, cell.y);
            m_cgefwPtr->cell_attr_.max_y = std::max(m_cgefwPtr->cell_attr_.max_y, cell.y);

            m_cgefwPtr->cell_attr_.min_area = std::min(m_cgefwPtr->cell_attr_.min_area, cell.area);
            m_cgefwPtr->cell_attr_.max_area = std::max(m_cgefwPtr->cell_attr_.max_area, cell.area);

            m_cgefwPtr->cell_attr_.min_gene_count = std::min(m_cgefwPtr->cell_attr_.min_gene_count, cell.gene_count);
            m_cgefwPtr->cell_attr_.max_gene_count = std::max(m_cgefwPtr->cell_attr_.max_gene_count, cell.gene_count);

            m_cgefwPtr->cell_attr_.min_exp_count = std::min(m_cgefwPtr->cell_attr_.min_exp_count, cell.exp_count);
            m_cgefwPtr->cell_attr_.max_exp_count = std::max(m_cgefwPtr->cell_attr_.max_exp_count, cell.exp_count);

            m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, cell.dnb_count);
            m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, cell.dnb_count);

            m_cgefwPtr->expression_num_ += cell.gene_count;
            m_cgefwPtr->exp_count_sum_ += cell.exp_count;
            m_cgefwPtr->dnb_count_sum_ += cell.dnb_count;
            m_cgefwPtr->area_sum_ += cell.area;
            m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
            blkcnt++;
        }
        vec_blkidx.emplace_back(offcnt);
        offcnt += blkcnt;
    }

    vec_blkidx.emplace_back(newcid);
    m_cgefwPtr->cell_num_ = newcid;
    m_cgefwPtr->max_mid_count_ = maxExpmid;

    m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), m_cgefwPtr->cell_num_, m_effective_rect);
    m_cgefwPtr->storeCell(m_block_size[2] * m_block_size[3], vec_blkidx.data(), m_block_size);
    m_cgefwPtr->storeCellExp();
    m_cgefwPtr->storeCellTypeList_N();

    if (m_bexon) {
        m_cgefwPtr->storeCellExon(minExon, maxExon, vec_cellexon, maxExpExon, vec_cellexon_exp);
    }
}

void cellAdjust::writeGeneToCgef() {
    timer st(__FUNCTION__);
    m_cgefwPtr->gene_num_ = m_map_genedata.size();
    GeneData *gene_data_list = static_cast<GeneData *>(calloc(m_cgefwPtr->gene_num_, sizeof(GeneData)));

    vector<GeneExpData> gene_exp_list;
    gene_exp_list.reserve(m_cgefwPtr->gene_num_);

    uint32_t *gene_exon_ptr = (uint32_t *)calloc(m_cgefwPtr->gene_num_, 4);
    vector<uint16_t> vec_exonExp;
    vec_exonExp.reserve(m_cgefwPtr->gene_num_);

    uint16_t maxExpExon = 0;
    uint32_t offset = 0, cellcnt = 0;
    uint32_t minExon = UINT_MAX, maxExon = 0;

    uint32_t minExp = UINT_MAX, maxExp = 0, minCell = UINT_MAX, maxCell = 0;
    int i = 0;
    auto itor = m_map_genedata.begin();
    for (; itor != m_map_genedata.end(); itor++, i++) {
        memcpy(gene_data_list[i].gene_name, m_genePtr[itor->first].gene_name, 64);

        uint16_t maxmid = 0;
        uint32_t expsum = 0, exonsum = 0;
        for (geneData &gd : itor->second) {
            gene_exp_list.emplace_back(gd.cid, gd.mid);
            expsum += gd.mid;
            exonsum += gd.exon;
            maxmid = std::max(maxmid, gd.mid);
            vec_exonExp.emplace_back(gd.exon);
            maxExpExon = std::max(maxExpExon, gd.exon);
        }

        gene_data_list[i].cell_count = itor->second.size();
        gene_data_list[i].exp_count = expsum;
        gene_data_list[i].max_mid_count = maxmid;
        gene_data_list[i].offset = offset;
        offset += itor->second.size();
        m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, maxmid);

        minExp = std::min(minExp, expsum);
        maxExp = std::max(maxExp, expsum);
        minCell = std::min(minCell, cellcnt);
        maxCell = std::max(maxCell, cellcnt);

        minExon = std::min(minExon, exonsum);
        maxExon = std::max(maxExon, exonsum);
    }

    m_cgefwPtr->expression_num_ = gene_exp_list.size();
    m_cgefwPtr->storeGeneAndGeneExp(minExp, maxExp, minCell, maxCell, gene_data_list, gene_exp_list);

    if (m_bexon) {
        m_cgefwPtr->storeGeneExon(minExon, maxExon, gene_exon_ptr, maxExpExon, vec_exonExp);
    }
    free(gene_data_list);
    free(gene_exon_ptr);
}

//////////////////////////////////sap////////////////////////////////////////
void cellAdjust::getSapRegion(const string &strinput, int bin, int thcnt, vector<vector<int>> &vecpos,
                              vector<sapBgefData> &vecdata, float &region_area) {
    timer st(__FUNCTION__);
    m_bgeffile_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    char dataName[32] = {0};
    sprintf(dataName, "/wholeExp/bin%d", bin);

    hsize_t dims[2];
    hid_t gene_did = H5Dopen(m_bgeffile_id, dataName, H5P_DEFAULT);
    if (gene_did < 0) {
        log_error << errorCode::E_MISSINGFILEINFO << "can't find " << dataName;
        exit(-1);
    }
    hid_t gene_sid = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(gene_sid, dims, nullptr);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(BinStat));
    H5Tinsert(memtype, "MIDcount", HOFFSET(BinStat, mid_count), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "genecount", HOFFSET(BinStat, gene_count), H5T_NATIVE_USHORT);

    m_parry = (BinStat *)malloc(dims[0] * dims[1] * sizeof(BinStat));
    H5Dread(gene_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_parry);
    H5Tclose(memtype);

    hid_t attr = H5Aopen(gene_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_x);
    attr = H5Aopen(gene_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_y);
    attr = H5Aopen(gene_did, "lenX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_x);
    attr = H5Aopen(gene_did, "lenY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_y);

    uint32_t resolution = 0;
    attr = H5Aopen(gene_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &resolution);
    log_info << util::Format("minx:{0} miny:{1} lenx:{2} leny:{3} resolution:{4}", m_min_x, m_min_y, m_max_x, m_max_y,
                             resolution);
    H5Aclose(attr);
    H5Sclose(gene_sid);
    H5Dclose(gene_did);

    vector<vector<cv::Point>> vecbors;
    for (vector<int> &vec : vecpos) {
        vector<cv::Point> vtmp;
        vtmp.reserve(vec.size() / 2);
        for (int i = 0; i < vec.size(); i += 2) {
            vtmp.emplace_back(vec[i], vec[i + 1]);
        }
        vecbors.emplace_back(std::move(vtmp));
    }
    cv::Mat fill_points = cv::Mat::zeros(m_max_y, m_max_x, CV_8UC1);
    cv::drawContours(fill_points, vecbors, -1, 1, -1);
    region_area = (float)cv::countNonZero(fill_points) * std::pow((float)resolution / (float)1000, 2);

    if (bin != 1) {
        int id = 0, m, n;
        for (uint32_t x = 0; x < dims[0]; x++) {      // cols
            for (uint32_t y = 0; y < dims[1]; y++) {  // rows
                id = x * dims[1] + y;
                m = x * bin;
                n = y * bin;
                if (fill_points.at<uchar>(n, m)) {
                    if (m_parry[id].gene_count) {
                        vecdata.emplace_back(m_parry[id].gene_count, m_parry[id].mid_count, m, n);
                    }
                }
            }
        }
    } else {
        ThreadPool thpool(thcnt);
        for (int i = 0; i < thcnt; i++) {
            getsapdataTask *ptask = new getsapdataTask(i, thcnt, fill_points, m_parry, vecdata);
            thpool.addTask(ptask);
        }
        thpool.waitTaskDone();
    }
    free(m_parry);
}

void cellAdjust::getSapRegionIndex(const string &strinput, int bin, int thcnt, vector<vector<int>> &vecpos,
                                   vector<vector<int>> &vecdata) {
    timer st(__FUNCTION__);
    m_bgeffile_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    char dataName[32] = {0};
    sprintf(dataName, "/wholeExp/bin%d", bin);

    hsize_t dims[2];
    hid_t gene_did = H5Dopen(m_bgeffile_id, dataName, H5P_DEFAULT);
    if (gene_did < 0) {
        log_error << errorCode::E_MISSINGFILEINFO << "can't find " << dataName;
        exit(-1);
    }
    hid_t gene_sid = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(gene_sid, dims, nullptr);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(BinStat));
    H5Tinsert(memtype, "MIDcount", HOFFSET(BinStat, mid_count), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "genecount", HOFFSET(BinStat, gene_count), H5T_NATIVE_USHORT);

    m_parry = (BinStat *)malloc(dims[0] * dims[1] * sizeof(BinStat));
    H5Dread(gene_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_parry);
    H5Tclose(memtype);

    hid_t attr = H5Aopen(gene_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_x);
    attr = H5Aopen(gene_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_y);
    attr = H5Aopen(gene_did, "lenX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_x);
    attr = H5Aopen(gene_did, "lenY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_y);

    uint32_t resolution = 0;
    attr = H5Aopen(gene_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &resolution);
    log_info << util::Format("minx:{0} miny:{1} lenx:{2} leny:{3} resolution:{4}", m_min_x, m_min_y, m_max_x, m_max_y,
                             resolution);
    H5Aclose(attr);
    H5Sclose(gene_sid);
    H5Dclose(gene_did);

    vector<vector<cv::Point>> vecbors;
    for (vector<int> &vec : vecpos) {
        vector<cv::Point> vtmp;
        vtmp.reserve(vec.size() / 2);
        for (int i = 0; i < vec.size(); i += 2) {
            vtmp.emplace_back(vec[i], vec[i + 1]);
        }
        vecbors.emplace_back(std::move(vtmp));
    }
    cv::Mat fill_points = cv::Mat::zeros(m_max_y, m_max_x, CV_8UC1);
    cv::drawContours(fill_points, vecbors, -1, 1, -1);

    std::vector<int> x_list;
    std::vector<int> y_list;
    // if (bin != 1) {
    int id = 0, m, n;
    for (uint32_t x = 0; x < dims[0]; x++) {      // cols
        for (uint32_t y = 0; y < dims[1]; y++) {  // rows
            id = x * dims[1] + y;
            m = x * bin;
            n = y * bin;
            if (fill_points.at<uchar>(n, m)) {
                if (m_parry[id].gene_count) {
                    x_list.emplace_back(m);
                    y_list.emplace_back(n);
                }
            }
        }
    }
    vecdata.emplace_back(x_list);
    vecdata.emplace_back(y_list);
    // } else {
    //     ThreadPool thpool(thcnt);
    //     for (int i = 0; i < thcnt; i++) {
    //         getsapdataTask *ptask = new getsapdataTask(i, thcnt, fill_points, m_parry, vecdata);
    //         thpool.addTask(ptask);
    //     }
    //     thpool.waitTaskDone();
    // }
    free(m_parry);
}

void cellAdjust::getRegionCelldataSap(vector<vector<int>> &m_vecpos) {
    if (m_vecpos.empty()) {
        std::cout << "No region data input!" << std::endl;
    }

    int num = m_vecpos.size();
    int rows = 0;
    int cols = 0;
    std::vector<cv::Point> relativepoint;
    std::vector<std::vector<cv::Point>> relativepoints;
    std::vector<cv::Point> non_zerovecpoint;

    for (int i = 0; i < num; i++) {
        relativepoint.clear();

        int cnt = m_vecpos[i].size();
        int *pos = m_vecpos[i].data();
        for (int j = 0; j < cnt; j += 2) {
            relativepoint.emplace_back(pos[j] - cellbin_minx, pos[j + 1] - cellbin_miny);
        }

        relativepoints.emplace_back(std::move(relativepoint));
    }

    rows = cellbin_maxy - cellbin_miny + 1;
    cols = cellbin_maxx - cellbin_minx + 1;
    multilabel_img = cv::Mat::zeros(rows, cols, CV_8UC1);
    cv::fillPoly(multilabel_img, relativepoints, 1);
    cv::findNonZero(multilabel_img, non_zerovecpoint);

    lasso_total_area_ = cv::countNonZero(multilabel_img);
}

void cellAdjust::getSapCellbinRegion(sapCgefData &vecdata) {
    std::vector<uint32_t> vec_cid;
    std::cout << "cell count is : " << m_cellcnt << std::endl;
    for (uint32_t i = 0; i < m_cellcnt; i++) {
        if (multilabel_img.at<uchar>((m_cell_arrayptr[i].y - cellbin_miny), (m_cell_arrayptr[i].x - cellbin_minx))) {
            std::cout << "i is : " << i << std::endl;
            vec_cid.push_back(i);
        }
    }

    std::vector<CellData> cell_list_;
    unsigned int expression_num_ = 0;
    unsigned long long int exp_count_sum_ = 0;
    unsigned long long int dnb_count_sum_ = 0;
    unsigned long long int area_sum_ = 0;
    uint32_t offset = 0;
    uint32_t newcid = 0;

    float pixel_ratio = std::pow((float)m_resolution / (float)1000, 2);

    if (vec_cid.empty()) {
        vecdata.total_area = (float)lasso_total_area_ * pixel_ratio;
        vecdata.cell_count = 0;
        vecdata.median_area = 0;
        vecdata.median_gene_count = 0;
        vecdata.median_exp_count = 0;
        vecdata.median_dnb_count = 0;
        vecdata.average_gene_count = 0;
        vecdata.average_exp_count = 0;
        vecdata.average_dnb_count = 0;
        vecdata.average_area = 0;
        return;
    }
    for (uint32_t cid : vec_cid) {
        CellData cell {newcid++,
                       m_cell_arrayptr[cid].x,
                       m_cell_arrayptr[cid].y,
                       offset,
                       m_cell_arrayptr[cid].gene_count,
                       m_cell_arrayptr[cid].exp_count,
                       m_cell_arrayptr[cid].dnb_count,
                       m_cell_arrayptr[cid].area,
                       m_cell_arrayptr[cid].cell_type_id};
        offset += cell.gene_count;

        expression_num_ += cell.gene_count;
        exp_count_sum_ += cell.exp_count;
        dnb_count_sum_ += cell.dnb_count;
        area_sum_ += cell.area;
        cell_list_.emplace_back(std::move(cell));
    }
    vecdata.cell_count = newcid;
    // data.total_area = (float)(area_sum_) * pixel_ratio;
    vecdata.total_area = (float)lasso_total_area_ * pixel_ratio;

    // median
    auto *index = (unsigned int *)malloc(newcid * sizeof(unsigned int));
    std::iota(index, index + newcid, 0);
    std::sort(index, index + newcid,
              [this, &cell_list_](int a, int b) { return cell_list_[a].area < cell_list_[b].area; });

    unsigned int mi = newcid / 2;
    if (newcid % 2 != 0) {
        vecdata.median_area = (float)cell_list_[index[mi]].area * pixel_ratio;
    } else {
        vecdata.median_area = (float(cell_list_[index[mi]].area + cell_list_[index[mi - 1]].area) / 2) * pixel_ratio;
    }

    std::iota(index, index + newcid, 0);
    std::sort(index, index + newcid,
              [this, &cell_list_](int a, int b) { return cell_list_[a].gene_count < cell_list_[b].gene_count; });
    if (newcid % 2 != 0) {
        vecdata.median_gene_count = cell_list_[index[mi]].gene_count;
    } else {
        vecdata.median_gene_count = float(cell_list_[index[mi]].gene_count + cell_list_[index[mi - 1]].gene_count) / 2;
    }

    std::iota(index, index + newcid, 0);
    std::sort(index, index + newcid,
              [this, &cell_list_](int a, int b) { return cell_list_[a].exp_count < cell_list_[b].exp_count; });
    if (newcid % 2 != 0) {
        vecdata.median_exp_count = cell_list_[index[mi]].exp_count;
    } else {
        vecdata.median_exp_count = float(cell_list_[index[mi]].exp_count + cell_list_[index[mi - 1]].exp_count) / 2;
    }

    std::iota(index, index + newcid, 0);
    std::sort(index, index + newcid,
              [this, &cell_list_](int a, int b) { return cell_list_[a].dnb_count < cell_list_[b].dnb_count; });
    if (newcid % 2 != 0) {
        vecdata.median_dnb_count = cell_list_[index[mi]].dnb_count;
    } else {
        vecdata.median_dnb_count = float(cell_list_[index[mi]].dnb_count + cell_list_[index[mi - 1]].dnb_count) / 2;
    }

    // average
    vecdata.average_gene_count = static_cast<float>(expression_num_) / static_cast<float>(newcid);
    vecdata.average_exp_count = static_cast<float>(exp_count_sum_) / static_cast<float>(newcid);
    vecdata.average_dnb_count = static_cast<float>(dnb_count_sum_) / static_cast<float>(newcid);
    vecdata.average_area = (static_cast<float>(area_sum_) / static_cast<float>(newcid)) * pixel_ratio;
}

void cellAdjust::getMultiLabelInfoFromBgef(const string &strinput, vector<vector<int>> &vecpos,
                                           vector<LabelGeneData> &vecdata, uint32_t &total_mid, int bin, int thcnt) {
    timer st(__FUNCTION__);
    hid_t m_bgeffile_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (m_bgeffile_id < 0) {
        std::cout << "can't open spatial bin gef file. "
                  << "\n";
        return;
    }

    // get gene info
    hsize_t dims[1];
    char dataName[32] = {0};
    sprintf(dataName, "/geneExp/bin%d/gene", bin);
    hid_t gene_did = H5Dopen(m_bgeffile_id, dataName, H5P_DEFAULT);
    if (m_bgeffile_id < 0) {
        std::cout << "can't open " << dataName << "dataset. "
                  << "\n";
        return;
    }
    hid_t gene_sid = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(gene_sid, dims, nullptr);

    uint32_t gene_cnt = dims[0];
    Gene *genePtr = (Gene *)malloc(dims[0] * sizeof(Gene));

    hid_t genememtype, strtype;
    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 64);

    genememtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(genememtype, "gene", HOFFSET(Gene, gene), strtype);
    H5Tinsert(genememtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(genememtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);
    H5Dread(gene_did, genememtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, genePtr);
    H5Tclose(genememtype);
    H5Sclose(gene_sid);
    H5Dclose(gene_did);

    // get expression info
    hsize_t exp_dims[1];
    char exp_name[32] = {0};
    sprintf(exp_name, "/geneExp/bin%d/expression", bin);
    hid_t exp_did = H5Dopen(m_bgeffile_id, exp_name, H5P_DEFAULT);
    hid_t exp_sid = H5Dget_space(exp_did);
    H5Sget_simple_extent_dims(exp_sid, exp_dims, nullptr);

    uint64_t geneexp_cnt = exp_dims[0];
    std::cout << "exp count : " << exp_dims[0] << std::endl;
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    Expression *expPtr = (Expression *)calloc(exp_dims[0], sizeof(Expression));
    H5Dread(exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, expPtr);

    // get attribution
    hid_t attr = H5Aopen(exp_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_x);
    attr = H5Aopen(exp_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_y);
    attr = H5Aopen(exp_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_x);
    attr = H5Aopen(exp_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_y);
    attr = H5Aopen(exp_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_resolution);
    std::cout << "minx:" << m_min_x << "miny:" << m_min_y << "maxx:" << m_max_x << "maxy:" << m_max_y << '\n';

    H5Aclose(attr);
    H5Tclose(memtype);
    H5Sclose(exp_sid);
    H5Dclose(exp_did);
    H5Tclose(strtype);

    cv::Mat fill_points;
    {
        timer st("find region... ");
        // find gene in region
        vector<vector<cv::Point>> vecbors;
        for (vector<int> &vec : vecpos) {
            vector<cv::Point> vtmp;
            vtmp.reserve(vec.size() / 2);
            for (int i = 0; i < vec.size(); i += 2) {
                vtmp.emplace_back(vec[i], vec[i + 1]);
            }
            vecbors.emplace_back(std::move(vtmp));
        }
        fill_points = cv::Mat::zeros((m_max_y + 1), (m_max_x + 1), CV_8UC1);
        cv::drawContours(fill_points, vecbors, -1, 1, -1);
    }

    {
        timer st("process cost time... ");
        total_mid = 0;

        ThreadPool thpool(thcnt);
        for (int i = 0; i < thcnt; i++) {
            getLabelInfoTask *ptask =
                new getLabelInfoTask(i, thcnt, fill_points, genePtr, expPtr, vecdata, total_mid, gene_cnt);
            thpool.addTask(ptask);
        }
        thpool.waitTaskDone();

        std::cout << vecdata.size() << std::endl;

        sort(vecdata.begin(), vecdata.end(), [](const LabelGeneData &f1, const LabelGeneData &f2) {
            if (f1.mid_cnt > f2.mid_cnt) {
                return true;
            }
            return false;
        });
        std::cout << "sort... " << std::endl;
    }

    free(genePtr);
    free(expPtr);
}

void cellAdjust::getMultiLabelInfoFromCgef(const string &strcgef, vector<vector<int>> &vecpos,
                                           vector<LabelCellData> &vecdata, vector<LabelCellData> &total_data) {
    timer st(__FUNCTION__);

    // read cgef file
    hid_t file_id = H5Fopen(strcgef.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t cell_did = H5Dopen(file_id, "/cellBin/cell", H5P_DEFAULT);
    hsize_t dims[1];
    hid_t cell_sid = H5Dget_space(cell_did);
    H5Sget_simple_extent_dims(cell_sid, dims, nullptr);
    m_cellcnt = dims[0];
    hid_t memtype = getMemtypeOfCellData();
    m_cell_arrayptr = (CellData *)malloc(dims[0] * sizeof(CellData));
    H5Dread(cell_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cell_arrayptr);

    hid_t attr = H5Aopen(cell_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &cellbin_minx);
    attr = H5Aopen(cell_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &cellbin_miny);
    attr = H5Aopen(cell_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &cellbin_maxx);
    attr = H5Aopen(cell_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &cellbin_maxy);
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", cellbin_minx, cellbin_miny, cellbin_maxx, cellbin_maxy);

    attr = H5Aopen(file_id, "offsetX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetX);
    attr = H5Aopen(file_id, "offsetY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetY);
    attr = H5Aopen(file_id, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_resolution);
    log_info << "offsetx:" << m_offsetX << "offsety:" << m_offsetY << "resolution:" << m_resolution;
    H5Tclose(memtype);
    H5Sclose(cell_sid);
    H5Dclose(cell_did);
    H5Aclose(attr);
    H5Fclose(file_id);

    // get region
    if (vecpos.empty()) {
        log_info << "No region data input!";
    }

    int num = vecpos.size();
    std::vector<cv::Point> non_zerovecpoint;
    std::vector<cv::Point> relativepoint;
    int rows = 0;
    int cols = 0;
    std::vector<std::vector<cv::Point>> relativepoints;

    for (int i = 0; i < num; i++) {
        relativepoint.clear();

        int cnt = vecpos[i].size();
        int *pos = vecpos[i].data();
        for (int j = 0; j < cnt; j += 2) {
            relativepoint.emplace_back(pos[j] - cellbin_minx, pos[j + 1] - cellbin_miny);
        }

        relativepoints.emplace_back(std::move(relativepoint));
    }

    rows = cellbin_maxy - cellbin_miny + 1;
    cols = cellbin_maxx - cellbin_minx + 1;
    multilabel_img = cv::Mat::zeros(rows, cols, CV_8UC1);
    cv::fillPoly(multilabel_img, relativepoints, 1);
    cv::findNonZero(multilabel_img, non_zerovecpoint);

    lasso_total_area_ = cv::countNonZero(multilabel_img);

    // get statistics datas
    std::vector<uint32_t> vec_cid;

    for (uint32_t i = 0; i < m_cellcnt; i++) {
        if (multilabel_img.at<uchar>((m_cell_arrayptr[i].y - cellbin_miny), (m_cell_arrayptr[i].x - cellbin_minx))) {
            vec_cid.push_back(i);
        }
    }

    if (vec_cid.empty()) {
        std::cout << "lasso region is empty. ";
        return;
    }

    float pixel_ratio = std::pow((float)m_resolution / (float)1000, 2);
    std::map<uint16_t, LabelCellData> total_cell_data;

    for (uint32_t cid : vec_cid) {
        unsigned short cls_id = m_cell_arrayptr[cid].cluster_id;

        if (total_cell_data.find(cls_id) != total_cell_data.end()) {
            total_cell_data[cls_id].cluster_id = m_cell_arrayptr[cid].cluster_id;
            total_cell_data[cls_id].cell_id++;
            total_cell_data[cls_id].area += m_cell_arrayptr[cid].area;
            total_cell_data[cls_id].mid_cnt += m_cell_arrayptr[cid].exp_count;
        } else {
            LabelCellData tmp = {m_cell_arrayptr[cid].cluster_id, m_cell_arrayptr[cid].exp_count,
                                 (float)m_cell_arrayptr[cid].area, 1};

            total_cell_data.insert(std::pair<uint16_t, LabelCellData>(cls_id, tmp));
        }

        float cell_area = (float)m_cell_arrayptr[cid].area * pixel_ratio;
        vecdata.emplace_back(m_cell_arrayptr[cid].cluster_id, m_cell_arrayptr[cid].exp_count, cell_area,
                             m_cell_arrayptr[cid].id);
    }
    for (auto itor : total_cell_data) {
        total_data.emplace_back(itor.second.cluster_id, itor.second.mid_cnt, ((float)itor.second.area * pixel_ratio),
                                itor.second.cell_id);
    }

    std::cout << "start sort ... " << std::endl;
    sort(vecdata.begin(), vecdata.end(), [](const LabelCellData &foo1, const LabelCellData &foo2) {
        if (foo1.cluster_id < foo2.cluster_id) {
            return true;
        } else if (foo1.cluster_id == foo2.cluster_id) {
            if (foo1.cell_id < foo2.cell_id) {
                return true;
            } else {
                return false;
            }
        }
        return false;
    });
}

void cellAdjust::GetPositionIndexByClusterId(const char *input_file, std::vector<int> cls_id,
                                             std::vector<std::vector<int>> &clusterpos_list) {
    timer st(__FUNCTION__);

    hid_t anndata_fileId = H5Fopen(input_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (anndata_fileId < 0) {
        std::cout << "can't open " << input_file << "\n";
        return;
    }

    // get cluster id
    char dataName[32] = {0};
    sprintf(dataName, "/obs/leiden");
    hid_t cluster_did = H5Dopen(anndata_fileId, dataName, H5P_DEFAULT);
    if (cluster_did < 0) {
        std::cout << "can't find " << dataName << "\n";
        return;
    }
    hsize_t dims[1];
    hid_t cluster_sid = H5Dget_space(cluster_did);
    H5Sget_simple_extent_dims(cluster_sid, dims, nullptr);

    int *m_parry = (int *)malloc(dims[0] * sizeof(int));
    H5Dread(cluster_did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_parry);

    vector<int> v(m_parry, m_parry + dims[0]);

    H5Sclose(cluster_sid);
    H5Dclose(cluster_did);

    // get spatial position
    sprintf(dataName, "/obsm/spatial");
    hid_t spatial_did = H5Dopen(anndata_fileId, dataName, H5P_DEFAULT);
    if (spatial_did < 0) {
        std::cout << "can't find " << dataName << "\n";
        return;
    }
    hsize_t spatial_dims[2];
    hid_t spatial_sid = H5Dget_space(spatial_did);
    H5Sget_simple_extent_dims(spatial_sid, spatial_dims, nullptr);

    Pos *pos_data = (Pos *)malloc(spatial_dims[0] * sizeof(Pos));
    H5Dread(spatial_did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_data);

    H5Sclose(spatial_sid);
    H5Dclose(spatial_did);
    // close file
    H5Fclose(anndata_fileId);

    // std::vector<ClusterPosition> clusterpos_list;
    std::vector<int> x_cod;
    std::vector<int> y_cod;
    for (int i = 0; i < cls_id.size(); i++) {
        for (int j = 0; j < dims[0]; j++) {
            if (m_parry[j] == cls_id[i]) {
                x_cod.emplace_back(static_cast<int>(pos_data[j].x));
                y_cod.emplace_back(static_cast<int>(pos_data[j].y));
            }
        }
    }
    clusterpos_list.push_back(x_cod);
    clusterpos_list.push_back(y_cod);
}

int cellAdjust::GenerateFilterBgefFileByMidCount(const std::string input_file, const std::string output_file,
                                                 int bin_size, std::vector<MidCntFilter> filter_genes,
                                                 bool only_filter) {
    // check input param
    if (filter_genes.empty()) {
        log_info << "input filter genes is empty. ";
        return -1;
    }

    hid_t file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        log_info << "open bgef file error. ";
        return -1;
    }
    if (H5Lexists(file_id, util::Format("/geneExp/bin{0}", bin_size).c_str(), H5P_DEFAULT) < 0) {
        H5Fclose(file_id);
        log_info << "open /geneExp/bin" << bin_size << " failed. ";
        return -1;
    }
    H5Fclose(file_id);
    // do generate
    process_rate_ = 0;
    log_info << filter_genes.size();
    BgefOptions::GetInstance()->input_file_ = input_file;
    BgefOptions::GetInstance()->output_file_ = output_file;

    // generate_bgef_thread_ = std::thread(&cellAdjust::DoGenerate, this, bin_size, filter_genes, only_filter);
    DoGenerate(bin_size, filter_genes, only_filter);
    return 0;
}

int cellAdjust::GenerateFilterBgefDuration() {
    // get rate of progress
    if ((process_rate_ == -1) || (process_rate_ == 3)) {
        if (generate_bgef_thread_.joinable()) generate_bgef_thread_.join();
        return process_rate_;
    }
    return process_rate_;
}

void cellAdjust::FilterGeneInfo(int bin_size, std::vector<MidCntFilter> filter_genes,
                                std::map<std::string, std::set<uint64_t>> &filter_data, bool only_filter) {
    timer st(__FUNCTION__);

    std::unordered_map<std::string, MidCntFilter> filter_gene_map;
    for (auto itor : filter_genes) {
        filter_gene_map.insert(std::unordered_map<std::string, MidCntFilter>::value_type(itor.gene_name, itor));
    }

    if (1 == bin_size) {
        m_bgefopts = BgefOptions::GetInstance();
        m_bgefopts->map_gene_exp_.clear();

        BgefReader bgef_reader(BgefOptions::GetInstance()->input_file_, 1);
        m_bexon = bgef_reader.isExonExist();
        log_info << "before filter size is : " << bgef_reader.getExpressionNum();

        Gene *gene = bgef_reader.getGene();
        Expression *expression = bgef_reader.getExpression();

        if (only_filter) {
            for (unsigned int gene_id = 0; gene_id < bgef_reader.getGeneNum(); gene_id++) {
                if (filter_gene_map.find(gene[gene_id].gene) != filter_gene_map.end()) {
                    unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                    vector<Expression> exps;
                    for (unsigned int i = gene[gene_id].offset; i < end; i++) {
                        if ((expression[i].count >= filter_gene_map[gene[gene_id].gene].min_mid) ||
                            (expression[i].count <= filter_gene_map[gene[gene_id].gene].max_mid)) {
                            exps.emplace_back(expression[i]);
                        }
                    }
                    if (!exps.empty()) {
                        m_bgefopts->map_gene_exp_.insert(
                            unordered_map<string, vector<Expression>>::value_type(gene[gene_id].gene, exps));
                    }
                }
            }
        } else {
            for (unsigned int gene_id = 0; gene_id < bgef_reader.getGeneNum(); gene_id++) {
                if (filter_gene_map.find(gene[gene_id].gene) != filter_gene_map.end()) {
                    unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                    vector<Expression> exps;
                    for (unsigned int i = gene[gene_id].offset; i < end; i++) {
                        if ((expression[i].count >= filter_gene_map[gene[gene_id].gene].min_mid) ||
                            (expression[i].count <= filter_gene_map[gene[gene_id].gene].max_mid)) {
                            exps.emplace_back(expression[i]);
                        }
                    }
                    if (!exps.empty()) {
                        m_bgefopts->map_gene_exp_.insert(
                            unordered_map<string, vector<Expression>>::value_type(gene[gene_id].gene, exps));
                    }
                } else {
                    vector<Expression> exps;
                    exps.reserve(gene[gene_id].count);
                    unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                    for (unsigned int i = gene[gene_id].offset; i < end; i++) {
                        exps.emplace_back(expression[i]);
                    }
                    m_bgefopts->map_gene_exp_.insert(
                        unordered_map<string, vector<Expression>>::value_type(gene[gene_id].gene, exps));
                }
            }
        }

    } else {
        // build filter map
        std::unordered_map<std::string, MidCntFilter> filter_gene_map;
        for (auto itor : filter_genes) {
            filter_gene_map.insert(std::unordered_map<std::string, MidCntFilter>::value_type(itor.gene_name, itor));
        }

        // begin filter
        BgefReader bgef_reader(BgefOptions::GetInstance()->input_file_, bin_size);
        Gene *gene = bgef_reader.getGene();
        Expression *expression = bgef_reader.getExpression();
        m_bexon = bgef_reader.isExonExist();

        // 这个filter_data可能很大，性能热点
        for (unsigned int gene_id = 0; gene_id < bgef_reader.getGeneNum(); gene_id++) {
            if (filter_gene_map.find(gene[gene_id].gene) != filter_gene_map.end()) {
                std::set<uint64_t> tmp;
                uint64_t idx;
                unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                for (unsigned int i = gene[gene_id].offset; i < end; i++) {
                    if ((expression[i].count < filter_gene_map[gene[gene_id].gene].min_mid) ||
                        (expression[i].count > filter_gene_map[gene[gene_id].gene].max_mid)) {
                        idx = expression[i].x;
                        tmp.emplace((idx << 32 | expression[i].y));
                    }
                }
                if (!tmp.empty()) {
                    filter_data.insert(std::map<std::string, std::set<uint64_t>>::value_type(gene[gene_id].gene, tmp));
                }
            }
        }
    }
}

void cellAdjust::DoGenerate(int bin_size, std::vector<MidCntFilter> filter_genes, bool only_filter) {
    // get gene info at bin1 and binx
    timer st(__FUNCTION__);

    std::map<std::string, std::set<uint64_t>> filter_data;
    FilterGeneInfo(bin_size, filter_genes, filter_data, only_filter);
    log_info << "get filter info finish. ";

    process_rate_ = 1;
    m_bgefopts = BgefOptions::GetInstance();

    // filter bin1 info
    if (bin_size != 1) {
        if (filter_data.empty()) {
            log_info << "error filter gene is empty. ";
            process_rate_ = -1;
            return;
        }

        if (only_filter) {
            BgefReader bgef_reader(BgefOptions::GetInstance()->input_file_, 1);
            m_bexon = bgef_reader.isExonExist();
            log_info << "before filter size is : " << bgef_reader.getExpressionNum();
            Gene *gene = bgef_reader.getGene();
            Expression *expression = bgef_reader.getExpression();
            for (unsigned int gene_id = 0; gene_id < bgef_reader.getGeneNum(); gene_id++) {
                if (filter_data.find(gene[gene_id].gene) != filter_data.end()) {
                    unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                    vector<Expression> exps;
                    for (unsigned int i = gene[gene_id].offset; i < end; i++) {
                        uint64_t idx = (expression[i].x / bin_size) * bin_size;
                        idx = idx << 32 | ((expression[i].y / bin_size) * bin_size);
                        if (filter_data[gene[gene_id].gene].find(idx) != filter_data[gene[gene_id].gene].end()) {
                            continue;
                        }
                        exps.emplace_back(expression[i]);
                    }
                    if (!exps.empty()) {
                        m_bgefopts->map_gene_exp_.insert(
                            unordered_map<string, vector<Expression>>::value_type(gene[gene_id].gene, exps));
                    }
                }
            }
        } else {
            BgefReader bgef_reader(BgefOptions::GetInstance()->input_file_, 1);
            m_bexon = bgef_reader.isExonExist();
            log_info << "before filter size is : " << bgef_reader.getExpressionNum();
            Gene *gene = bgef_reader.getGene();
            Expression *expression = bgef_reader.getExpression();
            for (unsigned int gene_id = 0; gene_id < bgef_reader.getGeneNum(); gene_id++) {
                if (filter_data.find(gene[gene_id].gene) != filter_data.end()) {
                    unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                    vector<Expression> exps;
                    for (unsigned int i = gene[gene_id].offset; i < end; i++) {
                        uint64_t idx = (expression[i].x / bin_size) * bin_size;
                        idx = idx << 32 | ((expression[i].y / bin_size) * bin_size);
                        if (filter_data[gene[gene_id].gene].find(idx) != filter_data[gene[gene_id].gene].end()) {
                            continue;
                        }
                        exps.emplace_back(expression[i]);
                    }
                    if (!exps.empty()) {
                        m_bgefopts->map_gene_exp_.insert(
                            unordered_map<string, vector<Expression>>::value_type(gene[gene_id].gene, exps));
                    }
                } else {
                    vector<Expression> exps;
                    exps.reserve(gene[gene_id].count);
                    unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                    for (unsigned int i = gene[gene_id].offset; i < end; i++) {
                        exps.emplace_back(expression[i]);
                    }
                    m_bgefopts->map_gene_exp_.insert(
                        unordered_map<string, vector<Expression>>::value_type(gene[gene_id].gene, exps));
                }
            }
        }
    }

    // 4. generate bgef file according bin1 info
    if (m_bgefopts->map_gene_exp_.size() == 0) {
        log_error << "region has no information. please check. ";
        process_rate_ = -1;
        return;
    }

    hid_t file_id = H5Fopen(BgefOptions::GetInstance()->input_file_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file_id, "/geneExp", H5P_DEFAULT);
    if (gid < 0) {
        log_info << "can not find input spatial bin gef file... ";
        process_rate_ = -1;
        return;
    }
    std::vector<std::string> group_names;
    herr_t idx = H5Literate(gid, H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, &group_names);
    H5Gclose(gid);
    m_bgefopts->bin_sizes_.clear();
    for (string &str : group_names) {
        int bin = std::stoi(str.substr(3));
        m_bgefopts->bin_sizes_.push_back(bin);
    }

    if (H5Aexists(file_id, "omics") > 0) {
        hid_t omics_strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(omics_strtype, 32);
        hid_t fattr = H5Aopen(file_id, "omics", H5P_DEFAULT);
        H5Aread(fattr, omics_strtype, m_szomics);
        H5Tclose(omics_strtype);
    } else {
        strcpy(m_szomics, "Transcriptomics");
    }
    hid_t exp_did = H5Dopen(file_id, "/geneExp/bin1/expression", H5P_DEFAULT);
    hid_t attr = H5Aopen(exp_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_x);
    attr = H5Aopen(exp_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_y);
    attr = H5Aopen(exp_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_x);
    attr = H5Aopen(exp_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_y);
    attr = H5Aopen(exp_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_resolution);
    log_info << util::Format("minx:{0} miny:{1} maxx:{2} maxy:{3}", m_min_x, m_min_y, m_max_x, m_max_y);
    H5Aclose(attr);
    H5Dclose(exp_did);
    H5Fclose(file_id);

    m_bgefopts->m_genes_queue.init(m_bgefopts->map_gene_exp_.size());
    ThreadPool thpool(m_bgefopts->thread_ * 2);
    m_bgefopts->m_stromics.append(m_szomics);
    log_info << "thread count is : " << m_bgefopts->thread_;

    BgefWriter bgef_writer(BgefOptions::GetInstance()->output_file_, false, m_bexon, m_bgefopts->m_stromics);
    bgef_writer.setResolution(m_resolution);
    // bgef_writer.SetGefArea(gef_area);

    // do bin , write to bgef file
    int genecnt = 0;
    for (unsigned int bin : m_bgefopts->bin_sizes_) {
        timer st("Do bin... ");
        auto &dnb_matrix = m_bgefopts->dnbmatrix_;
        auto &dnbAttr = m_bgefopts->dnbmatrix_.dnb_attr;

        dnbAttr.min_x = (m_min_x / bin) * bin;
        dnbAttr.len_x = (m_max_x) / bin + 1;
        dnbAttr.min_y = (m_min_y / bin) * bin;
        dnbAttr.len_y = (m_max_y) / bin + 1;

        dnbAttr.max_x = (m_max_x / bin) * bin;
        dnbAttr.max_y = (m_max_y / bin) * bin;

        dnbAttr.max_gene = 0;
        dnbAttr.max_mid = 0;
        dnbAttr.number = 0;
        unsigned long matrix_len = (unsigned long)(dnbAttr.len_x) * dnbAttr.len_y;

        log_info << "bin " << bin << " matrix: min_x=" << dnbAttr.min_x << " len_x=" << dnbAttr.len_x
                 << " min_y=" << dnbAttr.min_y << " len_y=" << dnbAttr.len_y << " matrix_len=" << matrix_len;

        if (bin == 1) {
            dnb_matrix.pmatrix_us = (BinStatUS *)calloc(matrix_len, sizeof(BinStatUS));
            assert(dnb_matrix.pmatrix_us);
            if (!dnb_matrix.pmatrix_us) {
                log_error << errorCode::E_ALLOCMEMORYFAILED << "can not alloc memory for wholeExp matrix. ";
                process_rate_ = -1;
                return;
            }
            if (m_bexon) {
                dnb_matrix.pexon16 = (unsigned short *)calloc(matrix_len, 2);
                assert(dnb_matrix.pexon16);
                if (!dnb_matrix.pexon16) {
                    log_error << errorCode::E_ALLOCMEMORYFAILED << "can not alloc memory for wholeExp matrix. ";
                    process_rate_ = -1;
                    return;
                }
            }
        } else {
            dnb_matrix.pmatrix = (BinStat *)calloc(matrix_len, sizeof(BinStat));
            assert(dnb_matrix.pmatrix);
            if (!dnb_matrix.pmatrix) {
                log_error << errorCode::E_ALLOCMEMORYFAILED << "can not alloc memory for wholeExp matrix. ";
                process_rate_ = -1;
                return;
            }
            if (m_bexon) {
                dnb_matrix.pexon32 = (unsigned int *)calloc(matrix_len, 4);
                assert(dnb_matrix.pexon32);
                if (!dnb_matrix.pexon32) {
                    log_error << errorCode::E_ALLOCMEMORYFAILED << "can not alloc memory for wholeExp matrix. ";
                    process_rate_ = -1;
                    return;
                }
            }
        }

        process_rate_ = 2;

        for (int i = 0; i < m_bgefopts->thread_; i++) {
            auto *task = new DnbMergeTask(m_bgefopts->map_gene_exp_.size(), i, bin);
            thpool.addTask(task);
        }

        auto itor = m_bgefopts->map_gene_exp_.begin();
        for (; itor != m_bgefopts->map_gene_exp_.end(); itor++) {
            auto *task = new BinTask(bin, itor->first.c_str());
            thpool.addTask(task);
        }

        unsigned int offset = 0;
        unsigned int maxexp = 0;
        unsigned int maxexon = 0;
        genecnt = 0;
        map<string, vector<Expression>> gene_info;
        while (true) {
            GeneInfo *pgeneinfo = m_bgefopts->m_geneinfo_queue.getPtr();
            if (1 != bin) {
                gene_info.insert(map<string, vector<Expression>>::value_type(pgeneinfo->geneid, *pgeneinfo->vecptr));
            }

            maxexp = std::max(maxexp, pgeneinfo->maxexp);
            maxexon = std::max(maxexon, pgeneinfo->maxexon);

            if (bin == 100) {
                m_bgefopts->m_vec_bin100.emplace_back(pgeneinfo->geneid, pgeneinfo->umicnt, pgeneinfo->e10);
            }

            delete pgeneinfo;
            genecnt++;
            if (genecnt == m_bgefopts->map_gene_exp_.size()) {
                break;
            }
        }

        if (bin == 1) {
            for (auto itor : m_bgefopts->map_gene_exp_) {
                m_bgefopts->expressions_.insert(m_bgefopts->expressions_.end(), itor.second.begin(), itor.second.end());
                m_bgefopts->genes_.emplace_back(itor.first.c_str(), offset,
                                                static_cast<unsigned int>(itor.second.size()));
                offset += itor.second.size();
            }
        } else {
            for (auto itor : gene_info) {
                for (auto g : itor.second) {
                    g.x *= bin;
                    g.y *= bin;
                    m_bgefopts->expressions_.push_back(std::move(g));
                }

                m_bgefopts->genes_.emplace_back(itor.first.c_str(), offset,
                                                static_cast<unsigned int>(itor.second.size()));
                offset += itor.second.size();
            }
        }

        log_info << "after filter size is : " << m_bgefopts->expressions_.size();
        bgef_writer.storeGene(m_bgefopts->expressions_, m_bgefopts->genes_, dnb_matrix.dnb_attr, maxexp, bin);
        bgef_writer.storeGeneExon(m_bgefopts->expressions_, maxexon, bin);
        m_bgefopts->expressions_.clear();
        m_bgefopts->genes_.clear();
        std::vector<Expression>().swap(m_bgefopts->expressions_);
        std::vector<Gene>().swap(m_bgefopts->genes_);

        thpool.waitTaskDone();
        m_bgefopts->m_genes_queue.clear(bin);

        // write dnb
        {
            timer st("write dnb... ");

            if (bin == 100) {
                vector<GeneStat> &geneStat = m_bgefopts->m_vec_bin100;
                std::sort(geneStat.begin(), geneStat.end(), [](const GeneStat &p1, const GeneStat &p2) {
                    if (p1.mid_count > p2.mid_count) {
                        return true;
                    } else if (p1.mid_count == p2.mid_count) {
                        int ret = strcmp(p1.gene, p2.gene);
                        return ret < 0;
                    } else {
                        return false;
                    }
                });
                bgef_writer.storeStat(geneStat);
            }

            vector<unsigned int> vec_mid;
            unsigned long number = 0;

            if (bin == 1) {
                for (unsigned long i = 0; i < matrix_len; i++) {
                    if (dnb_matrix.pmatrix_us[i].gene_count) {
                        ++number;
                        vec_mid.push_back(dnb_matrix.pmatrix_us[i].mid_count);
                    }
                }
            } else {
                for (unsigned long i = 0; i < matrix_len; i++) {
                    if (dnb_matrix.pmatrix[i].gene_count) {
                        ++number;
                        vec_mid.push_back(dnb_matrix.pmatrix[i].mid_count);
                    }
                }
            }

            int sz = vec_mid.size();
            sort(vec_mid.begin(), vec_mid.end(), [](const unsigned int a, const unsigned int b) { return a < b; });
            if (bin > 50) {
                dnbAttr.max_mid = vec_mid[sz - 1];
            } else {
                int limit = sz * 0.999;
                dnbAttr.max_mid = vec_mid[limit];
            }

            dnbAttr.number = number;
            bgef_writer.storeDnb(dnb_matrix, bin);
            bgef_writer.storeWholeExon(dnb_matrix, bin);
        }

        if (bin == 1) {
            if (dnb_matrix.pmatrix_us != nullptr) {
                free(dnb_matrix.pmatrix_us);
                dnb_matrix.pmatrix_us = nullptr;
                if (m_bexon) {
                    free(dnb_matrix.pexon16);
                    dnb_matrix.pexon16 = nullptr;
                }
            }
        } else {
            if (dnb_matrix.pmatrix != nullptr) {
                free(dnb_matrix.pmatrix);
                dnb_matrix.pmatrix = nullptr;
                if (m_bexon) {
                    free(dnb_matrix.pexon32);
                    dnb_matrix.pexon32 = nullptr;
                }
            }
        }
    }
    process_rate_ = 3;
    return;
}