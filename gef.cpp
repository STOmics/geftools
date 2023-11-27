/*
 * @Author: zhaozijian
 * @Date: 2022-02-10 14:53:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-29 15:35:32
 * @Description: file content
 */
#include "gef.h"

#include "utils.h"

hid_t getMemtypeOfGeneData() {
    hid_t memtype;
    // hid_t str32_type_ = H5Tcopy(H5T_C_S1);
    // H5Tset_size(str32_type_, 32);
    hid_t str64_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str64_type_, 64);
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(GeneData));
    H5Tinsert(memtype, "geneName", HOFFSET(GeneData, gene_name), str64_type_);
    H5Tinsert(memtype, "offset", HOFFSET(GeneData, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "cellCount", HOFFSET(GeneData, cell_count), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "expCount", HOFFSET(GeneData, exp_count), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "maxMIDcount", HOFFSET(GeneData, max_mid_count), H5T_NATIVE_USHORT);
    return memtype;
}

hid_t getMemtypeOfGeneExpData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(GeneExpData));
    H5Tinsert(memtype, "cellID", HOFFSET(GeneExpData, cell_id), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(GeneExpData, count), H5T_NATIVE_USHORT);
    // H5Tinsert(memtype, "incnt", HOFFSET(GeneExpData, incnt), H5T_NATIVE_USHORT);
    return memtype;
}

hid_t getMemtypeOfCellData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellData));
    H5Tinsert(memtype, "id", HOFFSET(CellData, id), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "x", HOFFSET(CellData, x), H5T_NATIVE_INT);
    H5Tinsert(memtype, "y", HOFFSET(CellData, y), H5T_NATIVE_INT);
    H5Tinsert(memtype, "offset", HOFFSET(CellData, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "geneCount", HOFFSET(CellData, gene_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "expCount", HOFFSET(CellData, exp_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "dnbCount", HOFFSET(CellData, dnb_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "area", HOFFSET(CellData, area), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "cellTypeID", HOFFSET(CellData, cell_type_id), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "clusterID", HOFFSET(CellData, cluster_id), H5T_NATIVE_USHORT);
    // H5Tinsert(memtype, "incnt", HOFFSET(CellData, incnt), H5T_NATIVE_USHORT);

    return memtype;
}

hid_t getMemtypeOfCellExpData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(CellExpData, gene_id), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(CellExpData, count), H5T_NATIVE_USHORT);
    return memtype;
}

hid_t getMemtypeOfOlderCellExpData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(olderCellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(olderCellExpData, gene_id), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "count", HOFFSET(olderCellExpData, count), H5T_NATIVE_USHORT);
    return memtype;
}

bool isOlderCellExpDataVersion(hid_t fileId) {
    unsigned int geftoolVersion[3] = {0};
    if (H5Aexists(fileId, "geftool_ver") > 0) {
        hid_t attr = H5Aopen(fileId, "geftool_ver", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_UINT32, geftoolVersion);
        log_info << util::Format("version is {0}.{1}.{2} ", geftoolVersion[0], geftoolVersion[1], geftoolVersion[2]);
        H5Aclose(attr);
        
        if (geftoolVersion[0] == 0 && geftoolVersion[1] <= 7) {
            if (geftoolVersion[1] == 7 && geftoolVersion[2] >= 6) {
                return false;
            } else {
                return true;
            }
        } else {
            return false;
        }
    } else {
        return true;
    }
}

std::string getOmicsName(hid_t file_id) {
    std::string omics_type {""};
    std::string omics_name {""};
    if (H5Aexists(file_id, "omics") > 0) {
        hid_t f_attr = H5Aopen(file_id, "omics", H5P_DEFAULT);
        char szbuf[128] = {0};
        hid_t omics_strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(omics_strtype, 32);
        H5Aread(f_attr, omics_strtype, szbuf);
        omics_type.append(szbuf);
        H5Aclose(f_attr);
        H5Tclose(omics_strtype);
    } else {
        log_info << "can not find omics type from file. using default type: Transcriptomics. ";
        omics_name = "gene";
        return omics_name;
    }
    if (omics_type == "Transcriptomics") {
        omics_name = "gene";
    } else {
        omics_name = "protein";
    }
    return omics_name;
}

bool ParseOmicsType(const std::string &bgef_file, std::string &input_omics) {
    hid_t file_id = H5Fopen(bgef_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        log_info << "open bgef file error. ";
        return false;
    }
    std::string omics_type {""};
    if (H5Aexists(file_id, "omics") > 0) {
        hid_t f_attr = H5Aopen(file_id, "omics", H5P_DEFAULT);
        char szbuf[128] = {0};
        hid_t omics_strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(omics_strtype, 32);
        H5Aread(f_attr, omics_strtype, szbuf);
        omics_type.append(szbuf);
        if (omics_type != input_omics) {
            log_info << "\'-O\' information does not match the omics recorded in " << bgef_file
                     << ",please check input parameter or files. ";
            H5Aclose(f_attr);
            H5Tclose(omics_strtype);
            H5Fclose(file_id);
            return false;
        }
        H5Aclose(f_attr);
        H5Tclose(omics_strtype);
    } else {
        log_info << "can not find omics type from file. using default type: Transcriptomics. ";
        omics_type = "Transcriptomics";
        if (omics_type != input_omics) {
            log_info << "\'-O\' information does not match the omics recorded in " << bgef_file
                     << ",please check input parameter or files. ";
            H5Fclose(file_id);
            return false;
        }
    }
    H5Fclose(file_id);
    return true;
}

std::string getOmicsType(const std::string &file_path, const std::string &input_omics) {
    std::string omics {""};

    hid_t file_id = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        log_error << errorCode::E_FILEOPENERROR << "open bgef file error. ";
        return omics;
    }
    if (H5Aexists(file_id, "omics") > 0) {
        hid_t f_attr = H5Aopen(file_id, "omics", H5P_DEFAULT);
        char szbuf[128] = {0};
        hid_t omics_strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(omics_strtype, 32);
        H5Aread(f_attr, omics_strtype, szbuf);
        omics.append(szbuf);
        if (omics != input_omics) {
            log_error << errorCode::E_INVALIDPARAM << "information does not match the omics recorded in " << file_path
                      << ",please check input parameter or files. ";
            H5Aclose(f_attr);
            H5Tclose(omics_strtype);
            H5Fclose(file_id);
            return "";
        }
        H5Aclose(f_attr);
        H5Tclose(omics_strtype);
        H5Fclose(file_id);
        return omics;
    } else {
        if (input_omics == "Transcriptomics") {
            log_info << "can not find omics type from file. using default type: Transcriptomics. ";
            omics = "Transcriptomics";
            H5Fclose(file_id);
            return omics;
        }

        log_error << errorCode::E_INVALIDPARAM << " can not find omics type from file. ";
        H5Fclose(file_id);
        return omics;
    }
}