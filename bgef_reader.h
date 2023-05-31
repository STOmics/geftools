/** @file common_bin.h
    @brief Declare a BgefReader class for reading common bin GEF.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__COMMON_BIN_H_
#define GEFTOOLS__COMMON_BIN_H_

#include <numeric>
#include <vector>
#include <map>
#include "hdf5.h"
#include "utils.h"
#include "gef.h"
#include "bgef_options.h"

class GEFTOOLS_API BgefReader {
  private:
    int bin_size_ = 0;
    unsigned int gene_num_ = 0;
    unsigned int cell_num_ = 0;
    vector<Coordinate> cell_pos_;
    unsigned int * cell_indices_ = nullptr;
    unsigned long long expression_num_ = 0;
    ExpressionAttr expression_attr_{};
    bool expression_attr_init_ = false;
    unsigned int whole_exp_matrix_shape_[2] = {0};
    Gene* genes_ = nullptr;
    Expression* expressions_ = nullptr;
    Expression* reduce_expressions_ = nullptr;
    cv::Mat whole_exp_matrix_t_;
    int version_{};
    int verbose_ = true;
    int n_thread_ = 1;
    BgefOptions *opts_ = nullptr;
    unsigned int *m_exonPtr = nullptr;
    unsigned int max_exon_;

    hid_t file_id_;
    hid_t exp_dataspace_id_{};
    hid_t exp_dataset_id_{};
    hid_t gene_dataspace_id_{};
    hid_t gene_dataset_id_{};
    hid_t whole_exp_dataspace_id_{0};
    hid_t whole_exp_dataset_id_{0};
    hid_t m_exon_did = 0;

    void openExpressionSpace(int bin_size);
    void openGeneSpace(int bin_size);
    void openWholeExpSpace();
    
    void buildCellInfo();
    void buildCellInfo2();
    bool m_bexon = false;
  public:
   float gef_area_;
  public:
    BgefReader(const string &filename, int bin_size, int n_thread = 1, bool verbose = false);
    virtual ~BgefReader();
    int getVersion() const;
    int getBinSize() const;
    unsigned int getGeneNum() const;
    unsigned int getCellNum();
    Gene *getGene();
    /**
     * @brief Get the number of expression.
     */
    unsigned int getExpressionNum() const;
    ExpressionAttr &getExpressionAttr();

    void getGeneExpression(unordered_map<string, vector<Expression>> & gene_exp_map, const vector<int>& regions);
    void getGeneExpression(unordered_map<string, vector<Expression>> & gene_exp_map);

    /**
     * @brief Get the shape of wholeExp matrix.
     * @return [rows, cols]
     */
    const unsigned int *getWholeExpMatrixShape();

    /**
     * @brief Get gene name list.
     *
     * @param gene_list
     */
    void getGeneNameList(vector<string> & gene_list);

    /**
     *  @brief Get cell name list.
     * @param cell_name_list
     */
    void getCellNameList(unsigned long long int * cell_name_list);

    unsigned long long int * getCellPos();

    Expression * getExpression();

    Expression * getReduceExpression();

    void cacheWholeExpMatrix();

    cv::Mat getWholeExpMatrix(cv::Rect roi);

    /**
     * @brief Read WholeExp data to matrix.
     * @param offset_x The starting position on the x-axis to be read.
     * @param offset_y The starting position on the y-axis to be read.
     * @param rows    Number of rows to read.
     * @param cols    Number of cols to read.
     * @param key     MIDcount or genecount.
     * @param matrix  When the value is greater than 255, it will be set to 255.
     */
    void readWholeExpMatrix(unsigned int offset_x,
                           unsigned int offset_y,
                           unsigned int rows,
                           unsigned int cols,
                           string & key,
                           unsigned char *matrix);

    /**
     * @brief Read WholeExp data to matrix.
     * @param key     MIDcount or genecount.
     * @param matrix  When the value is greater than 255, it will be set to 255.
     */
    void readWholeExpMatrix(string & key,
                            unsigned char *matrix);


    /**
     * Gets indices of cell for building csr_matrix.
     * @param cell_ind
     * @param count
     * @deprecated For special reasons, this function may be removed in future versions
     */
    //vector<unsigned long long int> getSparseMatrixIndicesOfExp(unsigned int * cell_ind, unsigned int * count);
    void getSparseMatrixIndicesOfExp(vector<unsigned long long> &uniq_cells, unsigned int * cell_ind, unsigned int * count);

    /**
     * @brief Gets indices of gene for building csr_matrix.
     * @param gene_ind
     * @param gene_names
     * @deprecated For special reasons, this function may be removed in future versions
     */
    void getSparseMatrixIndicesOfGene(unsigned int * gene_ind, char * gene_names);

    /*
     * @brief deprecated
     * @deprecated For special reasons, this function may be removed in future versions
     */
    vector<string> getSparseMatrixIndicesOfGene(unsigned int * gene_index);

    /**
     * @brief Gets indices for building csr_matrix.
     *
     * Examples:
     * @code
     * # Python
     * from scipy import sparse
     * sparse.csr_matrix((data, indices, indptr), shape=(cell_num, gene_num))
     * @endcode
     * @param indices  CSR format index array of the matrix. Cell id array, the column indices,
     * is the same size as count.
     * @param indptr   CSR format index pointer array of the matrix. indptr length = gene_num_ + 1 .
     * @param count    CSR format data array of the matrix. Expression count.
     * @return
     */
    int getSparseMatrixIndices(unsigned int * indices, unsigned int * indptr, unsigned int * count);

    /**
     * @brief Gets indices for building csr_matrix.
     *
     * @param cell_ind     CSR format index array of the matrix. same size as count.
     * @param gene_ind     CSR format index array of the matrix. same size as count.
     * @param count        CSR format data array of the matrix. Expression count.
     */
    int getSparseMatrixIndices2(unsigned int * cell_ind, unsigned int * gene_ind, unsigned int * count);

    /**
     * @brief Get geneID and expCount of this gene for each bin
     * @return key is bin id: x << 32 | y, value is a vector of gene_exp compoud value (geneID << 16 | geneExpCount).
     */
    void getBinGeneExpMap(map<unsigned long long int, pair<unsigned int, unsigned short>> &bin_exp_map,
                          DnbExpression * dnb_exp_info);

    void getGeneAndCount(unsigned short * gene_ind, unsigned short * count);

    //unsigned int toGem(string & filename, string &sn);

    /**
     * @brief Free memory for cache variables
     */
    void clear();

    static bool expressionComp(const DnbExpression& p1, const DnbExpression& p2);

    int generateGeneExp(int bin_size, int n_thread);

    void generateWholeExp(int size, int thread);

    void getfiltereddata(vector<int> &region, vector<string> &genelist,
                        vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                        vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, 
                        vector<unsigned int> &count);

    void getfiltereddata_exon(vector<int> &region, vector<string> &genelist,
                        vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                        vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, 
                        vector<unsigned int> &count, vector<unsigned int> &exon);

    void getOffset(int *data);
    void getExpAttr(int *data);

    bool isExonExist()
    {
        return m_exonPtr != nullptr;
    }
    unsigned int *getGeneExon();
    unsigned int getGeneExonAttr();

    Expression *getExpression_abs();
    bool isContainExon(){return m_bexon;}
    void openExonSpace(int bin_size);
};

#endif //GEFTOOLS__COMMON_BIN_H_
