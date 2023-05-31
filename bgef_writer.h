/** @file bgef_writer.h
    @brief Declare a BgefWriter class for writing common bin gef.

    Created by huangzhibo on 2021/01/24.
*/

#ifndef GEFTOOLS_BGEF_WRITER_H
#define GEFTOOLS_BGEF_WRITER_H

#include "hdf5.h"
#include "utils.h"
#include "gef.h"

static constexpr unsigned int version = 2;

class GEFTOOLS_API BgefWriter {
  private:
    int binsize_;
    hid_t str32_type_;
    hid_t str64_type_;
    hid_t file_id_;
    hid_t gene_exp_group_id_;
    hid_t whole_exp_group_id_;
    hid_t m_wholeExpExon_id;

    unsigned int resolution_;
    bool verbose_ = false;
    bool m_bexon = false;
    bool raw_gef_ = false;
  public:
    BgefWriter(const string& output_filename, bool verbose, bool bexon, const string& stromics);
    BgefWriter(const string& output_filename, unsigned int raw_gef_version);
    ~BgefWriter();

    bool storeGene(vector<Expression> &exps, vector<Gene> &genes, DnbAttr &dnbAttr, unsigned int maxexp, int binsize);
    bool storeDnb(DnbMatrix & dnb_matrix, int binsize);
    bool storeStat(vector<GeneStat>& geneStat) const;

    unsigned int getResolution() const;

    void setResolution(unsigned int resolution);
    bool storeGeneExon(vector<Expression>& exps, unsigned int maxexon, int binsize);
    bool storeWholeExon(DnbMatrix & dnb_matrix, int binsize);

    void StoreRawGef(Expression *exps, unsigned int exp_size, ExpressionAttr &exp_attr, Gene *genes, 
                     unsigned int gene_cnt, unsigned int *exons, unsigned int maxexon);
                     
    void SetGefArea(float &area);
};

#endif