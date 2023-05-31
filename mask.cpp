//
// Created by huangzhibo on 2021/12/14.
//

#include "mask.h"
Mask::Mask(const string& file, const int block_size[], const unsigned int mask_size[]){
    cv::Mat img = cv::imread(file,-1);
    if( img.empty() ) 
    { 
        cerr << "Mask is empty!" << endl;
        reportErrorCode2File(errorCode::E_FILEOPENERROR, 
                            "Mask is empty!");
        exit(-1);
    }

    if(mask_size[0]> 0 && mask_size[1] > 0)
    {
        if (img.rows != mask_size[0] || img.cols != mask_size[1]) {
            if (img.rows == mask_size[1] && img.cols == mask_size[0]) {
                img = img.t();
            } else {
                cerr << "The size of mask picture is inconsistent with the size of expression" << endl;
                reportErrorCode2File(errorCode::E_FILEMISMATCH, 
                            "The size of mask picture is inconsistent with the size of expression");
                exit(2);
            }
        }else if(img.rows == img.cols){
            cerr << "[WARN] Mask rows == cols, the mask coordinates are not automatically adjusted (transposed)"
                    "to be consistent with the expression coordinates." << endl;
        }
    }

    rows_ = img.rows;
    cols_ = img.cols;
    // x_block_size, y_block_size, x_block_num, y_block_num
    block_size_[0] = block_size[0];
    block_size_[1] = block_size[1];
    block_size_[2] = ceil(cols_ * 1.0 / block_size[0]);
    block_size_[3] = ceil(rows_ * 1.0 / block_size[1]);

    //findContours从二值图像中检索轮廓，并返回检测到的轮廓的个数
    cv::findContours(
            img,
            contours_,
            hierarchy_,
            cv::RETR_EXTERNAL,
            cv::CHAIN_APPROX_SIMPLE
    );
    block_num_ = block_size_[2] * block_size_[3];

    for(auto & contour : contours_){
        GefTools::Polygon p;
        bool cell_polygon_is_good = p.applyContour(contour);
        if(cell_polygon_is_good){
            p.setBlockId(block_size_);
            polygons_.emplace_back(std::move(p));
            min_x_ = p.getMinX() < min_x_ ? p.getMinX() : min_x_;
            max_x_ = p.getMaxX() > max_x_ ? p.getMaxX() : max_x_;
            min_y_ = p.getMinY() < min_y_ ? p.getMinY() : min_y_;
            max_y_ = p.getMaxY() > max_y_ ? p.getMaxY() : max_y_;
        }
    }

    preBlockSort();

    cell_num_ = polygons_.size();
}

Mask::~Mask() {
    if(block_index_ != nullptr)
        free(block_index_);
}

const vector<GefTools::Polygon> &Mask::getPolygons() const {
    return polygons_;
}

unsigned int Mask::getCellNum() const {
    return cell_num_;
}

void Mask::showMaskInWindow() {
    cv::Mat cnt_img = cv::Mat::zeros(rows_, cols_, CV_8UC3);
    drawContours( cnt_img, contours_, -1, cv::Scalar(128,255,255),
                  3, cv::LINE_AA, hierarchy_, 3 );
    imshow("Mask Contours", cnt_img);
    cv::waitKey(0);
}

void Mask::getBorders(short * border_array) {
    for (unsigned int i = 0; i < cell_num_; i++) {
        GefTools::Polygon polygon = polygons_[i];
        vector<cv::Point> border = polygon.getBorder();
        const cv::Point& center = polygon.getCenter();
        unsigned int index1 = i * 64;

        auto border_size = static_cast<short>(border.size());

        for (short j = 0; j < BORDERCNT; ++j) {
            unsigned int index2 = index1 + (j << 1);
            if(j >= border_size){
                border_array[index2] = SHRT_MAX;
                border_array[index2+1] = SHRT_MAX;
            }else{
                cv::Point p = border[j];
                border_array[index2] = static_cast<short>(p.x - center.x);
                border_array[index2+1] = static_cast<short>(p.y - center.y);
            }
        }
    }
}

/**
 * @brief Gets the coordinates of effective rectangle region.
 * @param effective_rect  Should be a int array with 4 size. They represent min_x, min_y, max_x, max_y in turn.
 */
void Mask::getEffectiveRectangle(int* effective_rect) const {
    effective_rect[0] = min_x_;
    effective_rect[1] = min_y_;
    effective_rect[2] = max_x_;
    effective_rect[3] = max_y_;
}

void Mask::preBlockSort() {
    sort(polygons_.begin(), polygons_.end(), polygonComp);
}

bool Mask::polygonComp(const GefTools::Polygon& p1, const GefTools::Polygon& p2) {
    return p1.getBlockId() < p2.getBlockId();
}

unsigned int Mask::getBlockNum() const {
    return block_num_;
}

unsigned int * Mask::getBlockIndex() {
    if(block_index_ != nullptr)
        return block_index_;

    unsigned int block_index_size = block_num_ + 1;

    block_index_ = (unsigned int *) calloc(block_index_size, sizeof(unsigned int));

    for(unsigned int i = 0; i< cell_num_; i++){
        GefTools::Polygon p = polygons_[i];
        block_index_[p.getBlockId()] += 1;
    }

    block_index_[block_num_] = cell_num_;
    for(unsigned int i = block_num_; i > 0; i--){
        block_index_[i-1] = block_index_[i] - block_index_[i-1];
    }
    return block_index_;
}

const unsigned int *Mask::getBlockSize() const {
    return block_size_;
}

