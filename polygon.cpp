#include "polygon.h"
#include "utils.h"

using namespace GefTools;

bool Polygon::applyContour(const vector<cv::Point>& contour){
    original_contour_size_ = static_cast<short>(contour.size());

    if(contour.size() > BORDERCNT){
        double epsilon = 0.01 * cv::arcLength(contour, true);
        approxPolyDP(contour, border_, epsilon, true);
    }else{
        border_ = contour;
    }

    border_size_ = static_cast<short>(border_.size());
    if(border_size_ <= 2 ) return false;
    assert(border_size_ < 33);

    cv::Moments mu = cv::moments(border_, true);
    if(mu.m00 == 0) return false;

    center_ = cv::Point(static_cast<int>(mu.m10/mu.m00), static_cast<int>(mu.m01/mu.m00));

    area_ = mu.m00;

    for(const auto& p : border_){
        min_x_ = p.x < min_x_ ? p.x : min_x_;
        max_x_ = p.x > max_x_ ? p.x : max_x_;
        min_y_ = p.y < min_y_ ? p.y : min_y_;
        max_y_ = p.y > max_y_ ? p.y : max_y_;
    }

    for(const auto& p : border_){
        cv::Point relative_point = cv::Point(p.x - min_x_, p.y - min_y_);
        relative_border_.emplace_back(relative_point);
    }

    cols_ = max_x_ - min_x_ + 1;
    rows_ = max_y_ - min_y_ + 1;

    return true;
}

const vector<cv::Point> &Polygon::getBorder() const {
    return border_;
}

const cv::Point &Polygon::getCenter() const {
    return center_;
}

double Polygon::getArea() const {
    return area_;
}

short Polygon::getBorderSize() const {
    return border_size_;
}

short Polygon::getOriginalContourSize() const {
    return original_contour_size_;
}

void Polygon::setMinMaxXY() {
    for(const auto& p : border_){
        min_x_ = p.x < min_x_ ? p.x : min_x_;
        max_x_ = p.x > max_x_ ? p.x : max_x_;
        min_y_ = p.y < min_y_ ? p.y : min_y_;
        max_y_ = p.y > max_y_ ? p.y : max_y_;
    }
    rows_ = max_x_ - min_x_ + 1;
    cols_ = max_y_ - min_y_ + 1;
}

int Polygon::getMinX() const {
    return min_x_;
}

int Polygon::getMaxX() const {
    return max_x_;
}

int Polygon::getMinY() const {
    return min_y_;
}

int Polygon::getMaxY() const {
    return max_y_;
}

int Polygon::getRows() const {
    return rows_;
}

int Polygon::getCols() const {
    return cols_;
}

cv::Mat Polygon::getFillPolyMat() const {
    cv::Mat fill_points = cv::Mat::zeros(rows_, cols_, CV_8UC1);
    fillPoly(fill_points, relative_border_, 1);
    return fill_points;
}

const vector<cv::Point> &Polygon::getRelativeBorder() const {
    return relative_border_;
}

// ceil up to the nearest (unsigned short) integer
unsigned short Polygon::getAreaUshort() const {
    return ceil(area_);
}

unsigned int Polygon::getBlockId() const {
    return block_id_;
}

void Polygon::setBlockId(const unsigned int* block_size) {
    unsigned int block_id = center_.x / block_size[0];
    unsigned int block_id_y = center_.y / block_size[1];
    block_id_ = block_id + block_id_y * block_size[2];
}