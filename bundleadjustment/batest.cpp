#include<iostream>
#include<opencv2/opencv.hpp>

 
 
using namespace std;
using namespace cv;
 
 
int main(int argc, char** argv) {
 
 
  Mat src = imread("/home/zzhfro/code/ImageBasedModellingEduV1.0/examples/data/kxm1.jpg");
 
  imshow("src", src);
 
 
  Mat src2 = imread("/home/zzhfro/code/ImageBasedModellingEduV1.0/examples/data/kxm2.jpg");

  imshow("src2", src2);
 
 
  //定义Sift的基本参数
  int numFeatures = 500;
  //创建detector存放到KeyPoints中
  Ptr<SIFT> detector = SIFT::create(numFeatures);
  vector<KeyPoint> keypoints, keypoints2;
  detector->detect(src, keypoints);
  detector->detect(src2, keypoints2);
  //打印Keypoints
  cout << "Keypoints:" << keypoints.size() << endl;
  cout << "Keypoints2:" << keypoints2.size() << endl;
 
 
  Mat drawsrc, drawsrc2;
  drawKeypoints(src, keypoints, drawsrc);

  imshow("drawsrc", drawsrc);
  drawKeypoints(src2, keypoints2, drawsrc2);
 
  imshow("drawsrc2", drawsrc2);
 
 
  //计算特征点描述符,特征向量提取
  Mat dstSIFT, dstSIFT2;
  Ptr<SiftDescriptorExtractor> descriptor = SiftDescriptorExtractor::create();
  descriptor->compute(src, keypoints, dstSIFT);
  descriptor->compute(src2, keypoints2, dstSIFT2);
  cout << dstSIFT.cols << endl;
  cout << dstSIFT2.rows << endl;
 
 
  //进行BFMatch暴力匹配
  BFMatcher matcher(NORM_L2);
  //定义匹配结果变量
  vector<DMatch> matches;
  //实现描述符之间的匹配
  matcher.match(dstSIFT, dstSIFT2, matches);
 
 
  //定义向量距离的最大值与最小值
  double max_dist = 0;
  double min_dist = 1000;
  for (int i = 1; i < dstSIFT.rows; ++i)
  {
    //通过循环更新距离，距离越小越匹配
    double dist = matches[i].distance;
    if (dist > max_dist)
      max_dist = dist;
    if (dist < min_dist)
      min_dist = dist;
  }
  cout << "min_dist=" << min_dist << endl;
  cout << "max_dist=" << max_dist << endl;
  //匹配结果筛选    
  vector<DMatch> goodMatches;
  for (int i = 0; i < matches.size(); ++i)
  {
    double dist = matches[i].distance;
    if (dist < 2 * min_dist)
      goodMatches.push_back(matches[i]);
  }
  cout << "goodMatches:" << goodMatches.size() << endl;
 
  for (const auto& match : goodMatches) {
        std::cout << "QueryIdx: " << match.queryIdx << ", ";
        std::cout << "TrainIdx: " << match.trainIdx << ", ";
        std::cout << "Distance: " << match.distance << std::endl;
        cv::Point2f pt1 = keypoints[match.queryIdx].pt; // 查询图像中的特征点位置
        cv::Point2f pt2 = keypoints2[match.trainIdx].pt; // 训练图像中的特征点位置

        std::cout << "Matched points: (" << pt1.x << ", " << pt1.y << ") - (" 
                  << pt2.x << ", " << pt2.y << ")" << std::endl;
    }
  Mat result;
  //匹配特征点天蓝色，单一特征点颜色随机
  drawMatches(src, keypoints, src2, keypoints2, goodMatches, result, 
    Scalar(255, 255, 0), Scalar::all(-1));
  
  imshow("Result", result);
 
 
 
 
  waitKey(0);
  return 0;
}