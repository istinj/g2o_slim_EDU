#include "utilities.h"
using namespace Eigen;
using namespace std;

Matrix4f v2t(const Vector6f& v){
   Matrix4f T = Matrix4f::Identity();
   Matrix3f Rx, Ry, Rz;
   Rx = AngleAxisf(v(3), Vector3f::UnitX());
   Ry = AngleAxisf(v(4), Vector3f::UnitY());
   Rz = AngleAxisf(v(5), Vector3f::UnitZ());
   T.block<3,3>(0,0) = Rx * Ry * Rz;
   T.block<3,1>(0,3) = v.block<3,1>(0,0);
   return T;
}

Eigen::Matrix3f skew(const Eigen::Vector3f& p)
{
   Eigen::Matrix3f s;
   s <<
         0,  -p.z(), p.y(),
         p.z(), 0,  -p.x(),
         -p.y(), p.x(), 0;
   return s;
}

void cholesky(const Matrix6f& input_, Matrix6f& L_){
   int dim = input_.rows();
   L_.setZero();
   for (int i = 0; i < dim; i++)
      for (int j = 0; j < (i+1); j++) {
         float s = 0;
         for (int k = 0; k < j; k++)
            s += L_(i, k) * L_(j, k);
         L_(i, j) = (i == j) ?
               sqrt((float)input_(i, j) - s) :
               (1.0 / L_(j,j) * (input_(i,j) - s));
      }
   L_.transposeInPlace();
}

void loadMatrix(const std::string& name_, Eigen::MatrixXf& data_){
   cout << BOLDYELLOW << "\t" << "Opening file " << name_ << RESET << endl;
   fstream file(name_);
   string line;
   int size = 0;

   getline(file, line);
   stringstream ss(line);
   float num;
   while(ss >> num)
      size++;


   data_.resize(size, size);
   int i = 0;
   while(getline(file, line)){
      stringstream ss(line);
      for(int j = 0; j < size; ++j){
         ss >> data_(i,j);
      }
      ++i;
   }
   cout << BOLDGREEN << "\t" << "Matrix loaded successfully"  << RESET << endl;
   return;
}
