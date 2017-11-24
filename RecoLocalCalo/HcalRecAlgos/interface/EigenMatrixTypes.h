#include <Eigen/Dense>

constexpr int MaxSVSize=10;
constexpr int MaxFSVSize=19;
constexpr int MaxPVSize=10;

typedef Eigen::Matrix<double,MaxSVSize,1> SampleVector;
typedef Eigen::Matrix<double,MaxFSVSize,1> FullSampleVector;
typedef Eigen::Matrix<double,Eigen::Dynamic,1,0,MaxPVSize,1> PulseVector;
typedef Eigen::Matrix<int,Eigen::Dynamic,1,0,MaxPVSize,1> BXVector;
typedef Eigen::Matrix<double,MaxSVSize,MaxSVSize> SampleMatrix;
typedef Eigen::Matrix<double,MaxFSVSize,MaxFSVSize> FullSampleMatrix;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,MaxPVSize,MaxPVSize> PulseMatrix;
typedef Eigen::Matrix<double,10,Eigen::Dynamic,0,MaxSVSize,MaxPVSize> SamplePulseMatrix;
typedef Eigen::LLT<SampleMatrix> SampleDecompLLT;
typedef Eigen::LLT<PulseMatrix> PulseDecompLLT;
typedef Eigen::LDLT<PulseMatrix> PulseDecompLDLT;

typedef Eigen::Matrix<double,1,1> SingleMatrix;
typedef Eigen::Matrix<double,1,1> SingleVector;

//typedef Eigen::Matrix<double,10,1> SampleVector;
//typedef Eigen::Matrix<double,12,1> FullSampleVector;
//typedef Eigen::Matrix<double,Eigen::Dynamic,1,0,10,1> PulseVector;
//typedef Eigen::Matrix<char,Eigen::Dynamic,1,0,10,1> BXVector;
//typedef Eigen::Matrix<double,10,10> SampleMatrix;
//typedef Eigen::Matrix<double,12,12> FullSampleMatrix;
//typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,10,10> PulseMatrix;
//typedef Eigen::Matrix<double,10,Eigen::Dynamic,0,10,10> SamplePulseMatrix;
//typedef Eigen::LLT<SampleMatrix> SampleDecompLLT;
//typedef Eigen::LLT<PulseMatrix> PulseDecompLLT;
//typedef Eigen::LDLT<PulseMatrix> PulseDecompLDLT;
//
//typedef Eigen::Matrix<double,1,1> SingleMatrix;
//typedef Eigen::Matrix<double,1,1> SingleVector;
