#include "RecoLocalCalo/HcalRecAlgos/interface/DoMahiAlgo.h"
#include <iostream>
#include <fstream> 

void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP);

void DoMahiAlgo::setPulseShapeTemplate(bool useCSV, std::string filename="") {
  _useCSV = useCSV;

  if (_useCSV && filename!="") {
    std::ifstream ifs;
    ifs.open(filename.c_str());
    assert(ifs.is_open());
    std::string line;

    int i = 0;
    while(getline(ifs,line)) {
      if(line[0]=='#') continue;
      
      std::string tmpStr;
      std::stringstream ss(line);
      ss >> tmpStr; 
      minCharge_[i] = std::atoi(tmpStr.c_str());
      ss >> tmpStr;
      maxCharge_[i] = std::atoi(tmpStr.c_str());
      for (int k=0; k<10; k++) { ss >> tmpStr; pulseFrac_[i][k] = std::atof(tmpStr.c_str()); }
      for (int k=0; k<10; k++) { ss >> tmpStr; pulseFracDeriv_[i][k] = std::atof(tmpStr.c_str()); }

      i++;

    }
  }

}

void DoMahiAlgo::getPulseShape(float q, HcalDetId detID, float t, SampleVector &pulseShape, float sigma=0) {

  if (!_useCSV) {
    pulseShapeObj.computeLAGShape(q, detID, t, pulseShape, sigma);
  }
  else {
    int chargeBin = -1;
    for (int i=0; i<58; i++) {
      if (q>minCharge_[i] && q<maxCharge_[i]) chargeBin=i;
    }
    if (q>maxCharge_[57]) chargeBin=57; 
    if (chargeBin==-1) chargeBin=0;
    
    int dt= (t/25);
    for (int i=0; i<10; i++) {
      if (i-dt<0 ||i-dt>9) pulseShape.coeffRef(i-dt) = 0;
      else {
	if (pulseFrac_[chargeBin][i] > 1e-6) {
	  pulseShape.coeffRef(i-dt) = pulseFrac_[chargeBin][i-dt] + (t-dt*25)*pulseFracDeriv_[chargeBin][i-dt];
	}
	else {
	  pulseShape.coeffRef(i-dt)=0;
	}
      }
    }
  }

}


void DoMahiAlgo::Apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalDetId & detID, const HcalCalibrations & calibs, std::vector<double> & correctedOutput) {

  const unsigned int cssize = cs.size();

  SampleVector charges;
  SampleVector gains;

  double tsTOT = 0, tstrig = 0; // in fC
  for(unsigned int ip=0; ip<cssize; ++ip){
    if( ip >= (unsigned)10) continue; // Too many samples than what we wanna fit (10 is enough...) -> skip them
    const int capid = capidvec[ip];
    double charge = cs[ip];
    double ped = calibs.pedestal(capid);
    double gain = calibs.respcorrgain(capid);

    charges.coeffRef(ip) = charge - ped;
    gains.coeffRef(ip) = gain;

    tsTOT += charge - ped;
    if( ip ==4 || ip==5 ){
      tstrig += charge - ped;
    }
  }

  std::vector<double> fitParsVec;

  _detID = detID;

  bool status =false;
  if(tstrig >= 0) {
    status = DoFit(charges, gains, fitParsVec);
  }

  if (!status) {
    fitParsVec.clear();
    fitParsVec.push_back(0.);
    fitParsVec.push_back(0.);
    fitParsVec.push_back(0.);
    fitParsVec.push_back(999.);
  }

  correctedOutput.swap(fitParsVec);

}  

bool DoMahiAlgo::DoFit(SampleVector amplitudes, SampleVector gains, std::vector<double> &correctedOutput) {

  _nP = 0;
  //ECAL does it better -- to be fixed
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#L151-L171
  _bxs.resize(3);
  _bxs << -1,0,1;
  _nPulseTot = _bxs.rows();

  //_nPulseTot = _bxs.rows();
  //_bxs = bxs;
  //_bxs.resize(3);
  //_bxs << -1,0,1;
  //_detID = HcalDetId(detID.rawId());
    
  _amplitudes = amplitudes;

  _pulseMat.resize(Eigen::NoChange,_nPulseTot);
  _ampVec = PulseVector::Zero(_nPulseTot);
  _errVec = PulseVector::Zero(_nPulseTot);

  _ampVec.coeffRef(0) = _amplitudes.coeff(3);
  _ampVec.coeffRef(1) = _amplitudes.coeff(4);
  _ampVec.coeffRef(2) = _amplitudes.coeff(5);

  _chiSq = 9999;

  aTaMat.resize(_nPulseTot, _nPulseTot);
  aTbVec.resize(_nPulseTot);
  wVec.resize(_nPulseTot);

  pulseShape = PulseVector::Zero(10);
  //pulseShapeObj.computeLAGShape(_amplitudes.coeff(3), 0, -25, pulseShape,0);
  //pulseShapeObj.computeLAGShape(_amplitudes.coeff(4), 0, 0, pulseShape,0);
  //pulseShapeObj.computeLAGShape(_amplitudes.coeff(5), 0, 25, pulseShape,0);
  _pulseMat.col(0) = pulseShape.segment<10>(0);
  _pulseMat.col(1) = pulseShape.segment<10>(0);
  _pulseMat.col(2) = pulseShape.segment<10>(0);

  //std::cout << "initial pulseMat" << std::endl;
  //std::cout << _pulseMat << std::endl;

  bool status = Minimize(); 
  _ampVecMin = _ampVec;
  _bxsMin = _bxs;

  if (!status) return status;

  bool foundintime = false;
  unsigned int ipulseintime = 0;
  unsigned int ipulseprevtime = 0;
  unsigned int ipulsenexttime = 0;

  for (unsigned int ipulse=0; ipulse<_nPulseTot; ++ipulse) {
    if (_bxs.coeff(ipulse)==0) {
      ipulseintime = ipulse;
      foundintime = true;
    }
    else if (_bxs.coeff(ipulse)==-1) {
      ipulseprevtime = ipulse;
    }
    else if (_bxs.coeff(ipulse)==1) {
      ipulsenexttime = ipulse;
    }
  }
  if (!foundintime) return status;

  //std::cout << "------" << std::endl;
  //std::cout << "input: " ;
  //for (int i=0; i<10; i++) {
  //std::cout << _amplitudes.coeff(i) << ", ";
  //}
  //std::cout << std::endl;

  //std::vector<double> ans;

  //std::cout << "output: ";// << std::endl;
  //std::cout << _ampVec.coeff(ipulseprevtime) << ", " << _ampVec.coeff(ipulseintime) << ", " << _ampVec.coeff(ipulsenexttime) << std::endl;
  /*std::cout << "output we care about: ";

  pulseShapeObj.computeLAGShape(_ampVec.coeff(ipulseprevtime), _detID, -25, pulseShape,0);
  for (int i=0; i<10; i++) {
    ans.push_back(_ampVec.coeff(ipulseprevtime)*pulseShape.coeff(i));
    std::cout << _ampVec.coeff(ipulseprevtime)*pulseShape.coeff(i) << ", ";
  }
  std::cout << std::endl;

  pulseShapeObj.computeLAGShape(_ampVec.coeff(ipulseintime), _detID, 0, pulseShape,0);
  for (int i=0; i<10; i++) {
    ans[i]+=_ampVec.coeff(ipulseintime)*pulseShape.coeff(i);
    std::cout << _ampVec.coeff(ipulseintime)*pulseShape.coeff(i) << ", ";
  }
  std::cout << std::endl;

  pulseShapeObj.computeLAGShape(_ampVec.coeff(ipulsenexttime), _detID, 25, pulseShape,0);
  for (int i=0; i<10; i++) {
    ans[i]+=_ampVec.coeff(ipulsenexttime)*pulseShape.coeff(i);
    std::cout << _ampVec.coeff(ipulsenexttime)*pulseShape.coeff(i) << ", ";
  }
  std::cout << std::endl;
  for (int i=0; i<10; i++) {
    std::cout << ans[i] << ", ";
  }
  std::cout << std::endl;
  */
  //std::cout << "chi2: " ;
  //std::cout << _chiSq << std::endl;
  //std::cout << "-----------" << std::endl;

  //const unsigned int ipulseintimemin = ipulseintime;
  /*
  return _ampVec.coeff(ipulseintime);

  //return true;
  */

  double gain=gains.coeff(4); // this is the same for each TS

  correctedOutput.clear();
  //  correctedOutput.push_back(_ampVec.coeff(ipulseintime)); //charge
  correctedOutput.push_back(_ampVec.coeff(ipulseintime)*gain); //energy
  correctedOutput.push_back(_ampVec.coeff(ipulsenexttime)*gain); //energy TEMPORARY
  correctedOutput.push_back(_ampVec.coeff(ipulseprevtime)*gain); //energy TEMPORARY
  //correctedOutput.push_back(-999); //time
  ///correctedOutput.push_back(-999); //pedestal
  correctedOutput.push_back(_chiSq); //chi2

  return true;

}

bool DoMahiAlgo::Minimize() {
  //std::cout << "start minimize" << std::endl;
  int iter = 0;
  int maxIters = 500;
  bool status = false;

  while (true) {
    if (iter>=maxIters) {
      std::cout << "max number of iterations reached! " << std::endl;
      //std::cout << _chiSq << std::endl;
      break;
    }
    //std::cout << "iteration "<< iter << std::endl;

    status = UpdateCov();
    if (!status) break;

    status = NNLS();
    if (!status) break;

    double newChiSq = CalculateChiSq();
    double deltaChiSq = newChiSq - _chiSq;
    
    _chiSq = newChiSq;

    //std::cout << "chiSq = " << _chiSq << ", " << deltaChiSq << std::endl;
    
    if (std::abs(deltaChiSq)<1e-3) break;

    iter++;

  }
  
  return true;

}

bool DoMahiAlgo::NNLS() {
  //std::cout << "start NNLS" << std::endl;
  const unsigned int npulse = _bxs.rows();

  //std::cout << _ampVec << std::endl;
  for (uint i=0; i<npulse; i++) {
    double nomT = _bxs.coeff(i)*25;
    getPulseShape(_ampVec.coeff(i), _detID, nomT, pulseShape,0);
    _pulseMat.col(i) = pulseShape.segment<10>(0);
  }

  //std::cout << "new pulsemat" << std::endl;
  //std::cout << _pulseMat << std::endl;
  
  invcovp = _covDecomp.matrixL().solve(_pulseMat);
  aTaMat = invcovp.transpose()*invcovp;
  aTbVec = invcovp.transpose()*_covDecomp.matrixL().solve(_amplitudes);

  int iter = 0;
  Index idxwmax = 0;
  double wmax = 0.0;
  double threshold = 1e-11;
  
  while (true) {    
    if (iter>0 || _nP==0) {
      if ( _nP==npulse ) break;

      const unsigned int nActive = npulse - _nP;
      updateWork = aTbVec - aTaMat*_ampVec;

      Index idxwmaxprev = idxwmax;
      double wmaxprev = wmax;
      wmax = updateWork.tail(nActive).maxCoeff(&idxwmax);

      if (wmax<threshold || (idxwmax==idxwmaxprev && wmax==wmaxprev)) {
	break;
      }
      
      if (iter>=500) {
	std::cout << "Max Iterations reached!" << std::endl;
	break;
      }

      //unconstrain parameter
      Index idxp = _nP + idxwmax;
      
      aTaMat.col(_nP).swap(aTaMat.col(idxp));
      aTaMat.row(_nP).swap(aTaMat.row(idxp));
      _pulseMat.col(_nP).swap(_pulseMat.col(idxp));
      std::swap(aTbVec.coeffRef(_nP),aTbVec.coeffRef(idxp));
      std::swap(_ampVec.coeffRef(_nP),_ampVec.coeffRef(idxp));
      std::swap(_bxs.coeffRef(_nP),_bxs.coeffRef(idxp));

      wVec.tail(nActive) = updateWork.tail(nActive); 
      
      ++_nP;

    }

    while (true) {
      if (_nP==0) break;     

      ampvecpermtest = _ampVec;

      eigen_solve_submatrix(aTaMat,aTbVec,ampvecpermtest,_nP);

      //check solution    
      //std::cout << "nP..... " << _nP << std::endl;
      auto ampvecpermhead = ampvecpermtest.head(_nP);

      if ( ampvecpermhead.minCoeff()>0. ) {
	_ampVec.head(_nP) = ampvecpermhead.head(_nP);
	//std::cout << "eep?" << std::endl;
	break;
      }

      //update parameter vector
      Index minratioidx=0;

      // no realizable optimization here (because it autovectorizes!)
      double minratio = std::numeric_limits<double>::max();
      for (unsigned int ipulse=0; ipulse<_nP; ++ipulse) {
	if (ampvecpermtest.coeff(ipulse)<=0.) {
	  const double c_ampvec = _ampVec.coeff(ipulse);
	  const double ratio = c_ampvec/(c_ampvec-ampvecpermtest.coeff(ipulse));
	  if (ratio<minratio) {
	    minratio = ratio;
	    minratioidx = ipulse;
	  }
	}
      }
      _ampVec.head(_nP) += minratio*(ampvecpermhead - _ampVec.head(_nP));
      
      //avoid numerical problems with later ==0. check
      _ampVec.coeffRef(minratioidx) = 0.;
      
      aTaMat.col(_nP-1).swap(aTaMat.col(minratioidx));
      aTaMat.row(_nP-1).swap(aTaMat.row(minratioidx));
      _pulseMat.col(_nP-1).swap(_pulseMat.col(minratioidx));
      std::swap(aTbVec.coeffRef(_nP-1),aTbVec.coeffRef(minratioidx));
      std::swap(_ampVec.coeffRef(_nP-1),_ampVec.coeffRef(minratioidx));
      std::swap(_bxs.coeffRef(_nP-1),_bxs.coeffRef(minratioidx));
      --_nP;

    }
   
    ++iter;
    
    //adaptive convergence threshold to avoid infinite loops but still
    //ensure best value is used
    if (iter%10==0) {
      threshold *= 10.;
    }

    break;
  }

  return true;

}

bool DoMahiAlgo::UpdateCov() {
  //std::cout << "start update Cov" << std::endl;
  const double pederr2 = 1;
  _invCovMat = pederr2*SampleMatrix::Constant(1);

  for (int i=0; i<10; i++) {
    //std::cout << _amplitudes.coeff(i) << std::endl;
    double ifC=_amplitudes.coeff(i);
    double sigma = 0;
    if(ifC < 75) sigma = (0.577 + 0.0686*ifC)/3.; 
    else sigma = (2.75  + 0.0373*ifC + 3e-6*ifC*ifC)/3.; 
    _invCovMat(i, i) += (1+sigma)*(1+sigma);
  }

  for (int k=0; k<_ampVec.size(); k++) {
    double ifC=_ampVec.coeff(k);
    if (ifC==0) continue;
    
    double nomT = _bxs.coeff(k)*25;
    int maxTS = 4+_bxs.coeff(k);

    getPulseShape(ifC, _detID, nomT, pulseShape,0); 
    getPulseShape(ifC, _detID, nomT+deltaT, pulseShapeP,0);  
    getPulseShape(ifC, _detID, nomT-deltaT, pulseShapeM,0);

    for (int xx=0; xx<10; xx++) {
      pulseShape.coeffRef(xx) = pulseShape.coeff(xx)/pulseShape.coeff(maxTS);
      pulseShapeP.coeffRef(xx) = pulseShapeP.coeff(xx)/pulseShapeP.coeff(maxTS);
      pulseShapeM.coeffRef(xx) = pulseShapeM.coeff(xx)/pulseShapeM.coeff(maxTS);
    }

    //if (nomT!=0) std::cout << _invCovMat << std::endl << "---" << std::endl;
    for (int i=0; i<10; i++) {
      for (int j=0; j<i+1; j++) {
	double tmp=0.5*((pulseShapeP.coeff(i)-pulseShape.coeff(i))*(pulseShapeP.coeff(j)-pulseShape.coeff(j)) 
			+ (pulseShapeM.coeff(i)-pulseShape.coeff(i))*(pulseShapeM.coeff(j)-pulseShape.coeff(j)));
	_invCovMat(i,j) += ifC*ifC*tmp;
	_invCovMat(j,i) += ifC*ifC*tmp;
	
      }
    }
    //if (nomT!=0)std::cout << _invCovMat << std::endl;

  }
  //std::cout << "cov" << std::endl;
  //std::cout << _invCovMat << std::endl;
  //std::cout << "..." << std::endl;
  _covDecomp.compute(_invCovMat);

  return true;
  
}

double DoMahiAlgo::CalculateChiSq() {
  return _covDecomp.matrixL().solve(_pulseMat*_ampVec - _amplitudes).squaredNorm();
}

void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP) {
  //std::cout << "start solving " << NP << std::endl;
  //std::cout << mat << std::endl;
  using namespace Eigen;
  switch( NP ) { // pulse matrix is always square.
  case 10:
    {
      Matrix<double,10,10> temp = mat;
      outvec.head<10>() = temp.ldlt().solve(invec.head<10>());
    }
    break;
  case 9:
    {
      Matrix<double,9,9> temp = mat.topLeftCorner<9,9>();
      outvec.head<9>() = temp.ldlt().solve(invec.head<9>());
    }
    break;
  case 8:
    {
      Matrix<double,8,8> temp = mat.topLeftCorner<8,8>();
      outvec.head<8>() = temp.ldlt().solve(invec.head<8>());
    }
    break;
  case 7:
    {
      Matrix<double,7,7> temp = mat.topLeftCorner<7,7>();
      outvec.head<7>() = temp.ldlt().solve(invec.head<7>());
    }
    break;
  case 6:
    {
      Matrix<double,6,6> temp = mat.topLeftCorner<6,6>();
      outvec.head<6>() = temp.ldlt().solve(invec.head<6>());
    }
    break;
  case 5:
    {
      Matrix<double,5,5> temp = mat.topLeftCorner<5,5>();
      outvec.head<5>() = temp.ldlt().solve(invec.head<5>());
    }
    break;
  case 4:
    {
      Matrix<double,4,4> temp = mat.topLeftCorner<4,4>();
      outvec.head<4>() = temp.ldlt().solve(invec.head<4>());
    }
    break;
  case 3: 
  {
  Matrix<double,3,3> temp = mat.topLeftCorner<3,3>();
  outvec.head<3>() = temp.ldlt().solve(invec.head<3>());
  }
    break;
  case 2:
    {
      Matrix<double,2,2> temp = mat.topLeftCorner<2,2>();
      outvec.head<2>() = temp.ldlt().solve(invec.head<2>());
    }
    break;
  case 1:
    {
      Matrix<double,1,1> temp = mat.topLeftCorner<1,1>();
      outvec.head<1>() = temp.ldlt().solve(invec.head<1>());
    }
    break;
  default:
    //throw cms::Exception("MultFitWeirdState")
    std::cout << "Weird number of pulses encountered in multifit, module is configured incorrectly!" << std::endl;
    }
}

DoMahiAlgo::DoMahiAlgo() {
  pulseShapeObj = PulseShapes();
  deltaT = 2.5;
}
