#include "RecoLocalCalo/HcalRecAlgos/interface/DoMahiAlgo.h" 
#include <iostream>
#include <fstream> 

void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP);

DoMahiAlgo::DoMahiAlgo() { 
  std::cout << "hai" << std::endl;
}

void DoMahiAlgo::phase1Apply(const HBHEChannelInfo& channelData) {
  //float& reconstructedEnergy,
  //float& chi2) {
  
  std::cout << "phase1Apply" << std::endl;
  
  const unsigned cssize = channelData.nSamples();
  std::cout << cssize << std::endl;

  std::vector<float> reconstructedVals;

  SampleVector charges;
  SampleVector gains;

  double tsTOT = 0, tstrig = 0; // in fC
  for(unsigned int ip=0; ip<cssize; ++ip){
    if( ip >= (unsigned)10) continue; // Too many samples than what we wanna fit (10 is enough...) -> skip them
    //    const int capid = capidvec[ip];
    double charge = channelData.tsRawCharge(ip);
    double ped = channelData.tsPedestal(ip);
    double gain = channelData.tsGain(ip);

    charges.coeffRef(ip) = charge - ped;
    gains.coeffRef(ip) = gain;

    tsTOT += charge - ped;
    if( ip ==4 || ip==5 ){
      tstrig += charge - ped;
    }
  }


  _detID = channelData.id();

  //bool status =false;
  if(tstrig >= 0) {
    //status = 
    DoFit(charges, gains, reconstructedVals); 
  }
  
}

bool DoMahiAlgo::DoFit(SampleVector amplitudes, SampleVector gains, std::vector<float> &correctedOutput) {
  std::cout << "DoFit" << std::endl;

  _nP = 0;
  //ECAL does it better -- to be fixed
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#L151-L171
  _bxs.resize(3);
  _bxs << -1,0,1;
  _nPulseTot = _bxs.rows();

  _amplitudes = amplitudes;

  _pulseMat.resize(Eigen::NoChange,_nPulseTot);
  _ampVec = PulseVector::Zero(_nPulseTot);
  _errVec = PulseVector::Zero(_nPulseTot);

  _ampVec.coeffRef(0) = 0;//_amplitudes.coeff(3);
  _ampVec.coeffRef(1) = 0;//_amplitudes.coeff(4);
  _ampVec.coeffRef(2) = 0;//_amplitudes.coeff(5);

  _chiSq = 9999;

  aTaMat.resize(_nPulseTot, _nPulseTot);
  aTbVec.resize(_nPulseTot);
  wVec.resize(_nPulseTot);
  
  pulseShape = PulseVector::Zero(12);
  
  const double xx[4]={0.0, 1.0, 0.0, 3};
  (*pfunctor_)(&xx[0]);
  
  for (int i=-1; i<11; i++) {
    pulseShape.coeffRef(i+1) = psfPtr_->getPulseShape(i);
  }

  _pulseMat.col(0) = pulseShape.segment<10>(2);
  _pulseMat.col(1) = pulseShape.segment<10>(1);
  _pulseMat.col(2) = pulseShape.segment<10>(0);

  std::cout << "initial pulseMat" << std::endl;
  std::cout << _pulseMat << std::endl;

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

  std::cout << "------" << std::endl;
  std::cout << "input: " ;
  for (int i=0; i<10; i++) {
    std::cout << _amplitudes.coeff(i) << ", ";
  }
  std::cout << std::endl;
  //std::cout << ipulseintime << ", " << ipulseprevtime << ", " << ipulsenexttime << std::endl;

  std::cout << "output: ";
  std::cout << _ampVec.coeff(ipulseprevtime) << ", " << _ampVec.coeff(ipulseintime) << ", " << _ampVec.coeff(ipulsenexttime) << std::endl;

  //std::vector<double> ans;
  

  return status;
}

bool DoMahiAlgo::Minimize() {
  std::cout << "Minimize" << std::endl;

  int iter = 0;
  int maxIters = 500;
  bool status = false;

  while (true) {
    if (iter>=maxIters) {
      std::cout << "max number of iterations reached! " << std::endl;
      //std::cout << _chiSq << std::endl;
      break;
    }
    
    status=UpdateCov();
    if (!status) break;
    
    status = NNLS();
    if (!status) break;
    
    double newChiSq=CalculateChiSq();
    double deltaChiSq = newChiSq - _chiSq;
    
    _chiSq = newChiSq;
    
    //std::cout << "chiSq = " << _chiSq << ", " << deltaChiSq << std::endl;
    
    if (std::abs(deltaChiSq)<1e-3) break;
    
    iter++;
    
  }
  
  return status;
}

bool DoMahiAlgo::UpdateCov() {
  std::cout << "UpdateCov" << std::endl;

  const double pederr2 = 1;
  _invCovMat = pederr2*SampleMatrix::Constant(1);

  FullSampleVector pulseShapeM = FullSampleVector::Zero(12);
  FullSampleVector pulseShapeP = FullSampleVector::Zero(12);

  const double xxm[4]={-2.5, 1.0, 0.0, 3};
  const double xxp[4]={ 2.5, 1.0, 0.0, 3};
  (*pfunctor_)(&xxm[0]);
  for (int i=-1; i<11; i++) {
    pulseShapeM.coeffRef(i+1) = psfPtr_->getPulseShape(i);
  }

  (*pfunctor_)(&xxp[0]);
  for (int i=-1; i<11; i++) {
    pulseShapeP.coeffRef(i+1) = psfPtr_->getPulseShape(i);
  }

  for (int k=0; k<_ampVec.size(); k++) {

    double ifC=_ampVec.coeff(k);
    if (ifC==0) continue;

    //1-int(_bxs.coeff(i)

    int startV=1-int(_bxs.coeff(k));

    //std::cout << "startV: " << startV << std::endl;

    for (int i=startV; i<10+startV; i++) {
      for (int j=startV; j<i+1; j++) {
	//std::cout << i << ", " << j << std::endl;
	double tmp=0.5*((pulseShapeP.coeff(i)-pulseShape.coeff(i))*(pulseShapeP.coeff(j)-pulseShape.coeff(j)) 
			+ (pulseShapeM.coeff(i)-pulseShape.coeff(i))*(pulseShapeM.coeff(j)-pulseShape.coeff(j)));
	_invCovMat(i-startV,j-startV) += ifC*ifC*tmp;
	_invCovMat(j-startV,i-startV) += ifC*ifC*tmp;
	
      }
    }
    
    
  }
  
  _covDecomp.compute(_invCovMat);
  
  return true;
}

bool DoMahiAlgo::NNLS() {
  std::cout << "NNLS" << std::endl;

  const unsigned int npulse = _bxs.rows();
  
  std::cout << _ampVec << std::endl;
  for (uint i=0; i<npulse; i++) {
    //std::cout << double(_bxs.coeff(i)) << std::endl;
    //double nomT = _bxs.coeff(i)*25;
    //pulseShape=PulseVector::Zero(12);
    //getPulseShape(_ampVec.coeff(i), _detID, nomT, pulseShape,0);
    _pulseMat.col(i) = pulseShape.segment<10>(1-int(_bxs.coeff(i)));
  }
  
  std::cout << "new pulsemat" << std::endl;
  std::cout << _pulseMat << std::endl;
  
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
      std::cout << "nP..... " << _nP << std::endl;
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


double DoMahiAlgo::CalculateChiSq() {
  std::cout << "CalculateChiSq" << std::endl;

  return _covDecomp.matrixL().solve(_pulseMat*_ampVec - _amplitudes).squaredNorm();

  //return 0.0;
}

void DoMahiAlgo::setPulseShapeTemplate(const HcalPulseShapes::Shape& ps) {

  if (!(&ps == currentPulseShape_ ))
    {
      resetPulseShapeTemplate(ps);
      currentPulseShape_ = &ps;
    }
}

void DoMahiAlgo::resetPulseShapeTemplate(const HcalPulseShapes::Shape& ps) { 
  ++ cntsetPulseShape;
  psfPtr_.reset(new FitterFuncs::PulseShapeFunctor(ps,false,false,false,false,1,0,2.5,0,0.00065,1));
  pfunctor_    = std::unique_ptr<ROOT::Math::Functor>( new ROOT::Math::Functor(psfPtr_.get(),&FitterFuncs::PulseShapeFunctor::singlePulseShapeFunc, 3) );
  //p1functor_    = std::unique_ptr<ROOT::Math::Functor>( new ROOT::Math::Functor(psfPtr_.get(),&FitterFuncs::PulseShapeFunctor::singlePulseShapeFunc, 3) );
  //p2functor_    = std::unique_ptr<ROOT::Math::Functor>( new ROOT::Math::Functor(psfPtr_.get(),&FitterFuncs::PulseShapeFunctor::singlePulseShapeFunc, 3) );

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
