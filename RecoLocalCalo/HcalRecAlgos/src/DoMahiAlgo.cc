#include "RecoLocalCalo/HcalRecAlgos/interface/DoMahiAlgo.h" 
#include <iostream>
#include <fstream> 

void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP);

DoMahiAlgo::DoMahiAlgo() { 
}

void DoMahiAlgo::setDebug(int val) {
  doDebug=val;
  if (doDebug== 1) std::cout << "print debugging info" << std::endl;

}

void DoMahiAlgo::phase1Apply(const HBHEChannelInfo& channelData,
			     float& reconstructedEnergy,
			     float& chi2) {
  
  const unsigned cssize = channelData.nSamples();
  _detID = channelData.id();

  if (channelData.hasTimeInfo()) doDynamicPulseCov=false;
  else doDynamicPulseCov=true;

  //Dark current value for this channel
  double darkCurrent = psfPtr_->getSiPMDarkCurrent(channelData.darkCurrent(), 
						   channelData.fcByPE(),
						   channelData.lambda());

  //Average pedestal width (for covariance matrix constraint)
  _pedConstraint = 0.25*( channelData.tsPedestalWidth(0)*channelData.tsPedestalWidth(0)+
			  channelData.tsPedestalWidth(1)*channelData.tsPedestalWidth(1)+
			  channelData.tsPedestalWidth(2)*channelData.tsPedestalWidth(2)+
			  channelData.tsPedestalWidth(3)*channelData.tsPedestalWidth(3) );

  std::vector<float> reconstructedVals;
  SampleVector charges;
  
  double tsTOT = 0, tstrig = 0; // in fC
  for(unsigned int ip=0; ip<cssize; ++ip){
    if( ip >= (unsigned)HcalConst::maxSamples) continue; 
    double charge = channelData.tsRawCharge(ip);
    double ped = channelData.tsPedestal(ip);

    charges.coeffRef(ip) = charge - ped;

    //ADC granularity
    double noiseADC = (1./sqrt(12))*channelData.tsDFcPerADC(ip);

    //Dark current (for SiPMs)
    double noiseDC=0;
    if(channelData.hasTimeInfo() && (charge-ped)>channelData.tsPedestalWidth(ip)) {
      noiseDC = darkCurrent;
    }

    //Photostatistics
    double noisePhoto = 0;
    if ( (charge-ped)>channelData.tsPedestalWidth(ip)) {
      noisePhoto = sqrt((charge-ped)*channelData.fcByPE());
    }

    //Electronic pedestal
    double pedWidth = channelData.tsPedestalWidth(ip);

    //Total uncertainty from all sources
    _noiseTerms.coeffRef(ip) = noiseADC*noiseADC + noiseDC*noiseDC + noisePhoto*noisePhoto + pedWidth*pedWidth;

    tsTOT += charge - ped;
    if( ip ==4 || ip==5 ){
      tstrig += (charge - ped)*channelData.tsGain(4);
    }
  }

  if (doDebug==1) {
    std::cout << "----debugging info" << std::endl;
    std::cout << charges << std::endl;
  }

  bool status =false;
  if(tstrig >= 0) {
    //if (doDebug==1) std::cout << "one pulse fit " << std::endl;
    //status = DoFit(charges,reconstructedVals,1); 

    //if (reconstructedVals[1]>15 && tstrig < 45000) {

    if (doDebug==1) std::cout << "three pulse fit " << std::endl;
    //reconstructedVals.clear();
    
    status = DoFit(charges,reconstructedVals,3); 
    
    //}
    
  }
  
  if (!status) {
    if (doDebug==1) std::cout << "give up? " << std::endl;
    reconstructedVals.clear();
    reconstructedVals.push_back(0.);
    reconstructedVals.push_back(888.);
  }
  
  reconstructedEnergy = reconstructedVals[0]*channelData.tsGain(0);
  chi2 = reconstructedVals[1];

}

bool DoMahiAlgo::DoFit(SampleVector amplitudes, std::vector<float> &correctedOutput, int nbx) {
  if (doDebug==1) std::cout << "DoFit" << std::endl;

  _nP = 0;
  //ECAL does it better -- to be fixed
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#L151-L171
  
  if (nbx==3) {
    _bxs.resize(3);
    _bxs << -1,0,1;
  }
  else {
    _bxs.resize(1);
    _bxs << 0;
  }

  _nPulseTot = _bxs.rows();  

  _amplitudes = amplitudes;

  _pulseMat.resize(Eigen::NoChange,_nPulseTot);
  _ampVec = PulseVector::Zero(_nPulseTot);
  _errVec = PulseVector::Zero(_nPulseTot);

  //need to fix
  if (_nPulseTot==1) {
    _ampVec.coeffRef(0) = _amplitudes.coeff(4)*1.3;
  }
  else {
    _ampVec.coeffRef(0) = _amplitudes.coeff(3)*1.3;
    _ampVec.coeffRef(1) = _amplitudes.coeff(4)*1.3;
    _ampVec.coeffRef(2) = _amplitudes.coeff(5)*1.3;
  }

  _chiSq = 999;

  aTaMat.resize(_nPulseTot, _nPulseTot);
  aTbVec.resize(_nPulseTot);
  wVec.resize(_nPulseTot);

  bool status = UpdatePulseShape();

  status = Minimize(); 
  _ampVecMin = _ampVec;
  _bxsMin = _bxs;

  if (!status) return status;

  bool foundintime = false;
  unsigned int ipulseintime = 0;
  //unsigned int ipulseprevtime = 0;
  //unsigned int ipulsenexttime = 0;
  //
  for (unsigned int ipulse=0; ipulse<_nPulseTot; ++ipulse) {
    if (_bxs.coeff(ipulse)==0) {
      ipulseintime = ipulse;
      foundintime = true;
    }
  }
  //  else if (_bxs.coeff(ipulse)==-1) {
  //    ipulseprevtime = ipulse;
  //  }
  //  else if (_bxs.coeff(ipulse)==1) {
  //    ipulsenexttime = ipulse;
  //  }
  //}
  if (!foundintime) return status;

  correctedOutput.clear();
  correctedOutput.push_back(_ampVec.coeff(ipulseintime)); //charge
  correctedOutput.push_back(_chiSq); //chi2

  
  return status;
}

bool DoMahiAlgo::Minimize() {
  //std::cout << "Minimize" << std::endl;

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
    
    if (doDebug==1) std::cout << "chiSq = " << _chiSq << ", " << deltaChiSq << std::endl;
    
    if (std::abs(deltaChiSq)<1e-3) break;
    
    iter++;
    
  }
  
  return status;
}

bool DoMahiAlgo::UpdatePulseShape() {
  if (doDebug==1) std::cout << "UpdatePulseShape" << std::endl;

  pulseCov = FullSampleMatrix::Constant(0);

  pulseShape = PulseVector::Zero(12);

  //FullSampleVector pulseShapeM = FullSampleVector::Zero(12);
  //FullSampleVector pulseShapeP = FullSampleVector::Zero(12);

  const double xx[4]={0.0, 1.0, 0.0, 3};
  (*pfunctor_)(&xx[0]);
  
  for (int i=-1; i<11; i++) {
    pulseShape.coeffRef(i+1) = psfPtr_->getPulseShape(i);
  }

  if (_nPulseTot==1) {
    _pulseMat.col(0) = pulseShape.segment<10>(1);
  }

  else {
    _pulseMat.col(0) = pulseShape.segment<10>(2);
    _pulseMat.col(1) = pulseShape.segment<10>(1);
    _pulseMat.col(2) = pulseShape.segment<10>(0);
  }

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
  
  for (int i=0; i<11; i++) {
    for (int j=0; j<i+1; j++) {
      
      double tmp=0.5*((pulseShapeP.coeff(i)-pulseShape.coeff(i))*(pulseShapeP.coeff(j)-pulseShape.coeff(j))
		      + (pulseShapeM.coeff(i)-pulseShape.coeff(i))*(pulseShapeM.coeff(j)-pulseShape.coeff(j)));

      pulseCov(i,j) += tmp;
      pulseCov(j,i) += tmp;
      
    }
  }

  //if (doDebug==1) std::cout << "from updatepulseshape" << std::endl;
  //if (doDebug==1) std::cout << pulseCov.block(2,2,10,10) << std::endl;
  //
  //if (doDebug==1) std::cout << "from updatepulseshape2" << std::endl;
  //if (doDebug==1) std::cout << pulseCov << std::endl;

  return true;  
}



bool DoMahiAlgo::UpdateCov() {
  if (doDebug==1) std::cout << "UpdateCov" << std::endl;

  bool status=true;

  _invCovMat = _noiseTerms.asDiagonal();
  _invCovMat +=SampleMatrix::Constant(_pedConstraint);

  if (doDynamicPulseCov) status = UpdatePulseShape();

  for (int k=0; k< _ampVec.size(); k++) {
    _invCovMat += _ampVec.coeff(k)*_ampVec.coeff(k)*pulseCov.block(1-_bxs.coeff(k), 1-_bxs.coeff(k),10,10);
  }

  if (doDebug==1) {
    std::cout << "cov" << std::endl;
    std::cout << _invCovMat << std::endl;
    std::cout << "..." << std::endl;
  }
  
  _covDecomp.compute(_invCovMat);
  
  return true;
}

bool DoMahiAlgo::NNLS() {
  if (doDebug==1) std::cout << "NNLS" << std::endl;

  const unsigned int npulse = _bxs.rows();
  
  if (doDebug==1) {
    std::cout << "_ampVec" << std::endl;
    std::cout << _ampVec << std::endl;
  }

  for (uint i=0; i<npulse; i++) {
    _pulseMat.col(i) = pulseShape.segment<10>(1-int(_bxs.coeff(i)));
  }

  //if (doDebug==1) {
  //  std::cout << "new pulsemat" << std::endl;
  //  std::cout << _pulseMat << std::endl;
  //  std::cout << "ampvec" << std::endl;
  //  std::cout << _ampVec << std::endl;
  //}
  
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


double DoMahiAlgo::CalculateChiSq() {
  //std::cout << "CalculateChiSq" << std::endl;

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

}

void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP) {
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
