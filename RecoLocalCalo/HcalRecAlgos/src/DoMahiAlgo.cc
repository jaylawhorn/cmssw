#include "RecoLocalCalo/HcalRecAlgos/interface/DoMahiAlgo.h" 
#include <iostream>
#include <fstream> 


void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP);

DoMahiAlgo::DoMahiAlgo() { 
}

void DoMahiAlgo::setParameters(bool iDoPrefit, bool iFloatPedestal, bool iApplyTimeSlew, HcalTimeSlew::BiasSetting slewFlavor,
			       double iMeanTime, double iTimeSigmaHPD, double iTimeSigmaSiPM,
			       const std::vector <int> &iActiveBXs, int iNMaxIters) {

  doPrefit_      = iDoPrefit;
  floatPedestal_ = iFloatPedestal;
  applyTimeSlew_ = iApplyTimeSlew;
  slewFlavor_    = slewFlavor;

  meanTime_      = iMeanTime;
  timeSigmaHPD_  = iTimeSigmaHPD;
  timeSigmaSiPM_ = iTimeSigmaSiPM;
  //activeBXs_     = iActiveBXs;
  nMaxIters_     = iNMaxIters;
  nMaxItersNNLS_ = iNMaxIters;

  //ECAL does it better -- to be fixed
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#L151-L171

  //if (floatPedestal_) 
  //  _bxs.resize(iActiveBXs.size()+1);
  //else 
  _bxs.resize(iActiveBXs.size()+1);

  std::cout << iActiveBXs.size()+1 << std::endl;

  for (uint ibx=0; ibx<iActiveBXs.size(); ibx++) {
    _bxs.coeffRef(ibx) = iActiveBXs[ibx];
    std::cout << iActiveBXs[ibx] << ", ";
  }
  std::cout << std::endl;

  //if (floatPedestal_) 
  _bxs.coeffRef(iActiveBXs.size()) = 10; // for pedestal val

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
  
  if (channelData.hasTimeInfo()) {
    isHPD=false;
    dt_=timeSigmaSiPM_;
  }
  else {
    isHPD=true;
    dt_=timeSigmaHPD_;
  }

  //Dark current value for this channel (SiPM only)
  double darkCurrent =  psfPtr_->getSiPMDarkCurrent(channelData.darkCurrent(), 
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
    if((!isHPD) && (charge-ped)>channelData.tsPedestalWidth(ip)) {
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
    if( ip==HcalConst::soi || ip==HcalConst::soi+1 ){
      tstrig += (charge - ped);//*channelData.tsGain(4);
    }
  }

  if (doDebug>0) {
    std::cout << "----debugging info" << std::endl;
    std::cout << charges << std::endl;
  }

  if( tstrig*channelData.tsGain(0)>5 && tstrig*channelData.tsGain(0)<10) doDebug=2;

  bool status =false;
  if(tstrig >= 0) {
    if (doDebug>0) std::cout << "three pulse fit " << std::endl;
    
    status = DoFit(charges,reconstructedVals); 
    
  }
  
  if (!status) {
    if (doDebug>0) std::cout << "give up? " << std::endl;
    reconstructedVals.clear();
    reconstructedVals.push_back(0.);
    reconstructedVals.push_back(888.);
  }
  
  reconstructedEnergy = reconstructedVals[0]*channelData.tsGain(0);
  chi2 = reconstructedVals[1];

}

bool DoMahiAlgo::DoFit(SampleVector amplitudes, std::vector<float> &correctedOutput) {
  if (doDebug>0) std::cout << "DoFit" << std::endl;

  _nP = 0;

  _nPulseTot = _bxs.rows();

  _amplitudes = amplitudes;

  _pulseMat.resize(Eigen::NoChange,_nPulseTot);
  _ampVec = PulseVector::Zero(_nPulseTot);
  _errVec = PulseVector::Zero(_nPulseTot);

  bool status = true;

  status=CalculatePulseShapes();

  for (uint i=0; i<_nPulseTot; i++) {

    if (int(_bxs.coeff(i))==10) _pulseMat.col(i) = SampleVector::Constant(1.0/float(HcalConst::maxSamples));
    else _pulseMat.col(i) = pulseShapeArray.at( mapBXs[int(_bxs.coeff(i))] ).segment<HcalConst::maxSamples>(int(_bxs.coeff(i))+1);


    if (int(_bxs.coeff(i))==0) _ampVec.coeffRef(i) = _amplitudes.coeff(HcalConst::soi) / pulseShapeArray.at(mapBXs[int(_bxs.coeff(i))]).coeff(HcalConst::soi+1);
    else if (int(_bxs.coeff(i))==10) _ampVec.coeffRef(i)=1.5;
    else _ampVec.coeffRef(i) =0;
  }

  if (doDebug>0) {
    std::cout << std::endl << "pulseMat : " << std::endl;
    std::cout << _pulseMat << std::endl << std::endl;
    std::cout << "ampVec : " << std::endl;
    std::cout << _ampVec << std::endl;
  }


  _chiSq = 999;

  aTaMat.resize(_nPulseTot, _nPulseTot);
  aTbVec.resize(_nPulseTot);
  wVec.resize(_nPulseTot);

  status = Minimize(); 
  _ampVecMin = _ampVec;
  _bxsMin = _bxs;

  if (!status) return status;

  bool foundintime = false;
  unsigned int ipulseintime = 10;
  unsigned int ipulseprevtime = 10;
  unsigned int ipulsenexttime = 10;
  unsigned int iped=10;

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
    else if (_bxs.coeff(ipulse)==10) {
      iped = ipulse;
    }
  }
  if (!foundintime) return status;

  if (doDebug==2) {

    std::cout << "inside mahi: ";
    if (ipulseprevtime!=10) std::cout << "(prev) " << _ampVec.coeff(ipulseprevtime) << ", ";
    if (ipulseintime!=10)   std::cout << "(in  ) " << _ampVec.coeff(ipulseintime) << ", ";
    if (ipulsenexttime!=10) std::cout << "(next) " << _ampVec.coeff(ipulsenexttime) << ", ";
    if (iped!=10) std::cout << "(ped ) " << _ampVec.coeff(iped) << ", ";
    std::cout << std::endl;
  }

  correctedOutput.clear();
  correctedOutput.push_back(_ampVec.coeff(ipulseintime)); //charge
  correctedOutput.push_back(_chiSq); //chi2
  
  
  return status;
}

bool DoMahiAlgo::Minimize() {

  int iter = 0;
  int maxIters = nMaxIters_;
  bool status = false;

  while (true) {
    if (iter>=maxIters) {
      std::cout << "max number of iterations reached! " << std::endl;
      break;
    }
    
    status=UpdateCov();
    if (!status) break;
    
    status = NNLS();
    if (!status) break;
    
    double newChiSq=CalculateChiSq();
    double deltaChiSq = newChiSq - _chiSq;
    
    _chiSq = newChiSq;
    
    if (doDebug>0) std::cout << "chiSq = " << _chiSq << ", " << deltaChiSq << std::endl;
    
    if (std::abs(deltaChiSq)<1e-3) break;
    
    iter++;
    
  }
  
  return status;
}

bool DoMahiAlgo::CalculatePulseShapes() {
  if (doDebug>0) std::cout << "CalculatePulseShapes" << std::endl;

  for (uint ibx = 0; ibx<_nPulseTot; ibx++) {
    
    mapBXs[int(_bxs.coeff(ibx))] = ibx;
    if (doDebug>0) std::cout << "mapBXs[" << int(_bxs.coeff(ibx)) << "] = " << ibx << std::endl;
    //pulseShapeArray.push_back(FullSampleVector::Zero(FullSampleVectorSize));

    if (_bxs.coeff(ibx) == 10 ) {
      //pedestal
      pulseCovArray.push_back(FullSampleMatrix::Constant(_pedConstraint));
      pulseShapeArray.push_back(FullSampleVector::Constant(1/float(HcalConst::maxSamples)));
    }
    else {
      //real pulse
      pulseCovArray.push_back(FullSampleMatrix::Constant(0));
      pulseShapeArray.push_back(FullSampleVector::Zero(12));

      double xx[4]={meanTime_, 1.0, 0.0, 3};
      if (applyTimeSlew_) 
	xx[0]+=HcalTimeSlew::delay(std::max(1.0, _amplitudes.coeff(_bxs.coeff(ibx)+4)), slewFlavor_);

      (*pfunctor_)(&xx[0]);

      for (int i=-1; i<11; i++) {
	pulseShapeArray.at(ibx).coeffRef(i+1) = psfPtr_->getPulseShape(i);
      }

      FullSampleVector pulseShapeM = FullSampleVector::Zero(12);
      FullSampleVector pulseShapeP = FullSampleVector::Zero(12);

      xx[0]=meanTime_ - dt_;
      if (applyTimeSlew_) 
	xx[0]+=HcalTimeSlew::delay(std::max(1.0, _amplitudes.coeff(ibx+4)), slewFlavor_);

      (*pfunctor_)(&xx[0]);
      for (int i=-1; i<11; i++) {
	pulseShapeM.coeffRef(i+1) = psfPtr_->getPulseShape(i);
      }      

      xx[0]=meanTime_ + dt_;
      if (applyTimeSlew_) 
	xx[0]+=HcalTimeSlew::delay(std::max(1.0, _amplitudes.coeff(ibx+4)), slewFlavor_);

      (*pfunctor_)(&xx[0]);
      for (int i=-1; i<11; i++) {
	pulseShapeP.coeffRef(i+1) = psfPtr_->getPulseShape(i);
      }

      for (int i=0; i<11; i++) {
	for (int j=0; j<i+1; j++) {
	  
	  double tmp=0.5*((pulseShapeP.coeff(i)-pulseShapeArray.at(ibx).coeff(i))*(pulseShapeP.coeff(j)-pulseShapeArray.at(ibx).coeff(j))
			  + (pulseShapeM.coeff(i)-pulseShapeArray.at(ibx).coeff(i))*(pulseShapeM.coeff(j)-pulseShapeArray.at(ibx).coeff(j)));
	  
	  pulseCovArray.at(ibx)(i,j) += tmp;
	  pulseCovArray.at(ibx)(j,i) += tmp;
	  
	}
      }
      //if (doDebug>0) std::cout << "----------- " << ibx << std::endl;
      //if (doDebug>0) std::cout << pulseCovArray.at(ibx) << std::endl;

    }
    
  }
  
  return true;  
}



bool DoMahiAlgo::UpdateCov() {
  if (doDebug>0) std::cout << "UpdateCov" << std::endl;

  bool status=true;

  _invCovMat = _noiseTerms.asDiagonal();

  SampleMatrix tempInvCov = SampleMatrix::Constant(_pedConstraint);
  
  for (int k=0; k< _ampVec.size(); k++) {

    int ibx=int(_bxs.coeff(k));

    //if (doDebug>0) std::cout << "_ampVec.coeff(k)*_ampVec.coeff(k)*pulseCovArray.at(mapBXs[ibx]).block(1-ibx, 1-ibx, HcalConst::maxSamples, HcalConst::maxSamples)" << std::endl;
    if (ibx==10) tempInvCov+= _ampVec.coeff(k)*_ampVec.coeff(k)*SampleMatrix::Constant(_pedConstraint);
    else 
      tempInvCov+= _ampVec.coeff(k)*_ampVec.coeff(k)*pulseCovArray.at(mapBXs[ibx]).block(1-ibx, 1-ibx, HcalConst::maxSamples, HcalConst::maxSamples);

    //if (doDebug>0) std::cout << tempInvCov << std::endl << std::endl;

  }
  
  _invCovMat+=tempInvCov;  
  
  if (doDebug>0) {
    std::cout << "cov" << std::endl;
    std::cout << _invCovMat << std::endl;
    std::cout << "..." << std::endl;
  }
  
  _covDecomp.compute(_invCovMat);

  return status;
}

bool DoMahiAlgo::NNLS() {
  if (doDebug>0) std::cout << "NNLS" << std::endl;

  const unsigned int npulse = _bxs.rows();
  
  if (doDebug>0) {
    std::cout << "_ampVec" << std::endl;
    std::cout << _ampVec << std::endl;
  }

  //for (uint i=0; i<npulse; i++) {
  //  if (int(_bxs.coeff(i))==10) _pulseMat.col(i) = SampleVector::Constant(1.0/float(HcalConst::maxSamples));
  //  else _pulseMat.col(i) = pulseShapeArray.at( mapBXs[int(_bxs.coeff(i))] ).segment<HcalConst::maxSamples>(int(_bxs.coeff(i))+1);
  //}

  if (doDebug>0) {
    std::cout << "new pulsemat" << std::endl;
    std::cout << _pulseMat << std::endl;
    std::cout << "ampvec" << std::endl;
    std::cout << _ampVec << std::endl;
  }
  
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
      
      if (iter>=nMaxItersNNLS_) {
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
