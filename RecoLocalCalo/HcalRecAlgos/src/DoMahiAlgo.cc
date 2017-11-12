#include "RecoLocalCalo/HcalRecAlgos/interface/DoMahiAlgo.h" 
#include <iostream>
#include <fstream> 

void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP);

DoMahiAlgo::DoMahiAlgo() { 
}

void DoMahiAlgo::setParameters(bool iApplyTimeSlew, HcalTimeSlew::BiasSetting slewFlavor,
			       double iMeanTime, double iTimeSigmaHPD, double iTimeSigmaSiPM,
			       const std::vector <int> &iActiveBXs, int iNMaxIters) {

  applyTimeSlew_ = iApplyTimeSlew;
  slewFlavor_    = slewFlavor;

  meanTime_      = iMeanTime;
  timeSigmaHPD_  = iTimeSigmaHPD;
  timeSigmaSiPM_ = iTimeSigmaSiPM;
  activeBXs_     = iActiveBXs;
  nMaxIters_     = iNMaxIters;
  nMaxItersNNLS_ = iNMaxIters;

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
  
  if (channelData.hasTimeInfo()) dt_=timeSigmaSiPM_;
  else dt_=timeSigmaHPD_;
  
  //Dark current value for this channel (SiPM only)
  double darkCurrent =  psfPtr_->getSiPMDarkCurrent(channelData.darkCurrent(), 
						    channelData.fcByPE(),
						    channelData.lambda());

  //Average pedestal width (for covariance matrix constraint)
  _pedConstraint = 0.25*( channelData.tsPedestalWidth(0)*channelData.tsPedestalWidth(0)+
			  channelData.tsPedestalWidth(1)*channelData.tsPedestalWidth(1)+
			  channelData.tsPedestalWidth(2)*channelData.tsPedestalWidth(2)+
			  channelData.tsPedestalWidth(3)*channelData.tsPedestalWidth(3) );

  if (channelData.hasTimeInfo()) _pedConstraint+=darkCurrent*darkCurrent;

  fcByPe_ = channelData.fcByPE();

  //_pedConstraint*=0.1;

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
    //double noiseDC=0;
    //if((channelData.hasTimeInfo()) && (charge-ped)>channelData.tsPedestalWidth(ip)) {
    //  noiseDC = darkCurrent;
    //}

    //Photostatistics
    //double noisePhoto = 0;
    //if ( (charge-ped)>channelData.tsPedestalWidth(ip)) {
    //noisePhoto = sqrt((charge-ped)*channelData.fcByPE());
    //}

    //Electronic pedestal
    double pedWidth = channelData.tsPedestalWidth(ip);

    //Total uncertainty from all sources
    //_noiseTerms.coeffRef(ip) = noiseADC*noiseADC + noiseDC*noiseDC + noisePhoto*noisePhoto + pedWidth*pedWidth;
    //_noiseTerms.coeffRef(ip) = noiseADC*noiseADC + noisePhoto*noisePhoto + pedWidth*pedWidth;
    _noiseTerms.coeffRef(ip) = noiseADC*noiseADC + pedWidth*pedWidth;

    tsTOT += charge - ped;
    if( ip==HcalConst::soi || ip==HcalConst::soi+1 ){
      tstrig += (charge - ped);
    }

  }

  //if (tsTOT>100) doDebug=2;

  //  if (doDebug==2) {
  //std::cout << "----debugging info" << std::endl;
  //std::cout << charges.transpose() << std::endl;
  //}

  bool status =false;
  if(tstrig >= 0) {
    //std::cout << "----debugging info" << std::endl;
    //std::cout << charges.transpose() << std::endl;
    //std::cout << _noiseTerms.transpose() << std::endl;
  
    status = DoFit(charges,reconstructedVals); 
    
  }
  //else {
  //  std::cout << "not calling mahi" << std::endl;
  //}
  
  if (!status) {
    if (doDebug==1) std::cout << "give up? " << std::endl;
    reconstructedVals.clear();
    reconstructedVals.push_back(0.);
    reconstructedVals.push_back(888.);
  }

  reconstructedEnergy = reconstructedVals[0]*channelData.tsGain(0);
  chi2 = reconstructedVals[1];

}

bool DoMahiAlgo::DoFit(SampleVector amplitudes, std::vector<float> &correctedOutput) {
  if (doDebug==1) std::cout << "DoFit" << std::endl;

  _bxs.resize(activeBXs_.size());

  for (uint ibx=0; ibx<activeBXs_.size(); ibx++) {
    _bxs.coeffRef(ibx) = activeBXs_[ibx];
  }

  _nPulseTot = _bxs.rows();  
  _amplitudes = amplitudes;

  _pulseMat.resize(Eigen::NoChange,_nPulseTot);
  _ampVec = PulseVector::Zero(_nPulseTot);
  _errVec = PulseVector::Zero(_nPulseTot);

  _nP=0;

  bool status = true;

  pulseShapeArray.clear();
  pulseCovArray.clear();
  mapBXs.clear();

  //pulseOffset_ = -_bxs.minCoeff();
  pulseOffset_=4;
  int OFFSET=0;
  
  for (uint i=0; i<_nPulseTot; i++) {
    
    OFFSET=_bxs.coeff(i);
    mapBXs[OFFSET] = i;

    pulseShapeArray.push_back(FullSampleVector::Zero(FullSampleVectorSize));
    pulseCovArray.push_back(FullSampleMatrix::Constant(0));

    status = UpdatePulseShape(_amplitudes.coeff(HcalConst::soi+OFFSET), 
			      pulseShapeArray.at(i), 
			      pulseCovArray.at(i));


    if (OFFSET==0) {
      //double ampCorrection = 1.0 + double(pulseShapeArray.at(i).coeff(HcalConst::soi+5));
      _ampVec.coeffRef(i)= _amplitudes.coeff(HcalConst::soi+OFFSET);//*ampCorrection;
    }
    else {
      _ampVec.coeffRef(i)=0;
    }

    int SEGSTART=std::max(pulseOffset_ - OFFSET, 0);
    _pulseMat.col(i) = pulseShapeArray.at(i).segment<HcalConst::maxSamples>(SEGSTART);
  }

  _chiSq = 9999;
  
  aTaMat.resize(_nPulseTot, _nPulseTot);
  aTbVec.resize(_nPulseTot);
  wVec.resize(_nPulseTot);

  status = Minimize(); 
  _ampVecMin = _ampVec;
  _bxsMin = _bxs;

  if (!status) return status;

  bool foundintime = false;
  unsigned int ipulseintime = 4;

  for (unsigned int ipulse=0; ipulse<_nPulseTot; ++ipulse) {
    if (_bxs.coeff(ipulse)==0) {
      ipulseintime = ipulse;
      foundintime = true;
      //std::cout << "intime: " << _ampVec.coeff(ipulse) << " vs " << _amplitudes.coeff(4)<< "; " << std::endl;
    }
    //else if (_bxs.coeff(ipulse)==10) {
    //  std::cout << "pedestal: " << _ampVec.coeff(ipulse) << ", ";
    //else if (_bxs.coeff(ipulse)==-1) {
    //  std::cout << "prev: " << _ampVec.coeff(ipulse) << ", ";
    //}
    //else if (_bxs.coeff(ipulse)==1) {
    //  std::cout << "next: " << _ampVec.coeff(ipulse) << ", ";
    //}
  }
  //std::cout << std::endl;
  if (!foundintime) return status;

  correctedOutput.clear();
  correctedOutput.push_back(_ampVec.coeff(ipulseintime)); //charge
  correctedOutput.push_back(_chiSq); //chi2
  
  return status;
}

bool DoMahiAlgo::Minimize() {

  int iter = 0;
  //int maxIters = 500;
  bool status = false;

  while (true) {
    if (iter>=nMaxIters_) {
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
    
    //if (doDebug==1) 
    //std::cout << "chiSq = " << _chiSq << ", " << deltaChiSq << std::endl;

    if (std::abs(deltaChiSq)<1e-3) break;
    
    iter++;
    
  }
  
  //if (doDebug==1) 
  //std::cout << "miniter " << iter << std::endl;
  
  return status;
}

bool DoMahiAlgo::UpdatePulseShape(double itQ, FullSampleVector &pulseShape, FullSampleMatrix &pulseCov) {
  if (doDebug==1) std::cout << "UpdatePulseShape" << std::endl;

  float t0=meanTime_;
  if (applyTimeSlew_) 
    t0+=HcalTimeSlew::delay(std::max(1.0, itQ), slewFlavor_);


  const double xx[4]={t0, 1.0, 0.0, 3};
  (*pfunctor_)(&xx[0]);

  double normP=psfPtr_->getPulseShape(4);

  for (int i=-4; i<15; i++) {
    pulseShape.coeffRef(i+4) = psfPtr_->getPulseShape(i)/normP;
  }

  FullSampleVector pulseShapeM = FullSampleVector::Zero(FullSampleVectorSize);
  FullSampleVector pulseShapeP = FullSampleVector::Zero(FullSampleVectorSize);

  const double xxm[4]={t0-dt_, 1.0, 0.0, 3};
  const double xxp[4]={t0+dt_, 1.0, 0.0, 3};

  if (doDebug==1) {
    std::cout << itQ << ", " << t0 << ", " << -dt_+t0 << ", " << dt_+t0 << std::endl;
  }

  (*pfunctor_)(&xxm[0]);
  normP=psfPtr_->getPulseShape(4);
  for (int i=-4; i<15; i++) {
    pulseShapeM.coeffRef(i+4) = psfPtr_->getPulseShape(i)/normP;
  }
  (*pfunctor_)(&xxp[0]);
  normP=psfPtr_->getPulseShape(4);
  for (int i=-4; i<15; i++) {
    pulseShapeP.coeffRef(i+4) = psfPtr_->getPulseShape(i)/normP;
  }
  for (int i=0; i<19; i++) {
    for (int j=0; j<i+1; j++) {
      
      double tmp=0.5*((pulseShapeP.coeff(i)-pulseShape.coeff(i))*(pulseShapeP.coeff(j)-pulseShape.coeff(j))
		      + (pulseShapeM.coeff(i)-pulseShape.coeff(i))*(pulseShapeM.coeff(j)-pulseShape.coeff(j)));
      
      pulseCov(i,j) += tmp;
      pulseCov(j,i) += tmp;
      
    }
  }


  return true;  
}



bool DoMahiAlgo::UpdateCov() {
  if (doDebug==1) std::cout << "UpdateCov" << std::endl;

  bool status=true;

  _invCovMat = _noiseTerms.asDiagonal();
  _invCovMat += _pedConstraint*SampleMatrix::Ones();
  //
  for (int k=0; k< _ampVec.size(); k++) {
    if (_ampVec.coeff(k)==0) continue;
    
    int OFFSET=_bxs.coeff(k);
    _invCovMat += _ampVec.coeff(k)*_ampVec.coeff(k)*pulseCovArray.at(mapBXs[OFFSET]).block(4-OFFSET, 4-OFFSET, HcalConst::maxSamples, HcalConst::maxSamples);
    _invCovMat += _ampVec.coeff(k)*_ampVec.coeff(k)*fcByPe_*fcByPe_*SampleMatrix::Ones();
  }

  //if (doDebug==1) {
  //std::cout << "cov" << std::endl;
  //std::cout << _invCovMat << std::endl;
  //std::cout << "..." << std::endl;
  //}
  
  _covDecomp.compute(_invCovMat);
  
  return status;
}

bool DoMahiAlgo::NNLS() {
  //if (doDebug==1) 
  //std::cout << "NNLS" << std::endl;

  const unsigned int npulse = _bxs.rows();
  if (doDebug==1) {
    std::cout << "_ampVec" << std::endl;
    std::cout << _ampVec << std::endl;
  }

  for (uint i=0; i<npulse; i++) {
    int OFFSET=_bxs.coeff(i);
    if (OFFSET==10) {  
      _pulseMat.col(i) = pulseShapeArray.at(mapBXs[OFFSET]).segment<HcalConst::maxSamples>(0);
    }
    else {
      int SEGSTART=std::max(pulseOffset_-OFFSET, 0);
      _pulseMat.col(i) = pulseShapeArray.at(mapBXs[OFFSET]).segment<HcalConst::maxSamples>(SEGSTART);
    }
  }

  //if (doDebug==1) {
  //std::cout << "----- Problem statement: " << std::endl;
  //std::cout << "new ampvec " << std::endl;
  //std::cout << _ampVec.transpose() << std::endl;
  //std::cout << "new pulsemat" << std::endl;
  //std::cout << _ampVec.transpose()*_pulseMat.transpose() << std::endl; //*_ampVec.coeff(0) << std::endl;
  //std::cout << "ampvec" << std::endl;
  //std::cout << _ampVec.transpose() << std::endl;
  //}
  
  invcovp = _covDecomp.matrixL().solve(_pulseMat);
  aTaMat = invcovp.transpose()*invcovp;
  aTbVec = invcovp.transpose()*_covDecomp.matrixL().solve(_amplitudes);
  
  int iter = 0;
  Index idxwmax = 0;
  double wmax = 0.0;
  double threshold = 1e-11;

  while (true) {    
    //std::cout << "while loop " << std::endl;
    if (iter>0 || _nP==0) {
      //std::cout << "meep?" << std::endl;
      if ( _nP==npulse ) break;
      
      const unsigned int nActive = npulse - _nP;
      //std::cout << "nActive = " << nActive << std::endl;
      updateWork = aTbVec - aTaMat*_ampVec;
      //std::cout << "updateWork " << updateWork.transpose() << std::endl;

      Index idxwmaxprev = idxwmax;
      double wmaxprev = wmax;
      wmax = updateWork.tail(nActive).maxCoeff(&idxwmax);
      
      if (wmax<threshold || (idxwmax==idxwmaxprev && wmax==wmaxprev)) {
	//std::cout << "converged" << std::endl;
	break;
      }
      
      if (iter>=nMaxItersNNLS_) {
	//std::cout << "Max Iterations reached!" << std::endl;
	break;
      }

      //unconstrain parameter
      Index idxp = _nP + idxwmax;

      //if (iter==0) std::cout << "BXS SWAP " << _bxs.coeffRef(_nP) << ", " << _bxs.coeffRef(idxp) << std::endl;
      
      aTaMat.col(_nP).swap(aTaMat.col(idxp));
      aTaMat.row(_nP).swap(aTaMat.row(idxp));
      _pulseMat.col(_nP).swap(_pulseMat.col(idxp));
      std::swap(aTbVec.coeffRef(_nP),aTbVec.coeffRef(idxp));
      std::swap(_ampVec.coeffRef(_nP),_ampVec.coeffRef(idxp));
      std::swap(_bxs.coeffRef(_nP),_bxs.coeffRef(idxp));

      wVec.tail(nActive) = updateWork.tail(nActive); 
      //std::cout << "????" << std::endl;
      ++_nP;

    }

    //std::cout << "between the wires " << _nP << std::endl;

    while (true) {
      if (_nP==0) break;     

      ampvecpermtest = _ampVec;

      eigen_solve_submatrix(aTaMat,aTbVec,ampvecpermtest,_nP);
      //std::cout << ampvecpermtest.transpose() << std::endl;
      //check solution    
      //std::cout << "nP..... " << _nP << std::endl;
      //auto ampvecpermhead = ampvecpermtest.head(_nP);

      //if ( ampvecpermhead.minCoeff()>0. ) {
      //	//std::cout << "ampvecpermhead" << std::endl;
      //	std::cout << ampvecpermhead.minCoeff() << " ???" << std::endl;
      //	_ampVec.head(_nP) = ampvecpermhead.head(_nP);
      //	break;
      //}

      //check solution
      bool positive = true;
      for (unsigned int i = 0; i < _nP; ++i)
        positive &= (ampvecpermtest(i) > 0);
      if (positive) {
	//std::cout << "idek man, " << _nP << std::endl;
        _ampVec.head(_nP) = ampvecpermtest.head(_nP);
        break;
      } 

      //if (iter!=0) std::cout << "wow got to nnls part!" << std::endl;

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
      _ampVec.head(_nP) += minratio*(ampvecpermtest.head(_nP) - _ampVec.head(_nP));
      
      //avoid numerical problems with later ==0. check
      _ampVec.coeffRef(minratioidx) = 0.;
      
      aTaMat.col(_nP-1).swap(aTaMat.col(minratioidx));
      aTaMat.row(_nP-1).swap(aTaMat.row(minratioidx));
      _pulseMat.col(_nP-1).swap(_pulseMat.col(minratioidx));
      std::swap(aTbVec.coeffRef(_nP-1),aTbVec.coeffRef(minratioidx));
      std::swap(_ampVec.coeffRef(_nP-1),_ampVec.coeffRef(minratioidx));
      std::swap(_bxs.coeffRef(_nP-1),_bxs.coeffRef(minratioidx));
      --_nP;

      //std::cout << "end of second looppppp " << _nP << std::endl;

    }
   
    ++iter;

    //adaptive convergence threshold to avoid infinite loops but still
    //ensure best value is used
    if (iter%10==0) {
      threshold *= 10.;
    }

    //std::cout << "wtf" << std::endl;

    break;
  }
  //std::cout << "nnlsiter " << iter << std::endl;
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
