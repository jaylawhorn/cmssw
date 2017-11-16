#include "RecoLocalCalo/HcalRecAlgos/interface/DoMahiAlgo.h" 
#include <iostream>
#include <fstream> 

void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP);

DoMahiAlgo::DoMahiAlgo() :
  TSSize_(10),
  TSOffset_(4),
  FullTSSize_(19),
  FullTSofInterest_(8),
  FullTSOffset_(4) // FullTSofInterest - TSOffset
{}

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

  BXOffset_ = -(*std::min_element(activeBXs_.begin(), activeBXs_.end()));
  BXSize_   = activeBXs_.size();
}

void DoMahiAlgo::phase1Apply(const HBHEChannelInfo& channelData,
			     float& reconstructedEnergy,
			     float& chi2) {

  if (channelData.hasTimeInfo()) dt_=timeSigmaSiPM_;
  else dt_=timeSigmaHPD_;
  
  //Dark current value for this channel (SiPM only)
  float darkCurrent =  psfPtr_->getSiPMDarkCurrent(channelData.darkCurrent(), 
						   channelData.fcByPE(),
						   channelData.lambda());

  //Average pedestal width (for covariance matrix constraint)
  pedConstraint_ = 0.25*( channelData.tsPedestalWidth(0)*channelData.tsPedestalWidth(0)+
			  channelData.tsPedestalWidth(1)*channelData.tsPedestalWidth(1)+
			  channelData.tsPedestalWidth(2)*channelData.tsPedestalWidth(2)+
			  channelData.tsPedestalWidth(3)*channelData.tsPedestalWidth(3) );

  if (channelData.hasTimeInfo()) pedConstraint_+=darkCurrent*darkCurrent;

  fcByPe_ = channelData.fcByPE();

  std::vector<float> reconstructedVals;
  SampleVector charges;
  
  double tsTOT = 0, tstrig = 0; // in fC

  for(unsigned int iTS=0; iTS<TSSize_; ++iTS){
    float charge = channelData.tsRawCharge(iTS);
    float ped = channelData.tsPedestal(iTS);
    
    charges.coeffRef(iTS) = charge - ped;
    
    //ADC granularity
    float noiseADC = (1./sqrt(12))*channelData.tsDFcPerADC(iTS);

    //Dark current (for SiPMs)
    float noiseDC=0;
    if((channelData.hasTimeInfo()) && (charge-ped)>channelData.tsPedestalWidth(iTS)) {
      noiseDC = darkCurrent;
    }

    //Electronic pedestal
    float pedWidth = channelData.tsPedestalWidth(iTS);

    //Total uncertainty from all sources
    noiseTerms_.coeffRef(iTS) = noiseADC*noiseADC + pedWidth*pedWidth + noiseDC*noiseDC;

    tsTOT += charge - ped;
    if( iTS==TSOffset_) {
      tstrig += (charge - ped);
    }

  }

  bool status =false;
  if(tstrig >= 0) {
    status = DoFit(charges,reconstructedVals); 
  }
  
  if (!status) {
    reconstructedVals.clear();
    reconstructedVals.push_back(0.);
    reconstructedVals.push_back(888.);
  }

  reconstructedEnergy = reconstructedVals[0]*channelData.tsGain(0);
  chi2 = reconstructedVals[1];

}

bool DoMahiAlgo::DoFit(SampleVector amplitudes, std::vector<float> &correctedOutput) {

  bxs_.resize(BXSize_);
  for (unsigned int iBX=0; iBX<BXSize_; iBX++) {
    bxs_.coeffRef(iBX) = activeBXs_[iBX];
  }

  amplitudes_ = amplitudes;

  pulseMat_.resize(Eigen::NoChange,BXSize_);
  ampVec_ = PulseVector::Zero(BXSize_);
  errVec_ = PulseVector::Zero(BXSize_);

  nP_=0;

  bool status = true;

  int offset=0;
  for (unsigned int iBX=0; iBX<BXSize_; iBX++) {
    
    offset=bxs_.coeff(iBX);

    pulseShapeArray_[iBX] = FullSampleVector::Zero(MaxFSVSize);
    pulseCovArray_[iBX]   = FullSampleMatrix::Constant(0);

    if (offset==0) {
      fullPulseCorr_ = 1.0/double(pulseShapeArray_[iBX].coeff(FullTSofInterest_));
      ampVec_.coeffRef(iBX)= amplitudes_.coeff(TSOffset_);
    }
    else {
      ampVec_.coeffRef(iBX)=0;
    }

    status = UpdatePulseShape(amplitudes_.coeff(TSOffset_ + iBX), 
			      pulseShapeArray_[iBX], 
			      pulseCovArray_[iBX]);


    pulseMat_.col(iBX) = pulseShapeArray_[iBX].segment(FullTSOffset_ + iBX, TSSize_);
  }

  chiSq_ = 9999;
  
  aTaMat_.resize(BXSize_, BXSize_);
  aTbVec_.resize(BXSize_);
  wVec_.resize(BXSize_);

  status = Minimize(); 
  ampVecMin_ = ampVec_;
  bxsMin_ = bxs_;

  if (!status) return status;

  bool foundintime = false;
  unsigned int ipulseintime = 4;

  for (unsigned int iBX=0; iBX<BXSize_; ++iBX) {
    if (bxs_.coeff(iBX)==0) {
      ipulseintime = iBX;
      foundintime = true;
    }
  }
  if (!foundintime) return status;

  correctedOutput.clear();
  correctedOutput.push_back(ampVec_.coeff(ipulseintime)*fullPulseCorr_); //charge
  correctedOutput.push_back(chiSq_); //chi2
  
  return status;
}

bool DoMahiAlgo::Minimize() {
  int iter = 0;
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
    
    float newChiSq=CalculateChiSq();
    float deltaChiSq = newChiSq - chiSq_;
    
    chiSq_ = newChiSq;
    
    if (std::abs(deltaChiSq)<1e-3) break;

    iter++;
    
  }
  
  return status;
}

bool DoMahiAlgo::UpdatePulseShape(double itQ, FullSampleVector &pulseShape, FullSampleMatrix &pulseCov) {
  float t0=meanTime_;
  if (applyTimeSlew_) 
    t0+=HcalTimeSlew::delay(std::max(1.0, itQ), slewFlavor_);

  pulseN_.fill(0);
  pulseM_.fill(0);
  pulseP_.fill(0);

  const double xx[4]={t0, 1.0, 0.0, 3};
  (*pfunctor_)(&xx[0]);

  psfPtr_->getPulseShape(pulseN_);

  const double xxm[4]={t0-dt_, 1.0, 0.0, 3};
  const double xxp[4]={t0+dt_, 1.0, 0.0, 3};

  (*pfunctor_)(&xxm[0]);
  psfPtr_->getPulseShape(pulseM_);

  (*pfunctor_)(&xxp[0]);
  psfPtr_->getPulseShape(pulseP_);

  for (unsigned int iTS=FullTSOffset_; iTS<FullTSOffset_ + TSSize_; iTS++) {
    pulseShape.coeffRef(iTS) = pulseN_[iTS-FullTSOffset_];
    pulseM_[iTS-FullTSOffset_]-=pulseN_[iTS-FullTSOffset_];
    pulseP_[iTS-FullTSOffset_]-=pulseN_[iTS-FullTSOffset_];
  }

  for (unsigned int iTS=FullTSOffset_; iTS<FullTSOffset_+TSSize_; iTS++) {
    for (unsigned int jTS=FullTSOffset_; jTS<iTS+1; jTS++) {

      float tmp=0.5*(pulseP_[iTS-FullTSOffset_]*pulseP_[jTS-FullTSOffset_]
		     + pulseM_[iTS-FullTSOffset_]*pulseM_[jTS-FullTSOffset_]);
		     
      pulseCov(iTS,jTS) += tmp;
      pulseCov(jTS,iTS) += tmp;
      
    }
  }

  return true;  
}



bool DoMahiAlgo::UpdateCov() {
  bool status=true;

  invCovMat_ = noiseTerms_.asDiagonal();
  invCovMat_ += pedConstraint_*SampleMatrix::Ones();

  for (unsigned int iBX=0; iBX<BXSize_; iBX++) {
    if (ampVec_.coeff(iBX)==0) continue;

    unsigned int offset=bxs_.coeff(iBX);

    invCovMat_ += ampVec_.coeff(iBX)*ampVec_.coeff(iBX)
      *pulseCovArray_.at(offset+BXOffset_).block(FullTSOffset_-offset, FullTSOffset_-offset, TSSize_, TSSize_);
    invCovMat_ += ampVec_.coeff(iBX)*fcByPe_*SampleMatrix::Ones();

  }
  
  covDecomp_.compute(invCovMat_);
  
  return status;
}

bool DoMahiAlgo::NNLS() {
  const unsigned int npulse = BXSize_;
  
  for (unsigned int iBX=0; iBX<npulse; iBX++) {
    int offset=bxs_.coeff(iBX);
    pulseMat_.col(iBX) = pulseShapeArray_.at(offset+BXOffset_).segment(FullTSOffset_-offset, TSSize_);
  }

  invcovp_ = covDecomp_.matrixL().solve(pulseMat_);
  aTaMat_.noalias() = invcovp_.transpose().lazyProduct(invcovp_);
  aTbVec_.noalias() = invcovp_.transpose().lazyProduct(covDecomp_.matrixL().solve(amplitudes_));
  
  int iter = 0;
  Index idxwmax = 0;
  double wmax = 0.0;
  double threshold = 1e-11;

  while (true) {    
    if (iter>0 || nP_==0) {
      if ( nP_==npulse ) break;
      
      const unsigned int nActive = npulse - nP_;
      updateWork_ = aTbVec_ - aTaMat_*ampVec_;

      Index idxwmaxprev = idxwmax;
      double wmaxprev = wmax;
      wmax = updateWork_.tail(nActive).maxCoeff(&idxwmax);
      
      if (wmax<threshold || (idxwmax==idxwmaxprev && wmax==wmaxprev)) {
	break;
      }
      
      if (iter>=nMaxItersNNLS_) {
	break;
      }

      //unconstrain parameter
      Index idxp = nP_ + idxwmax;

      aTaMat_.col(nP_).swap(aTaMat_.col(idxp));
      aTaMat_.row(nP_).swap(aTaMat_.row(idxp));
      pulseMat_.col(nP_).swap(pulseMat_.col(idxp));
      std::swap(aTbVec_.coeffRef(nP_),aTbVec_.coeffRef(idxp));
      std::swap(ampVec_.coeffRef(nP_),ampVec_.coeffRef(idxp));
      std::swap(bxs_.coeffRef(nP_),bxs_.coeffRef(idxp));

      wVec_.tail(nActive) = updateWork_.tail(nActive); 
      ++nP_;

    }

    while (true) {
      if (nP_==0) break;     

      ampvecpermtest_ = ampVec_;

      eigen_solve_submatrix(aTaMat_,aTbVec_,ampvecpermtest_,nP_);

      //check solution
      bool positive = true;
      for (unsigned int i = 0; i < nP_; ++i)
        positive &= (ampvecpermtest_(i) > 0);
      if (positive) {
        ampVec_.head(nP_) = ampvecpermtest_.head(nP_);
        break;
      } 

      //update parameter vector
      Index minratioidx=0;

      // no realizable optimization here (because it autovectorizes!)
      double minratio = std::numeric_limits<double>::max();
      for (unsigned int ipulse=0; ipulse<nP_; ++ipulse) {
	if (ampvecpermtest_.coeff(ipulse)<=0.) {
	  const double c_ampvec = ampVec_.coeff(ipulse);
	  const double ratio = c_ampvec/(c_ampvec-ampvecpermtest_.coeff(ipulse));
	  if (ratio<minratio) {
	    minratio = ratio;
	    minratioidx = ipulse;
	  }
	}
      }
      ampVec_.head(nP_) += minratio*(ampvecpermtest_.head(nP_) - ampVec_.head(nP_));
      
      //avoid numerical problems with later ==0. check
      ampVec_.coeffRef(minratioidx) = 0.;
      
      aTaMat_.col(nP_-1).swap(aTaMat_.col(minratioidx));
      aTaMat_.row(nP_-1).swap(aTaMat_.row(minratioidx));
      pulseMat_.col(nP_-1).swap(pulseMat_.col(minratioidx));
      std::swap(aTbVec_.coeffRef(nP_-1),aTbVec_.coeffRef(minratioidx));
      std::swap(ampVec_.coeffRef(nP_-1),ampVec_.coeffRef(minratioidx));
      std::swap(bxs_.coeffRef(nP_-1),bxs_.coeffRef(minratioidx));
      --nP_;

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


float DoMahiAlgo::CalculateChiSq() {
  return covDecomp_.matrixL().solve(pulseMat_*ampVec_ - amplitudes_).squaredNorm();

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
  psfPtr_.reset(new FitterFuncs::PulseShapeFunctor(ps,false,false,false,false,1,0,2.5,0,0.00065,1,TSSize_));
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
