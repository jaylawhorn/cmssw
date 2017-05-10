#include "RecoLocalCalo/HcalRecAlgos/interface/DoMahiAlgo.h"
#include <iostream>
#include <fstream> 

//double pulse_temp[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};// TEST
//double digi_temp[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};// TEST
//double noise_temp[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};// TEST

void eigen_solve_submatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned NP);

void DoMahiAlgo::configurePulseShapes(NewPulseShapes pulseShapeInfo) {

  pulseShapeInfo_ = pulseShapeInfo;

}

void DoMahiAlgo::getPulseShape(float q, HcalDetId detID, float t, SampleVector &pulseShape, float sigma=0) {

  int offsetTS=0;
  float newT=t;
  if (abs(t)<12.5) {
    offsetTS=0;
    newT=t;
  }
  else if (abs(t-25)<12.5) {
    offsetTS=-1;
    newT=t-25;
  }
  else if (abs(t+25)<12.5) {
    offsetTS=1;
    newT=t+25;
  }

  float sum=pulseShapeInfo_.getPulseFracNorm(q,newT);

  for (int i=0; i<10; i++) {
    if (q<0 || i+offsetTS>9 || i+offsetTS<0) pulseShape[i]=0;
    else pulseShape[i] = pulseShapeInfo_.getPulseFrac(q,newT,i+offsetTS)/sum;
  }

}

void DoMahiAlgo::Apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalDetId & detID, const HcalCalibrations & calibs, std::vector<float> & correctedOutput) {

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

  std::vector<float> fitParsVec;

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

void DoMahiAlgo::phase1Apply(const HBHEChannelInfo& channelData,
			     std::vector<float>& reconstructedVals) {

  const unsigned cssize = channelData.nSamples();

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

  bool status =false;
  if(tstrig >= 0) {
    status = DoFit(charges, gains, reconstructedVals); 
  }

  if (!status) {
    reconstructedVals.clear();
    reconstructedVals.push_back(0.); //energy IT
    reconstructedVals.push_back(0.); // energy next
    reconstructedVals.push_back(0.); // energy previous
    reconstructedVals.push_back(-9999); // chi2
  }

  /*
  if(tsTOTen>20) std::cout << " --> ID=" << channelData.id() << " depth=" << channelData.id().depth() << " eta=" << channelData.id().ieta() << " phi=" << channelData.id().iphi() << std::endl;

  for(unsigned int ip=0; ip<cssize; ++ip){
    if( ip >= (unsigned) HcalConst::maxSamples ) continue; // Too many samples than what we wanna fit (10 is enough...) -> skip them                                                                                                   
    //      double chi2TS=(digi_temp[ip]-pulse_temp[ip])*(digi_temp[ip]-pulse_temp[ip])/(noise_temp[ip]*noise_temp[ip]);                                                                                                               
    //    if(tsTOTen>20) std::cout << "TS=" << ip << " == fittedPulsep[ip](GeV)=" << pulse_temp[ip] << " digi[ip](GeV)=" << digi_temp[ip] << std::endl;

  }
  */

}


bool DoMahiAlgo::DoFit(SampleVector amplitudes, SampleVector gains, std::vector<float> &correctedOutput) {

  bool verbose=false;

  _nP = 0;
  //ECAL does it better -- to be fixed
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#L151-L171
  _bxs.resize(3);
  _bxs << -1,0,1;
  _nPulseTot = _bxs.rows();

  //_detID = HcalDetId(detID.rawId());
    
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

  pulseShape = PulseVector::Zero(10);
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

  if(verbose) {
    std::cout << "------" << std::endl;
    std::cout << "input: " ;
    for (int i=0; i<10; i++) {
      std::cout << _amplitudes.coeff(i) << ", ";
    }
    std::cout << std::endl;

    std::cout << "output: ";
    std::cout << _ampVec.coeff(ipulseprevtime) << ", " << _ampVec.coeff(ipulseintime) << ", " << _ampVec.coeff(ipulsenexttime) << std::endl;
    std::cout << "output we care about: " << std::endl;
   
  }

  std::vector<double> ans;

  getPulseShape(_ampVec.coeff(ipulseprevtime), _detID, -25, pulseShape,0); 
  if(verbose) std::cout << "prev: ";
  for (int i=0; i<10; i++) {
    ans.push_back(_ampVec.coeff(ipulseprevtime)*pulseShape.coeff(i));
    if(verbose) std::cout << _ampVec.coeff(ipulseprevtime)*pulseShape.coeff(i) << ", ";
  }
  if(verbose) std::cout << std::endl;

  if(verbose) std::cout << "intime: ";
  getPulseShape(_ampVec.coeff(ipulseintime), _detID, 0, pulseShape,0);
  for (int i=0; i<10; i++) {
    ans[i]+= _ampVec.coeff(ipulseintime)*pulseShape.coeff(i);
    if(verbose) std::cout << _ampVec.coeff(ipulseintime)*pulseShape.coeff(i) << ", ";
  }
  if(verbose) std::cout << std::endl;

  if(verbose) std::cout << "next: ";
  getPulseShape(_ampVec.coeff(ipulsenexttime), _detID, 25, pulseShape,0);
  for (int i=0; i<10; i++) {
    ans[i]+= _ampVec.coeff(ipulsenexttime)*pulseShape.coeff(i);
    if(verbose) std::cout << _ampVec.coeff(ipulsenexttime)*pulseShape.coeff(i) << ", ";
  }

  if(verbose) {
    std::cout << std::endl;
    std::cout << "SUMMED ANSWER: " << std::endl;
    for (int i=0; i<10; i++) {
      std::cout << ans[i] << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "chi2: " ;
    std::cout << _chiSq << std::endl;
  }
  
  double gain=gains.coeff(4); // this is the same for each TS

  correctedOutput.clear();
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
    // this should be set with the new pulse shape
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
      
      //if (iter>=500) {
      if (iter>=50) {
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

  ////////
  // THIS IS THE ELECTONIC noise

  const double pederr2 = 1;
  _invCovMat = pederr2*SampleMatrix::Constant(1);

  ////////
  // THIS IS THE electronic noise + ADC granularity + photoStatistics

  for (int i=0; i<10; i++) {
    double ifC=_amplitudes.coeff(i);
    double sigma = 0;
    if(ifC < 75) sigma = (0.577 + 0.0686*ifC)/3.; 
    else sigma = (2.75  + 0.0373*ifC + 3e-6*ifC*ifC)/3.; 

    double sigma2 = ifC/sqrt(ifC/0.3305);
    if (ifC<1) sigma2=1/sqrt(1/0.3305);
    _invCovMat(i, i) += (1 + sigma*sigma + sigma2*sigma2);

  }

  ///////
  // pull all together now

  for (int k=0; k<_ampVec.size(); k++) {
    double ifC=_ampVec.coeff(k);
    if (ifC==0) continue;
    double nomT = _bxs.coeff(k)*25;
    int maxTS = 4+_bxs.coeff(k);

    // note: this should be set with the newPulseShape
    getPulseShape(ifC, _detID, nomT, pulseShape,0); 
    getPulseShape(ifC, _detID, nomT+deltaT, pulseShapeP,0); 
    getPulseShape(ifC, _detID, nomT-deltaT, pulseShapeM,0); 

    for (int xx=0; xx<10; xx++) {
      if (pulseShape.coeff(maxTS)==0) pulseShape.coeffRef(xx)=0;
      else pulseShape.coeffRef(xx) =  pulseShape.coeff(xx)/pulseShape.coeff(maxTS);

      if (pulseShapeP.coeff(maxTS)==0) pulseShape.coeffRef(xx)=0;
      else pulseShapeP.coeffRef(xx) = pulseShapeP.coeff(xx)/pulseShapeP.coeff(maxTS);

      if (pulseShapeM.coeff(maxTS)==0) pulseShape.coeffRef(xx)=0;
      else pulseShapeM.coeffRef(xx) = pulseShapeM.coeff(xx)/pulseShapeM.coeff(maxTS);
    }

    for (int i=0; i<10; i++) {
      for (int j=0; j<i+1; j++) {
	double tmp=0.5*((pulseShapeP.coeff(i)-pulseShape.coeff(i))*(pulseShapeP.coeff(j)-pulseShape.coeff(j)) 
			+ (pulseShapeM.coeff(i)-pulseShape.coeff(i))*(pulseShapeM.coeff(j)-pulseShape.coeff(j)));
	//std::cout << i << ", " << j << ", " << tmp << std::endl;
	_invCovMat(i,j) += ifC*ifC*tmp;
	_invCovMat(j,i) += ifC*ifC*tmp;
	
      }
    }
    //if (nomT==0)std::cout << _invCovMat << std::endl;

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
  //  pulseShapeObj = PulseShapes();
  deltaT = 2.5;
}
