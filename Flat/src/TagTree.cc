#include "../interface/TagTree.h"

TagTree::TagTree() {
//STARTCUSTOMCONST
  
  for (auto ibeta : ibetas) {
    for (auto N : Ns) {
      for (auto order : orders) {
        ECFParams p;
        p.ibeta = ibeta;
        p.N = N;
        p.order = order;
        ecfParams.push_back(p);
        fj1ECFNs[p] = -1;
      }
    }
  }


//ENDCUSTOMCONST
}

TagTree::~TagTree() {
//STARTCUSTOMDEST
//ENDCUSTOMDEST
}

void TagTree::Reset() {
//STARTCUSTOMRESET

  for (auto p : ecfParams) { 
    fj1ECFNs[p] = -1;
  }

//ENDCUSTOMRESET
    runNumber = 0;
    lumiNumber = 0;
    eventNumber = 0;
    npv = 0;
    pu = 0;
    mcWeight = -1;
    filter_maxRecoil = -1;
    filter_whichRecoil = -1;
    sf_ewkV = 1;
    sf_qcdV = 1;
    sf_ewkV2j = 1;
    sf_qcdV2j = 1;
    sf_qcdTT = 1;
    sf_pu = 1;
    sf_npv = 1;
    sf_tt = 1;
    sf_phoPurity = 1;
    pfmet = -1;
    puppimet = -1;
    partonPt = -1;
    partonEta = -1;
    partonIsHad = 0;
    partonSize = -1;
    partonIsReco = 0;
    partonPdgId = 0;
    nFatjet = 0;
    fj1Tau32 = -1;
    fj1Tau21 = -1;
    fj1Tau32SD = -1;
    fj1Tau21SD = -1;
    fj1MSD = -1;
    fj1Pt = -1;
    fj1Phi = -1;
    fj1Eta = -1;
    fj1M = -1;
    fj1MaxCSV = -1;
    fj1SubMaxCSV = -1;
    fj1MinCSV = -1;
    fj1DoubleCSV = -1;
    fj1IsTight = 0;
    fj1IsLoose = 0;
    fj1RawPt = -1;
    fj1NHF = 0;
    fj1HTTMass = -1;
    fj1HTTFRec = -1;
    fj1IsClean = 0;
}

void TagTree::WriteTree(TTree *t) {
  treePtr = t;
//STARTCUSTOMWRITE

  for (auto p : ecfParams) { 
    TString ecfn(makeECFString(p));
    Book("fj1"+ecfn,&(fj1ECFNs[p]),"fj1"+ecfn+"/F");
  }

//ENDCUSTOMWRITE
    Book("runNumber",&runNumber,"runNumber/I");
    Book("lumiNumber",&lumiNumber,"lumiNumber/I");
    Book("eventNumber",&eventNumber,"eventNumber/l");
    Book("npv",&npv,"npv/I");
    Book("pu",&pu,"pu/I");
    Book("mcWeight",&mcWeight,"mcWeight/F");
    Book("filter_maxRecoil",&filter_maxRecoil,"filter_maxRecoil/F");
    Book("filter_whichRecoil",&filter_whichRecoil,"filter_whichRecoil/F");
    Book("sf_ewkV",&sf_ewkV,"sf_ewkV/F");
    Book("sf_qcdV",&sf_qcdV,"sf_qcdV/F");
    Book("sf_ewkV2j",&sf_ewkV2j,"sf_ewkV2j/F");
    Book("sf_qcdV2j",&sf_qcdV2j,"sf_qcdV2j/F");
    Book("sf_qcdTT",&sf_qcdTT,"sf_qcdTT/F");
    Book("sf_pu",&sf_pu,"sf_pu/F");
    Book("sf_npv",&sf_npv,"sf_npv/F");
    Book("sf_tt",&sf_tt,"sf_tt/F");
    Book("sf_phoPurity",&sf_phoPurity,"sf_phoPurity/F");
    Book("pfmet",&pfmet,"pfmet/F");
    Book("puppimet",&puppimet,"puppimet/F");
    Book("partonPt",&partonPt,"partonPt/F");
    Book("partonEta",&partonEta,"partonEta/F");
    Book("partonIsHad",&partonIsHad,"partonIsHad/I");
    Book("partonSize",&partonSize,"partonSize/F");
    Book("partonIsReco",&partonIsReco,"partonIsReco/I");
    Book("partonPdgId",&partonPdgId,"partonPdgId/I");
    Book("nFatjet",&nFatjet,"nFatjet/I");
    Book("fj1Tau32",&fj1Tau32,"fj1Tau32/F");
    Book("fj1Tau21",&fj1Tau21,"fj1Tau21/F");
    Book("fj1Tau32SD",&fj1Tau32SD,"fj1Tau32SD/F");
    Book("fj1Tau21SD",&fj1Tau21SD,"fj1Tau21SD/F");
    Book("fj1MSD",&fj1MSD,"fj1MSD/F");
    Book("fj1Pt",&fj1Pt,"fj1Pt/F");
    Book("fj1Phi",&fj1Phi,"fj1Phi/F");
    Book("fj1Eta",&fj1Eta,"fj1Eta/F");
    Book("fj1M",&fj1M,"fj1M/F");
    Book("fj1MaxCSV",&fj1MaxCSV,"fj1MaxCSV/F");
    Book("fj1SubMaxCSV",&fj1SubMaxCSV,"fj1SubMaxCSV/F");
    Book("fj1MinCSV",&fj1MinCSV,"fj1MinCSV/F");
    Book("fj1DoubleCSV",&fj1DoubleCSV,"fj1DoubleCSV/F");
    Book("fj1IsTight",&fj1IsTight,"fj1IsTight/I");
    Book("fj1IsLoose",&fj1IsLoose,"fj1IsLoose/I");
    Book("fj1RawPt",&fj1RawPt,"fj1RawPt/F");
    Book("fj1NHF",&fj1NHF,"fj1NHF/I");
    Book("fj1HTTMass",&fj1HTTMass,"fj1HTTMass/F");
    Book("fj1HTTFRec",&fj1HTTFRec,"fj1HTTFRec/F");
    Book("fj1IsClean",&fj1IsClean,"fj1IsClean/I");
}

