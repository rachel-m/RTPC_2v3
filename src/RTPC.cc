// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class RTPC
// Geometry and materials of the RTPC
// 03/05/14 JRMA

#include "RTPC.hh"
#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_SpinEqRhs.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4NystromRK4.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
//#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "TOSCAField3D.hh"
#include "TOSCAField2D.hh"
#include "BonusBField.hh"
#include "Field3D.hh"
#include "ArraySD.hh"
//#include "TH3D.h"
enum {EUniField, ENonUniField, EToscaField3D, EBonusBField, EToscaField2D};

//----------------------------------------------------------------------------
RTPC::RTPC(DetectorConstruction* rectagg)
{
  fRtag = rectagg;
  fMaw = NULL;
  //  fPWT=NULL;
  fIsInteractive=1;
  fIsSrcPb = false;
  fIsOverlapVol = false;
  fNtarget = 1;
  fNgas = 2;
  fXmin = fBmin = fBmax = fBz = 0.0;
  fTx = fTy = fTz = fRst = fTst = fZst = fRHe1 = fRHe2 = fRw = fRbl = fTbl =
    fSbl = fZbl =
    fRbaf = fTbaf = fRG = fTG = fTend1 = fTend2 = fTsh = fZsh = 0;
  fNw = fNbl = 0;
  fORin = fORout = fOZ = fOTh = 0;
  fFieldMap = NULL;
  fBFieldType = EUniField;
  fMst = 0; // 0 = kapton, 1 = Al
  fOvHe = 0.0;
  fWMat = 0;
}

//----------------------------------------------------------------------------
RTPC::~RTPC()
{
 
}

//---------------------------------------------------------------------------
G4int RTPC::ReadParameters(G4String file)
{
  //
  // Initialise magnet setup
  //
  DefaultInit();
  char* keylist[] = { (char *)"Bfield-Box:", (char *)"Target-Dim:", 
		      (char *)"Uni-Field:", (char *)"Gas-Volume:", 
		      (char *)"End-Foil:", (char *)"Target-Mat:",
		      (char *)"Gas-Mat:", (char *)"BL-Dim:",
		      (char *)"GEM-Dim:",(char *)"Non-Uni-Field:",
		      (char *)"Tosca-Field:",(char *)"Baffle-Dim:",
		      (char *)"GEM-Pixels:", (char*)"Outer-Chamber:",
		      (char*)"Be-Window:", (char*)"Beam-Colli:",
		      (char*)"Bonus-Field:", (char*)"Tosca-2DField:",NULL };
  enum { EBfieldBox, ETargDim, EMagField, EGasVol, EEndFoil, ETargMat, EGasMat,
	 EBLDim, EGEMDim, ENonUniField, EToscaField, EBaffleDim, EGEMPixels,
	 EOutCham, EBeWind, EBeamColli, EBonusField, ETosca2Field,ENULL };
  G4int ikey, iread, ierr;
  ierr = 0;
  char line[256];
  char delim[64];
  FILE* pdata;
  if( (pdata = fopen(file.data(),"r")) == NULL ){
    printf("Error opening source parameter file: %s\n",file.data());
    return -1;
  }
  while( fgets(line,256,pdata) ){
    if( line[0] == '#' ) continue;       // comment #
    sscanf(line,"%s",delim);
    for(ikey=0; ikey<ENULL; ikey++)
      if(!strcmp(keylist[ikey],delim)) break;
    switch( ikey ){
    default:
      printf("Unknown setup key\n");
      return -1;
      break;
    case EBfieldBox:
      // 1/2 dimensions of magnetic field region
      iread = sscanf(line,"%*s%lf%lf%lf",&fTx,&fTy,&fTz);
      if( iread != 3 ) ierr++;
      break;
    case ETargDim:
      // dimensions target straw
      iread = sscanf(line,"%*s%lf%lf%lf%lf",&fRst,&fZst,&fTst,&fMst);
      if( iread < 3 ) ierr++;
      break;
    case EMagField:
      // Z component of uniform magnetic field
      iread = sscanf(line,"%*s%lf",&fBz);
      if( iread != 1 ) ierr++;
      fBFieldType = EUniField;
      break;
    case EGasVol:
      // Radii of active gas region
      iread =
	sscanf(line,"%*s%lf%lf%d%lf%d%lf",
	       &fRHe1,&fRHe2,&fIsSep,&fRw,&fNw,&fOvHe);
      if( iread < 5 ) ierr++;
      break;
    case EOutCham:
      // Outer chamber
      iread =
	sscanf(line,"%*s%lf%lf%lf%lf",&fORin,&fORout,&fOZ,&fOTh);
      if( iread != 4 ) ierr++;
      break;
    case EBeWind:
      // Be Window
      iread =
	sscanf(line,"%*s%lf%d",&fWTh,&fWMat);
      if( iread < 1 ) ierr++;
      break;
    case EBeamColli:
      // Beam collimators
      iread =
	sscanf(line,"%*s%lf%lf%lf%lf",&fBmClen1,&fBmCr1,&fBmClen2,&fBmCr2);
      if( iread != 4 ) ierr++;
      break;
    case EEndFoil:
      // End cover
      iread = sscanf(line,"%*s%lf%lf",&fTend1,&fTend2);
      if( iread != 2 ) ierr++;
      break;
    case ETargMat:
      // Target material and density
      iread = sscanf(line,"%*s%d%lf",&fNtarget,&fTDens);
      if( iread != 2 ) ierr++;
      break;
    case EGasMat:
      // He density
      iread = sscanf(line,"%*s%d%lf",&fNgas,&fHeDens);
      if( iread != 2 ) ierr++;
      break;
    case EBLDim:
      // Beam line radii
      //iread = sscanf(line,"%*s%lf%lf%d%lf%lf",&fRblI,&fRblO,&fNbl,&fTsh,&fZsh);
      iread = sscanf(line,"%*s%lf%lf%d%lf%lf",&fRbl,&fTbl,&fNbl,&fSbl,&fZbl);
      if( iread != 5) ierr++;
      break;
    case EGEMDim:
      // Be interior
      iread = sscanf(line,"%*s%lf%lf",&fRG,&fTG);
      if( iread != 2 ) ierr++;
      break;
    case ENonUniField:
      // Non uniform magnetic field
      iread = sscanf(line,"%*s%lf%lf%lf",&fXmin,&fBmin,&fBmax);
      if( iread != 3 ) ierr++;
      fBFieldType = ENonUniField;
      break;
    case EToscaField:
      // Tosca magnetic field
      fFieldMap = new char[64];
      iread = sscanf(line,"%*s%lf%lf%lf%s",
		     &fTXoff,&fTYoff,&fTZoff,fFieldMap);
      if( iread != 4 ) ierr++;
      fBFieldType = EToscaField3D;
      break;
    case ETosca2Field:
      // Tosca magnetic field
      fFieldMap = new char[64];
      iread = sscanf(line,"%*s%lf%lf%lf%s",
		     &fTXoff,&fTYoff,&fTZoff,fFieldMap);
      if( iread != 4 ) ierr++;
      fBFieldType = EToscaField2D;
      break;
    case EBonusField:
      // Bonus magnetic field
      fFieldMap = new char[64];
      iread = sscanf(line,"%*s%lf%lf%lf%s%lf",
		     &fTXoff,&fTYoff,&fTZoff,fFieldMap,&fBScaleFac);
      if( iread != 5 ) ierr++;
      fBFieldType = EBonusBField;
      break;
    case EBaffleDim:
      // Baffle dimensions
      iread = sscanf(line,"%*s%lf%lf",&fRbaf,&fTbaf);
      if( iread != 2 ) ierr++;
      break;
    case EGEMPixels:
      // no. pixels in z,phi for GEM pads
      iread = sscanf(line,"%*s%d%d",&fZPixG,&fPhiPixG);
      if( iread != 2 ) ierr++;
      break;
    }
  }
  if( ierr ) printf("Error detected in procedure RTPC::ReadParameters: %s\n",
		    line);
  return ierr;
}

//----------------------------------------------------------------------------
G4VPhysicalVolume* RTPC::Construct(G4LogicalVolume* maw){
  // Main tank build
  // Then build view ports and target
  //
  fMaw = maw;             // save the mother volume
  //
  G4double Rtot = fRHe2+fRG;
  G4double ZHe = fZst + fOvHe;
  //G4double Zbl = 0.5*(fTz - fZst);
  G4double Zbl = 0.5*(fTz - fOZ);
  G4double ZblOff = -fTz + Zbl;
  //G4double ZblOff = -fOZ + Zbl;
  fZbl = Zbl/fNbl - fTbl;
  G4double RendI;
  //if((fRbl+fTbl) > fRst) RendI = fRbl + fTbl;
  //else RendI = fRst;
  RendI = fRst;
  //G4double Rbld = fRblO - fRblI;     //diff beam-line radii
  // 
  G4Box* Bfield = new G4Box("Bfield", fTx,fTy,fTz);
  G4Tubs* CSpec = new G4Tubs("ExtCont",0.0,fORout,fOZ,0.0,360*deg);         
  G4Tubs* CSpecI = new G4Tubs("ExtContI",0.0,fORin,fOZ-2*fOTh,0.0,360*deg);
  G4Tubs* CSpecW = new G4Tubs("ExtContW",0.0,fRbl,fOTh,0.0,360*deg);
  G4Tubs* Window = new G4Tubs("Window",0.0,fRbl,fWTh,0.0,360*deg);
  //G4Tubs* RSpec = new G4Tubs("RecSpec",0.0,Rtot,fZst+2*fOvHe,0.0,360*deg);
  //G4int nbl = 2*(fNbl+1);
  G4int nbl = 6*fNbl;
  G4Tubs** Bl = new G4Tubs*[nbl];
  G4LogicalVolume** LBl = new G4LogicalVolume*[nbl];
  char name[64];
  G4Tubs* BlinI = new G4Tubs("BLinI",0.0,fRbl,Zbl,0.0,360*deg);
  G4Tubs* BlinO = new G4Tubs("BLinO",0.0,fRbl+fTbl,Zbl,0.0,360*deg);
  G4Tubs* Colli1 = new G4Tubs("Colli1",fBmCr1,fRbl,fBmClen1,0.0,360*deg);
  G4Tubs* Colli2 = new G4Tubs("Colli2",fBmCr2,fRbl,fBmClen2,0.0,360*deg);
  G4double rbl = fRbl;
  for(G4int ibl=0; ibl<nbl; ibl+=3){
    sprintf(name,"BlO-%d",ibl);
    Bl[ibl] = new G4Tubs(name,0.0,rbl+fTbl,fZbl,0.0,360*deg);
    sprintf(name,"BlI-%d",ibl);
    Bl[ibl+1] = new G4Tubs(name,0.0,rbl,fZbl,0.0,360*deg);
    sprintf(name,"BlE-%d",ibl);
    Bl[ibl+2] = new G4Tubs(name,rbl,rbl+fTbl+fSbl,fTbl,0.0,360*deg);
    rbl += fSbl;
  }
  G4Tubs* Shield;
  if(fTsh)
    Shield = new G4Tubs("Shield",fRbl+fTbl,Rtot,fTsh,0.0,360*deg);
  G4Tubs* Straw = new G4Tubs("Straw",fRst-fTst,fRst,ZHe,0.0,360*deg);
  G4Tubs* Straw1 = new G4Tubs("Straw1",0.0,fRst-fTst,fZst,0.0,360*deg);
  G4Tubs* WindT = new G4Tubs("WindT",0.0,fRst-fTst,fWTh,0.0,360*deg);
  G4Tubs* He1;
  G4Tubs* He2;
  if(fIsSep){
    He1 = new G4Tubs("He1",fRst,fRHe1,ZHe,0.0,360*deg);
    He2 = new G4Tubs("He2",fRHe1,fRHe2,ZHe,0.0,360*deg);
  }
  else
    He2 = new G4Tubs("He2",fRst,fRHe2,ZHe,0.0,360*deg);
  G4Tubs* Gem1 = new G4Tubs("Gem1",fRHe2,fRHe2+fTG,ZHe,0.0,360*deg);
  G4Tubs* Gem2 = new G4Tubs("Gem2",fRHe2+fTG,Rtot-fTG,ZHe,0.0,360*deg);
  G4Tubs* Gem3 = new G4Tubs("Gem3",Rtot-fTG,Rtot,ZHe,0.0,360*deg);
  G4Tubs* Wire = new G4Tubs("Wire",0.0,fRw,ZHe,0.0,360*deg);
  G4Tubs* Baffle;
  if(fTbaf) Baffle = new G4Tubs("Baffle",fRbaf-fTbaf,fRbaf,fZst,0.0,360*deg);
  G4Tubs* End1 = new G4Tubs("End",RendI,Rtot,fTend1,0,360*deg);
  G4Tubs* End2 = new G4Tubs("End",RendI,Rtot,fTend2,0,360*deg);
  //
  G4LogicalVolume* LBfield =
    new G4LogicalVolume(Bfield, fRtag->GetAir(),"LBfield",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LBfield,"PBfield",fMaw,false,0,
		    fIsOverlapVol);
  SetMagField(LBfield);

  G4LogicalVolume* LCSpec =
    new G4LogicalVolume(CSpec,fRtag->GetKapton(),"LCSpec",0,0,0);
  G4LogicalVolume* LCSpecI =
    new G4LogicalVolume(CSpecI,fRtag->GetVac(),"LCSpecI",0,0,0);
  //    new G4LogicalVolume(CSpecI,fRtag->GetHe(fHeDens),"LCSpecI",0,0,0);
  G4LogicalVolume* LCSpecW =
    new G4LogicalVolume(CSpecW,fRtag->GetVac(),"LCSpecW",0,0,0);
  //    new G4LogicalVolume(CSpecW,fRtag->GetHe(fHeDens),"LCSpecW",0,0,0);
  G4LogicalVolume* LWindow;
  if(fWMat == 0)
    LWindow = new G4LogicalVolume(Window,fRtag->GetBe(),"LWindow",0,0,0);
  else
    LWindow = new G4LogicalVolume(Window,fRtag->GetAl(),"LWindow",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LCSpec,"PCSpec",LBfield,false,0,
		    fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LCSpecI,"PCSpecI",LCSpec,false,0,
		    fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LWindow,"PWindow",LCSpecW,false,0,
		    fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,fOZ-fOTh),LCSpecW,"PCSpecW1",LCSpec,false,0,
		    fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,-fOZ+fOTh),LCSpecW,"PCSpecW2",LCSpec,false,0,
		    fIsOverlapVol);

  //G4LogicalVolume* LRSpec =
  //  new G4LogicalVolume(RSpec, fRtag->GetVac(),"LRSpec",0,0,0);
  // new G4PVPlacement(0,G4ThreeVector(0,0,0),LRSpec,"PRSpec",LCSpecI,false,0,
  //		    fIsOverlapVol);
  //
  // Create vacuum and Al beam pipe logical volumes
  // Put vacuum inside Al. 
  G4LogicalVolume* LBlinI = 
    new G4LogicalVolume(BlinI,fRtag->GetVac(),"LBlinI",0,0,0);
  G4LogicalVolume* LBlinO =
    new G4LogicalVolume(BlinO,fRtag->GetAl(),"LBlinO",0,0,0);
  G4LogicalVolume* LColli1 =
    new G4LogicalVolume(Colli1,fRtag->GetW(),"LColli1",0,0,0);
  G4LogicalVolume* LColli2 =
    new G4LogicalVolume(Colli2,fRtag->GetW(),"LColli2",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LColli1,"PColli1",LBlinI,false,
		    0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LBlinI,"PBlin",LBlinO,false,
		    0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,-fZst-2*fTend2-fBmClen2),LColli2,
		    "PColli2",LCSpecI,false,0,fIsOverlapVol);
  for(G4int ibl=0; ibl<nbl; ibl+=3){
    sprintf(name,"LBlI-%d",ibl);
    LBl[ibl+1] = new G4LogicalVolume(Bl[ibl+1],fRtag->GetVac(),name,0,0,0);
    sprintf(name,"LBlO-%d",ibl);
    LBl[ibl] = new G4LogicalVolume(Bl[ibl],fRtag->GetAl(),name,0,0,0);
    sprintf(name,"LBlE-%d",ibl+1);
    LBl[ibl+2] = new G4LogicalVolume(Bl[ibl+2],fRtag->GetAl(),name,0,0,0);
    sprintf(name,"PB-%d",ibl);
    new G4PVPlacement(0,G4ThreeVector(0,0,0),LBl[ibl+1],name,LBl[ibl],false,
		    0,fIsOverlapVol);
  }
  // Place beam-line elements
  new G4PVPlacement(0,G4ThreeVector(0,0,ZblOff),LBlinO,"PBlU",LBfield,
		    false,0,fIsOverlapVol);
  //G4double zz = fZst + fZbl;
  G4double zz = fOZ + fZbl;
  for(G4int ibl=0; ibl<nbl; ibl+=3){
    G4LogicalVolume* lbf = LBfield;
    if(ibl >= nbl/2) lbf = maw;
    sprintf(name,"PBl-%d",ibl);
    new G4PVPlacement(0,G4ThreeVector(0,0,zz),LBl[ibl],name,lbf,false,
		      0,fIsOverlapVol);
    zz = zz + fZbl + fTbl;
    sprintf(name,"PBlE-%d",ibl);
    new G4PVPlacement(0,G4ThreeVector(0,0,zz),LBl[ibl+2],name,lbf,false,
		      0,fIsOverlapVol);
    zz = zz + fZbl +fTbl;
  }
  // Downstream shield
  G4LogicalVolume* LShield = NULL;
  if( fTsh ){
    LShield =  new G4LogicalVolume(Shield,fRtag->GetAl(),"Shield",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,-ZblOff+fZsh),LShield,"PShield",
		      LBfield,false,0,fIsOverlapVol);
  }
  G4Material* targWall;
  if( fMst == 0 )
    targWall = fRtag->GetKapton();
  else
    targWall = fRtag->GetAl();
  G4LogicalVolume* LStraw =
    new G4LogicalVolume(Straw, targWall,"LStraw",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,fOvHe),LStraw,"PStraw",LCSpecI,false,0,
		    fIsOverlapVol);
  G4Material* target;
  if(fNtarget == 1) target = fRtag->GetCoH2(fTDens);
  else target = fRtag->GetCoD2(fTDens);
  G4cout << target << endl;
  G4LogicalVolume* LStraw1 =
    new G4LogicalVolume(Straw1, target, "LStraw1 ",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LStraw1,"PStraw1",LCSpecI,false,3,fIsOverlapVol);
  G4LogicalVolume* LWindT =
	new G4LogicalVolume(WindT,fRtag->GetBe(),"LWindT",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,-fZst+2*fWTh),
		    LWindT,"PWindT1",LStraw1,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,fZst-2*fWTh),
		    LWindT,"PWindT2",LStraw1,false,0,fIsOverlapVol);
  G4Material* gas;
  if(fNgas == 1) gas = fRtag->GetCoH2(fHeDens);
  else gas = fRtag->GetHe(fHeDens);
  G4cout << gas << endl;
  G4LogicalVolume* LHe1 = NULL;
  G4int he2ID = 0;
  if(fIsSep){
    LHe1 = new G4LogicalVolume(He1, gas,"LHe1",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,fOvHe),
		      LHe1,"PHe1",LCSpecI,false,0,fIsOverlapVol);
    he2ID = 1;
  }
  G4LogicalVolume* LHe2 = 
    new G4LogicalVolume(He2, gas,"LHe2",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,fOvHe),
		    LHe2,"PHe2",LCSpecI,false,he2ID,fIsOverlapVol);
  if(fTbaf){
    G4LogicalVolume* lhe = LHe1;
    if( !lhe ) lhe = LHe2; 
    G4LogicalVolume* LBaf = 
      new G4LogicalVolume(Baffle, fRtag->GetAl(),"LBaf",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,0),
		      LBaf,"PBaf",lhe,false,0,fIsOverlapVol);
  }    
  G4LogicalVolume* LGem1 = 
    new G4LogicalVolume(Gem1, fRtag->GetKapton(),"LGem1",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,fOvHe),
		    LGem1,"PGem1",LCSpecI,false,0,fIsOverlapVol);
  G4LogicalVolume* LGem2 = 
    new G4LogicalVolume(Gem2, fRtag->GetAr(),"LGem2",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,fOvHe),
		    LGem2,"PGem2",LCSpecI,false,2,fIsOverlapVol);
  G4LogicalVolume* LGem3 = 
    new G4LogicalVolume(Gem3, fRtag->GetKapton(),"LGem3",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,fOvHe),
		    LGem3,"PGem3",LCSpecI,false,0,fIsOverlapVol);
  G4LogicalVolume* LEnd1 = 
    new G4LogicalVolume(End1, fRtag->GetKapton(),"LEnd1",0,0,0);
  G4LogicalVolume* LEnd2 = 
    new G4LogicalVolume(End2, fRtag->GetKapton(),"LEnd2",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,fZst+fTend1+2*fOvHe),
		    LEnd1,"PEnd1",LCSpecI,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,-(fZst+fTend2)),
		    LEnd2,"PEnd2",LCSpecI,false,0,fIsOverlapVol);
  // Field Wires....fNw > 0
  G4LogicalVolume* LWire;
  if(fNw > 0){
    LWire = new G4LogicalVolume(Wire, fRtag->GetW(),"LWire",0,0,0);
    G4double r = fRHe1+fRw;
    G4double dth = 360*deg/fNw;
    char wnm[32];
    for(G4int i=0; i<fNw; i++){
      G4double th = dth * (i-1);
      G4double x = r*cos(th);
      G4double y = r*sin(th);
      sprintf(wnm,"Pwire%d",i);
      new G4PVPlacement(0,G4ThreeVector(x,y,0),
			LWire,wnm,LHe2,false,0,fIsOverlapVol);
    }
  }
  LStraw->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  LStraw1->SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));
  if(LHe1) LHe1->SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
  LHe2->SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
  LCSpec->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  LCSpecI->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  LWindow->SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));
  if( fNw )
    // LWire->SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
    LWire->SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
  // }
  // G4VisAttributes * LWireVisAtt = new G4VisAttributes(G4Colour::Green());
  // LWireVisAtt->SetForceSolid(true);
  // LWireVisAtt->SetLineWidth(10);

  // }
  //
  LBfield->SetVisAttributes (G4VisAttributes(G4Colour::Yellow()));
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4double parm[5];
  parm[0] = fZst;
  parm[1] = fRHe1;
  parm[2] = fRHe2;
  parm[3] = (G4double)fZPixG;
  parm[4] = (G4double)fPhiPixG;
  fArrSD = new ArraySD("ArraySD", 9, parm );
  SDman->AddNewDetector(fArrSD);
  if( LHe1 ) LHe1->SetSensitiveDetector( fArrSD );
  LHe2->SetSensitiveDetector( fArrSD );
  LGem2->SetSensitiveDetector( fArrSD );
  LStraw1->SetSensitiveDetector( fArrSD );
  return NULL;
}
//----------------------------------------------------------------------------
void RTPC::SetMagField( G4LogicalVolume* gapLog )
{
  G4FieldManager  *magFieldMgr;
  G4MagneticField* magField;
  G4double pos[3];
  G4double b[3];
  pos[0] = pos[1] = pos[2] = 0.0;
  switch(fBFieldType){
  case EUniField:
  default:
    magField = new G4UniformMagField( G4ThreeVector(0,0,fBz*tesla) );
    break;
  case ENonUniField:
    magField = new Field3D(fXmin, fBmin, fBmax);
    break;
  case EToscaField3D:
    magField = new TOSCAField3D(fFieldMap, fTXoff, fTYoff, fTZoff);
    break;
  case EBonusBField:
    magField = new BonusBField(fFieldMap, fTXoff, fTYoff, fTZoff, fBScaleFac);
    /*
    for(G4int i=0; i<10; i++){
      pos[0] = 0;
      pos[1] = i*0.6;
      pos[2] = i*0.6;
      magField->GetFieldValue(pos,b);
      printf("%g %g %g %g %g\n",pos[1],pos[2],b[0],b[1],b[2]);
      //GetFieldStrength(Double_t* pos, Double_t b)
    }
    */
    break;
  case EToscaField2D:
    magField = new TOSCAField2D(fFieldMap, fTXoff, fTYoff, fTZoff);
    for(G4int i=0;i<100;i++){
      pos[0] = 0.0;
      pos[1] = 0.0;
      pos[2] = -1500 + i*30;
      magField->GetFieldValue(pos,b);
      printf("%g, %g, %g, %g\n",pos[2],b[0],b[1],b[2]);
    }
    for(G4int i=0;i<100;i++){
      pos[0] = -500 + i*10;
      pos[1] = 0.0;
      pos[2] = 0.0;
      magField->GetFieldValue(pos,b);
      printf("%g, %g, %g, %g\n",pos[0],b[0],b[1],b[2]);
    }
    break;
  }
  magField->GetFieldValue(pos,b);
  printf("Value of Magnetic Field at (0,0.01,0) is Bx:%g  By:%g Bz:%g\n",
	 b[0],b[1],b[2]);
  //
  G4Mag_UsualEqRhs* Equation = new G4Mag_UsualEqRhs(magField);
  //G4MagIntegratorStepper* Stepper = new G4NystromRK4(Equation);
  G4MagIntegratorStepper* Stepper = new G4HelixImplicitEuler(Equation);
  //G4MagIntegratorStepper* Stepper = new G4ClassicalRK4(Equation);
  G4double minStep = 0.10*mm;
  G4ChordFinder* chordF = new G4ChordFinder(magField,minStep,Stepper);
  //
  magFieldMgr = new G4FieldManager(magField);
  magFieldMgr->SetDetectorField(magField);
  //magFieldMgr->CreateChordFinder(magField);
  magFieldMgr->SetChordFinder(chordF);
  gapLog->SetFieldManager(magFieldMgr,true);
}
