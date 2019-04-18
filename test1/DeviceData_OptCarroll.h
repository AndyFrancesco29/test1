//const int NumQDPopulations = 11;                             //[1, 1](Always odd) number of QD population created in EnergyGapQD[]
const DataType EnergyGapSCH = 1.2797;                         //[1, 1] Energy gap between CB and VB band edges in the SCH[eV]
const DataType EnergyGapWL = 1.1539;                         //[1, 1] Energy gap between CB and VB band edges in the WL[eV]
const DataType EnergyGapQDES2 = 1.114;
const DataType EnergyGapQDES1 = 1.0543;
const DataType EnergyGapQDGS = 0.9879;          //[NumStates, 1] Characteristic interband transition energies for the QD central(most likely) QD population(ES2, ES1, GS)[eV]
const DataType InhomDeltaEES2 = 0.025479654008641;
const DataType InhomDeltaEES1 = 0.019109740506480;
const DataType InhomDeltaEGS = 0.016986436005760;//[0.060; 0.045; 0.040] / sqrt(8 * log(2)); //[NumStates, 1] FWHM->sigma of inhomogeneous broadening of the QD ensemble(ES2, ES1, GS)[eV]
//const int CentralPop = (NumQDPopulations + 1) / 2;

const DataType EMass = 9.11e-31;                                             //[1, 1] Electron mass[kg]

// // Device parameters
const DataType Width = 5;               // Ridge width[um]
const DataType DeviceLength = 30;                         //[1, 1] total cavity length[um]
const DataType nr = 3.3445;                      //[1, 1] effective refractive index[> 0]
const DataType NumLayers = 15;                          //[1, 1] QD layer number[]
const DataType Alfa_i = 2;                     //[f(z)->1, 1] Intrinsic waveguide losses[cm^-1]
const DataType Temperature = 273 + 20;                     //[1, 1] temperature[K]
const DataType R0 = 0.32 * 2;                        //[1, 1] power reflectivity in z = 0 facet[>= 0 <= 1]
const DataType RL = 0.32 * 2;                        //[1, 1] power reflectivity in z = L facet[>= 0 <= 1]
const DataType hQD = 5;                           //[1, 1] QD height[nm]
const DataType Nd = 270;                         //[1, 1] QD surface density[um^-2]
const DataType DegeneracycbES2 = 6;
const DataType DegeneracycbES1 = 4;
const DataType DegeneracycbGS = 2;					//[NumStates, 1] QD state degeneracy(ES2, ES1, GS)[#]
const DataType GammaY = 0.95;                    //[f(z)->1, NumSlice] ridge transverse confinement factor[> 0]
const DataType Etai = 0.7;                         //[1, 1] Internal quantum efficiency[> 0]
const DataType H_WELL = 5e-3;                        //[1, 1] Quantum Well width[um]
const DataType gamma_x = 0.127;                       //[1, 1] field confinement factor in the QD layers in the growth direction[> 0]

//non radiative times
const DataType TauNRcbWL = 1.2;
const DataType OneOverTauNRcbWL = 1. / TauNRcbWL;                                              //[1, 1] Reciprocal non raditive recombination time for WL[1 / ns]
const DataType TauNRcbES2 = 1e50;
const DataType OneOverTauNRcbES2 = 1. / TauNRcbES2;                                              //[1, 1] Reciprocal non raditive recombination time for WL[1 / ns]
const DataType TauNRcbES1 = 1e50;
const DataType OneOverTauNRcbES1 = 1. / TauNRcbES1;                                              //[1, 1] Reciprocal non raditive recombination time for WL[1 / ns]
const DataType TauNRcbGS = 1e50;																//[NumStates + 2] Non radiative recombinations in SCH, QW and QD states(ES2, ES1, GS) in cb[ns]
const DataType OneOverTauNRcbGS = 1. / TauNRcbGS;                                              //[1, 1] Reciprocal non raditive recombination time for WL[1 / ns]
const DataType TauNRcbGSAC = 1e50;
const DataType OneOverTauNRcbGSAC = 1. / TauNRcbGSAC;

//[NumStates + 1] Capture and relaxation times between SCH, WL and QD states(ES2, ES1, GS) in cb[ns]

//capture times
const DataType TauCapcbWL2ES1 = 3e-3 + 3e-3;                                   //[1, 1] Capture time from WL to ES2[ns]
const DataType OneOverTauCapcbWL2ES1 = 1. / TauCapcbWL2ES1;                                    //[1, 1] 1 / Capture time from WL to ES2[ns^-1]
const DataType TauCapcbES12GS = 1e-3;                                                //[1, 1] Capture time from ES1 to GS[ns]
const DataType OneOverTauCapcbES12GS = 1. / TauCapcbES12GS;

//spontaneous emission times
const DataType TauSpontES1 = 1;
const DataType OneOverTauSpontES1 = 1 / TauSpontES1;
const DataType TauSpontGS = 2;                     //[NumStates, 1] Spontaneous emission time in the QD states(ES2, ES1, GS)[ns]
const DataType OneOverTauSpontGS = 1 / TauSpontGS;

//const DataType InhomDeltaE = InhomDeltaE;                 //[NumStates, 1] sigma of inhomogeneous broadening of the QD ensemble[eV]
//const DataType NumQDPopulations = NumQDPopulations;            //[1, 1] Number of QD population(Note: this number must be odd!!!)[#]
const DataType Gain = 0.1335e-15;//[NumStates, 1] QD states gain(ES1, GS)[cm ^ 2 eV].Dipole moment is sqrt(c*n_eff*eps0*'Gain' / (omega_0*degeneracy))
const DataType SpontaneousEmissionCoefficientQD = 3e-2;                 //[1, 1] Spontaneous emission coupling coefficients[> 0]
const DataType me_SCH = 0.067*EMass;                 //[1, 1] electron effective mass in the SCH barrier[kg]
const DataType me_wl = 0.02*EMass;                  //[1, 1] electron effective mass in the QW[kg]
const DataType me_QD = 0.02*EMass;                  //[1, 1] electron effective mass in the QD[kg]
const DataType mh_QD = 0.21*EMass * 2;                //[1, 1] QD hole effective mass[kg]
const DataType mh_wl = 0.21*EMass * 2;                //[1, 1] QW hole effective mass[kg]
const DataType mh_SCH = 0.47*EMass;                  //[1, 1] SCH hole effective mass[kg]

const DataType BWhomogeneousGS = 0.5*15.2e3*4.5*1.3;                  //[1, 1] GS homogeneous broadening, single side[rad / ns]

DataType dt = 30e-6; ;                      //[1, 1] Simulation time step[ns]
const DataType TStart = 0;                      //[1, 1] Simulation start time.This is the origin for the time vector.[ns]
const DataType TEnd = 1;                        //[1, 1] Simulation end time; it must be > TStart.[ns]
const DataType Pilot = 200;                     //[1, 1] or [f(e, t)->[NumSlices, 1]] Injected current[mA]
const int considerphotons = 1;             //[1, 1] Flag; 1: the simulation must also solve for photon / field; 0: the simulation calculates carrier distribution without photons[]

const DataType TimeStopEsp = 10000;
const DataType SaturableLosses = 2.0;
const DataType GainCompressionFactorGS = 0.;

