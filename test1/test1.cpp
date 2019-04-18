// test1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "CUDACOMPLEX.h"
#include "CUDACOMPLEX.cpp"
#include <time.h>
typedef float DataType;
#include "DeviceData_OptCarroll.h"
#define CUDART_PI_F 3.141592654f
//#define includeAC
int main()
{
	int sum = 0;
	int ArraySize = 10000;
	int threadnum = 1;
	int slicePerThread = ArraySize / threadnum;
	omp_set_num_threads(threadnum);
	clock_t start, finish;
	start = clock();
	CudaComplex *Prog_Vec, *Regr_Vec, *sumProg, *sumRegr, *Prog_AC, *Regr_AC;
	Prog_Vec = (CudaComplex*)malloc(ArraySize * sizeof(CudaComplex));
	Regr_Vec = (CudaComplex*)malloc(ArraySize * sizeof(CudaComplex));
	sumProg = (CudaComplex*)malloc(ArraySize * sizeof(CudaComplex));
	sumRegr = (CudaComplex*)malloc(ArraySize * sizeof(CudaComplex));
	Prog_AC = (CudaComplex*)malloc(ArraySize * sizeof(CudaComplex));
	Regr_AC = (CudaComplex*)malloc(ArraySize * sizeof(CudaComplex));

#pragma omp parallel
	{
		//printf("threadID:%d\n",omp_get_thread_num());
		const DataType DeltaEcb = (EnergyGapSCH - EnergyGapWL)*0.4 * 3;                                      //[NumStates, NumPops] Energy difference between CB states(SCH - WL; WL - ES2; ES2 - ES1; ES1 - GS)[eV]

		const DataType DeltaEcbWL2ES1 = (EnergyGapWL - EnergyGapQDES2)*0.77 / 1.4 / 2 + (EnergyGapQDES2 - EnergyGapQDES1)*0.35;
		const DataType DeltaEcbES12GS = (EnergyGapQDES1 - EnergyGapQDGS)*0.62*1.3;      //[1, 1] Energy difference between ES1 and GS,


		// Definition of commonly used phisical constants
		const DataType Constants_h = 4.1357e-6;                                          //[1, 1] Planck constant[eV ns]
		const DataType Constants_hcut = 4.1357e-6 / (2 * CUDART_PI_F);                               //[1, 1] Normalized Plank constant[eV ns rad^-1]
		const DataType Constants_e = 1.6e-19;                                            //[1, 1] Electron charge[C]
		const DataType Constants_Kb = 8.617343e-5;                                       //[1, 1] Boltzmann constant[eV K^-1]
		const DataType Constants_c = 2.997925e5;                                           //[1, 1] Vacuum speed of light[um ns^-1] 3e8m / s = 3e8*1e6um / 1e9ns = 3e5um / ns
		const DataType KbTemp = Constants_Kb * Temperature;                                            //[1, 1] Product of Boltzmann constant and temperature[eV]

			// Convert the WL effective mass to different units :
		// 1kg = 6.2415e24 eV*ns ^ 2 / um ^ 2 (Electron volt per nanosecond ^ 2 / micron ^ 2)
												   //[1, 1] Effective mass of electrons on the WL[eV*ns ^ 2 / um ^ 2]
		const DataType DOSWLcb = (me_wl*6.241e24*KbTemp) / (CUDART_PI_F*Constants_hcut*Constants_hcut);                               //[1, 1] WL effective density of states for electrons[1 / um ^ 2]

		// Use the group index, if available, to calculate the group velocity.
		const DataType GroupVelocity = Constants_c / nr;                                      //[1, 1] Group velocity(assign to group index the effective refractive index)[um / ns]


			//Terminal left(z = 0) and right(z = L) facets
		const DataType r0 = sqrt(R0);                                                           //[1, 1] Field reflection coefficient at the facet in z = 0[>= 0, <= 1]
		const DataType rL = sqrt(RL);                                                           //[1, 1] Field reflection coefficient at the facet in z = L[>= 0, <= 1]
		const DataType tL = sqrt(1 - RL);
		// Discretization step in space and time
		const size_t NumSlices = (size_t)(round(DeviceLength / (GroupVelocity*dt)));                //[1, 1] Fractional number of discretized slices, rounded to closest integer
		//Recalculate the space step accordingly
		const DataType dz = DataType(DeviceLength) / DataType(NumSlices);                                  //[1, 1] Spatial step[um]

		const DataType Halfdz = dz / 2;                                                                //[1, 1] Half Spatial step[um]
		dt = dz / GroupVelocity;                                                    //[1, 1] Time step[ns]


																		//[1, 1] Update simulation time step stored in Sim structure[ns]
		const size_t NumTimeSteps = (size_t)((TEnd - TStart) / dt);                    //[1, 1] Number of temporal steps[#]

		const DataType NetCurrentCoeff =                                                          //[1, NumSlices] Conversion coefficient from current[mA] to carriers in each slice[mA^-1 ns^-
			Etai / 1000 / 1e9 / NumSlices / Constants_e;
		// DOTS NUMBER
		const DataType DotTotali = Nd * NumLayers*Width*dz;                                      //[1, NumSlices] Total number of dots in each slice[> 0]
		const DataType DotTotaliES1 = DotTotali;                                                     //[1, NumSlices] Number of dots in each slice for each QD population in the ES1[> 0]
		const DataType DotTotaliGS = DotTotali;                                                      //[1, NumSlices] Number of dots in each slice for each QD population in the GS[> 0]
		const DataType DotNumberByDegeneracyES1 = DotTotaliES1 * DegeneracycbES1;                      //[1, NumSlices] Number of dots in the ES multiplied by the corresponding degeneracy[> 0]
		const DataType DotNumberByDegeneracyGS = DotTotaliGS * DegeneracycbGS;                         //[1, NumSlices] Number of dots in the GS multiplied by the corresponding degeneracy[> 0]
		// ENERGY LEVELS


		const DataType ReferenceWavelength = EnergyGapQDGS / Constants_h / Constants_c;                       //[1, 1] Simulation reference wavelength[um]
		const DataType EnergyGS = EnergyGapQDGS;

		//CHARACTERISTIC TIMES

		//escape times
		const DataType TauEsccbES12WL =                                                           //[1, 1] Escape time from ES2 to WL[ns]
			(DegeneracycbES1*Nd / DOSWLcb)*TauCapcbWL2ES1*exp(DeltaEcbWL2ES1 / (0.026 / 300 * Temperature)); //0.026 = > Boltzmann constant multiplied temperature @300K[eV]
		const DataType TauEsccbGS2ES1 =                                                           //[1, 1] Escape time from GS to ES1[ns]
			(DegeneracycbGS / DegeneracycbES1)*TauCapcbES12GS*
			exp(DeltaEcbES12GS / (0.026 / 300 * Temperature));
		const DataType OneOverTauEsccbGS2ES1 = 1 / TauEsccbGS2ES1;                                     //[1, 1] Reciprocal escape time from GS to WL[1 / ns]
		// CARRIERS INITIALIZATION
		DataType *NcbWL = (DataType*)malloc(slicePerThread * sizeof(DataType));                                                  //[1, NumSlices] Number of carriers in WL[>= 0]
		DataType *NcbES1 = (DataType*)malloc(slicePerThread * sizeof(DataType));                                                  //[NumPops, NumSlices] Number of carriers in ES1[>= 0]
		DataType *NcbGS = (DataType*)malloc(slicePerThread * sizeof(DataType));                                                   //[NumPops, NumSlices] Number of carriers in GS[>= 0]
		DataType *RhocbES1 = (DataType*)malloc(slicePerThread * sizeof(DataType));                                  //[NumPops, NumSlices] Occupation probability for ES1[>= 0, <= 1]
		DataType *RhocbGS = (DataType*)malloc(slicePerThread * sizeof(DataType));                                     //[NumPops, NumSlices] Occupation probability for GS[>= 0, <= 1]
		memset(NcbWL, 0., slicePerThread * sizeof(DataType));
		memset(NcbES1, 0., slicePerThread * sizeof(DataType));
		memset(NcbGS, 0., slicePerThread * sizeof(DataType));
		memset(RhocbES1, 0., slicePerThread * sizeof(DataType));
		memset(RhocbGS, 0., slicePerThread * sizeof(DataType));

		//CudaComplex NcbWLAC(0., 0.);                                                //[1, NumSlices] Number of carriers in WL[>= 0]
		//CudaComplex NcbES1AC(0., 0.);                                                  //[NumPops, NumSlices] Number of carriers in ES1[>= 0]
		//CudaComplex NcbGSAC(0., 0.);                                                   //[NumPops, NumSlices] Number of carriers in GS[>= 0]
		//CudaComplex RhocbES1AC(0., 0.);                                  //[NumPops, NumSlices] Occupation probability for ES1[>= 0, <= 1]
		//CudaComplex RhocbGSAC(0., 0.);                                     //[NumPops, NumSlices] Occupation probability for GS[>= 0, <= 1]


		const DataType NetCurrentCB = Pilot * NetCurrentCoeff;  //[1, NumSlices] Injection rate[ns^-1]

			//FIELD CONFINEMENT FACTORS
		const DataType GammaX = gamma_x;
		const DataType GammaYRef = GammaY;
		const DataType GammaXYRef = GammaX * GammaYRef;
		const DataType GammaXY = GammaXYRef;
		const DataType GammaYCorrection = 1.;
		const DataType StimulatedEmissionCoeffGSXYRef = 2 * GroupVelocity*GammaX / EnergyGS * GammaYRef;
		const DataType StimulatedEmissionCoeffGSXY = StimulatedEmissionCoeffGSXYRef;

		//Losses & Saturated material gain
		const DataType MaterialLossedHalfed = Alfa_i / 2 / 1e4;
		const DataType SaturableLossesHalfed = SaturableLosses / 2.0 / 1e4;
		const DataType SaturationCoefficient = 1.0e-3;
		//gain broadening
		const DataType HomoGS = BWhomogeneousGS;
		const DataType ExpMinusHomoGSSingledt = exp(-BWhomogeneousGS * dt);
		const DataType RatioExpGS = (1 - ExpMinusHomoGSSingledt) / (1 + ExpMinusHomoGSSingledt);
		const CudaComplex Homo_dt(-HomoGS * dt, 0);
		const CudaComplex FilterCoeffGS = exp(Homo_dt);

		const DataType GainNormFactor = 1e8*Nd / H_WELL / Constants_hcut;
		const DataType MaterialGainGS = Gain * GainNormFactor / BWhomogeneousGS;
		DataType **GainGS = (DataType**)malloc(sizeof(DataType*)*slicePerThread);
		DataType **ModalGainGS = (DataType**)malloc(sizeof(DataType*)*slicePerThread);
		DataType *ModalGainGSpp = (DataType*)malloc(slicePerThread * sizeof(DataType));
		CudaComplex **GainGSAC = (CudaComplex**)malloc(slicePerThread * sizeof(CudaComplex*));
		CudaComplex **ModalGainGSAC = (CudaComplex**)malloc(slicePerThread * sizeof(CudaComplex*));
		CudaComplex *ModalGainGSACpp = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));
		for (int i = 0; i < slicePerThread; i++) {
			*(GainGS + i) = (DataType*)malloc(sizeof(DataType) * 2);
			*(ModalGainGS + i) = (DataType*)malloc(sizeof(DataType) * 2);
			*(GainGSAC + i) = (CudaComplex*)malloc(sizeof(CudaComplex) * 2);
			*(ModalGainGSAC + i) = (CudaComplex*)malloc(sizeof(CudaComplex) * 2);
			GainGS[i][0] = 0.5*MaterialGainGS*(2 * RhocbGS[0] - 1);
			GainGS[i][1] = 0.5*MaterialGainGS*(2 * RhocbGS[0] - 1);
			ModalGainGS[i][0] = GammaXYRef * GainGS[i][0];
			ModalGainGS[i][1] = GammaXYRef * GainGS[i][0];
			GainGSAC[i][0] = CudaComplex(0., 0.);
			GainGSAC[i][1] = CudaComplex(0., 0.);
			ModalGainGSAC[i][0] = CudaComplex(0., 0.);
			ModalGainGSAC[i][1] = CudaComplex(0., 0.);
			ModalGainGSACpp[i] = CudaComplex(0., 0.);
		}
		const DataType GainCompressionGSCoeffRef = GainCompressionFactorGS / (Width*H_WELL*dz*1e-12)*(GammaXYRef / EnergyGS);;

		//spontaneous emission
		const DataType sqrtBetaOver2pi = sqrt(SpontaneousEmissionCoefficientQD / GammaYRef / (2 * CUDART_PI_F));
		//electric field and random noise
		const DataType Photons2Power = (Constants_e*(1 / dt) / 1e-9)*1e3;
		CudaComplex **EProg = (CudaComplex**)malloc(slicePerThread * sizeof(CudaComplex*));
		CudaComplex **ERegr = (CudaComplex**)malloc(slicePerThread * sizeof(CudaComplex*));

		DataType *PowerProg = (DataType*)malloc(slicePerThread * sizeof(DataType));
		DataType *PowerRegr = (DataType*)malloc(slicePerThread * sizeof(DataType));
		memset(PowerProg, 0., slicePerThread * sizeof(DataType));
		memset(PowerRegr, 0., slicePerThread * sizeof(DataType));


		CudaComplex **IfieldGSProg = (CudaComplex**)malloc(slicePerThread * sizeof(CudaComplex*));
		CudaComplex **IfieldGSRegr = (CudaComplex**)malloc(slicePerThread * sizeof(CudaComplex*));


		CudaComplex *PhaseNoiseGSProgNew = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));
		CudaComplex *PhaseNoiseGSRegrNew = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));


		CudaComplex *IspontGSProgNew = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));
		CudaComplex *IspontGSRegrNew = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));


		CudaComplex*EnergyGSProgProg = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));
		CudaComplex*EnergyGSRegrRegr = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));
		CudaComplex*EnergyGSProgRegr = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));
		CudaComplex*EnergyGSRegrProg = (CudaComplex*)malloc(slicePerThread * sizeof(CudaComplex));


		for (int i = 0; i < slicePerThread; i++) {
			*(EProg + i) = (CudaComplex*)malloc(sizeof(CudaComplex) * 2);
			*(ERegr + i) = (CudaComplex*)malloc(sizeof(CudaComplex) * 2);
			*(IfieldGSProg + i) = (CudaComplex*)malloc(sizeof(CudaComplex) * 2);
			*(IfieldGSRegr + i) = (CudaComplex*)malloc(sizeof(CudaComplex) * 2);
			EProg[i][0] = CudaComplex(0., 0.);
			EProg[i][1] = CudaComplex(0., 0.);
			ERegr[i][0] = CudaComplex(0., 0.);
			ERegr[i][1] = CudaComplex(0., 0.);
			IfieldGSProg[i][0] = CudaComplex(0., 0.);
			IfieldGSProg[i][1] = CudaComplex(0., 0.);
			IfieldGSRegr[i][0] = CudaComplex(0., 0.);
			IfieldGSRegr[i][1] = CudaComplex(0., 0.);


			PhaseNoiseGSProgNew[i] = CudaComplex(0., 0.);
			PhaseNoiseGSRegrNew[i] = CudaComplex(0., 0.);
			IspontGSProgNew[i] = CudaComplex(0., 0.);
			IspontGSRegrNew[i] = CudaComplex(0., 0.);
			EnergyGSProgProg[i] = CudaComplex(0., 0.);
			EnergyGSRegrRegr[i] = CudaComplex(0., 0.);
			EnergyGSProgRegr[i] = CudaComplex(0., 0.);
			EnergyGSRegrProg[i] = CudaComplex(0., 0.);

		}

		//initialization of indices before simulation
		int pp = 1;
		int cc = 0;
		CudaComplex C1, C2;
		CudaComplex TmpRegr, TmpRegrAC, TmpProg, TmpProgAC, sumEspRegr, sumEspProg;
		size_t idx = omp_get_thread_num();
		//printf("idx:%d\n",idx);


		for (size_t it = 1; it <= NumTimeSteps; it++) {
			for (int j = 0; j < slicePerThread; j++) {
				//printf("success1!\n");
				DataType current_time = it * dt;

				// INJECTED CURRENT

				DataType CurrentPower_eV = real(EnergyGSProgProg[j] + EnergyGSRegrRegr[j]);
				// CARRIER RATE EQUATION(see Ref. 3, eq. 6c - 6h)

				//////////% WL carriers dynamics(Ref. 3, eq. 6g - 6h) % //////////%
				DataType NonRadiativeRecombinationWL = NcbWL[j] * OneOverTauNRcbWL;                    //[1, NumSlices] Non radiative recomb.rate in WL[ns^-1]
				DataType CaptureFromWL2ES1 = ((1 - RhocbES1[j])*NcbWL[j])*OneOverTauCapcbWL2ES1;         //[NumPops, NumSlices] Capture rate from WL to ES1[ns^-1]
				DataType EscapeFromES12WL = NcbES1[j] / TauEsccbES12WL;                                //[NumPops, NumSlices] Escape rate from ES1 to WL[ns^-1]
				DataType RecombinationWL = NonRadiativeRecombinationWL;                            //[1, NumSlices] In the WL we have only non radiative recombinations[ns^-1]

#ifdef includeAC
				CaptureFromWL2ES1 = CaptureFromWL2ES1 - 2 * real(RhocbES1AC*conj(NcbWLAC))*OneOverTauCapcbWL2ES1;
				CudaComplex RecombinationWLAC = NcbWLAC * OneOverTauNRcbWL;
				CudaComplex EscapeFromES1AC2WLAC = NcbES1AC / TauEsccbES12WL;
				CudaComplex CaptureFromES1AC2WLAC = ((1 - RhocbES1)*NcbWLAC - RhocbES1AC * NcbWL)*OneOverTauCapcbWL2ES1;
				CudaComplex dNcbWLAC = -RecombinationWLAC + (EscapeFromES1AC2WLAC - CaptureFromES1AC2WLAC);
#endif
				DataType dNcbWL =                                                               //[1, NumSlices] Variation rate for WL carrier number[ns^-1]
					NetCurrentCB                                                     //[1, NumSlices] Injection rate[ns^-1]
					- RecombinationWL                                                 //[1, NumSlices] Non radiative recombinations in the WL[ns^-1]
					+ (EscapeFromES12WL                                               //[NumPops, NumSlices] Escape rate from ES1 to WL[ns^-1]
						- CaptureFromWL2ES1);                                                //[NumPops, NumSlices] Capture rate from WL to ES1[ns^-1]
					///////////% ES1 carriers dynamics(Ref. 3, eq. 6e-6f) % //////////%
				DataType SpontaneousEmissionES1 = NcbES1[j] * RhocbES1[j]*OneOverTauSpontES1;             //[NumPops, NumSlices] Spontaneous emission rate for ES1[ns^-1]
				DataType EscapeFromGS2ES1 = ((1 - RhocbES1[j])*NcbGS[j])*OneOverTauEsccbGS2ES1;          //[NumPops, NumSlices] Escape rate from GS to ES1[ns^-1]
				DataType CaptureFromES12GS = (NcbES1[j]*(1 - RhocbGS[j]))*OneOverTauCapcbES12GS;         //[NumPops, NumSlices] Capture rate from ES1 to GS[ns^-1]
				DataType NonRadiativeRecombES1 = NcbES1[j] * OneOverTauNRcbES1;                         //[NumPops, NumSlices] Non radiative recombination rate for ES1[ns^-1]
#ifdef includeAC
				SpontaneousEmissionES1 = SpontaneousEmissionES1 + 2 * real(RhocbES1AC*conj(NcbES1AC))*OneOverTauSpontES1;
				CudaComplex SpontaneousEmissionES1AC = (NcbES1*RhocbES1AC + NcbES1AC * RhocbES1)*OneOverTauSpontES1;
				EscapeFromGS2ES1 = EscapeFromGS2ES1 - 2 * real(RhocbES1AC*conj(NcbGSAC))*OneOverTauEsccbGS2ES1;
				CudaComplex EscapeFromGSAC2ES1AC = NcbGSAC * (1 - RhocbES1) - NcbGS * RhocbES1AC*OneOverTauEsccbGS2ES1;
				CaptureFromES12GS = CaptureFromES12GS - 2 * real(RhocbGSAC*conj(NcbES1AC))*OneOverTauCapcbES12GS;
				CudaComplex CaptureFromES1AC2GSAC = (NcbES1*(-RhocbGSAC) + NcbES1AC * (1 - RhocbGS))*OneOverTauCapcbES12GS;
				CudaComplex NonRadiativeRecombES1AC = NcbES1AC * OneOverTauNRcbES1;
				CudaComplex RecombinationsES1AC = NonRadiativeRecombES1AC + SpontaneousEmissionES1AC;
				CudaComplex dNcbES1AC = CaptureFromES1AC2WLAC - RecombinationsES1AC - EscapeFromES1AC2WLAC + EscapeFromGSAC2ES1AC - CaptureFromES1AC2GSAC;
#endif // includeAC
				DataType RecombinationsES1 =                                                    //[NumPops, NumSlices] Total recombination rate for ES1[ns^-1]
					NonRadiativeRecombES1 +                                            //[NumPops, NumSlices] Non radiative recombination rate for ES1[ns^-1]
					SpontaneousEmissionES1;                                             //[NumPops, NumSlices] Spontaneous emission rate for ES1[ns^-1]
				DataType dNcbES1 =                                                             //[NumPops, NumSlices] Variation rate for ES1[ns^-1]
					-RecombinationsES1                                               //[NumPops, NumSlices] Total recombination rate for ES1[ns^-1]
					+ CaptureFromWL2ES1                                               //[NumPops, NumSlices] Capture rate from ES1 to GS[ns^-1]
					- EscapeFromES12WL                                                //[NumPops, NumSlices] Escape rate from ES1 to WL[ns^-1]
					+ EscapeFromGS2ES1                                                //[NumPops, NumSlices] Escape rate from GS to ES1[ns^-1]
					- CaptureFromES12GS;                                                 //[NumPops, NumSlices] Capture rate from ES1 to GS[ns^-1]
					///// Variation of the GS carriers(Ref. 3, eq. 6c - 6d) % ////
				DataType SpontaneousEmissionGS = NcbGS[j] * RhocbGS[j]*OneOverTauSpontGS;                 //[NumPops, NumSlices] Spontaneous emission rate for GS[ns^-1]
				DataType StimulatedEmissionGS = StimulatedEmissionCoeffGSXY *
					(GainGS[j][cc] * real(EnergyGSProgProg[j] + EnergyGSRegrRegr[j]) + real(GainGSAC[j][cc] * conj(EnergyGSProgRegr[j])) + real(conj(GainGSAC[j][cc])*conj(EnergyGSRegrProg[j])));     //[NumPops, NumSlices] Stimulated emission rate for GS[ns^-1]

				DataType NonRadiativeRecombGS = NcbGS[j] * OneOverTauNRcbGS;                            //[NumPops, NumSlices] Non radiative recombination rate for GS[ns^-1]

#ifdef includeAC
				SpontaneousEmissionGS = SpontaneousEmissionGS + 2 * real(RhocbGSAC*conj(NcbGSAC))*OneOverTauSpontGS;
				CudaComplex SpontaneousEmissionGSAC = (NcbGS*RhocbGSAC + NcbGSAC * RhocbGS)*OneOverTauSpontGS;
				CudaComplex StimulatedEmissionGSAC = 0.5*StimulatedEmissionCoeffGSXY*
					(GainGS[cc] * (EnergyGSProgRegr + conj(EnergyGSRegrProg)) + GainGSAC[cc] * (EnergyGSProgProg + conj(EnergyGSRegrRegr)));
				CudaComplex NonRadiativeRecombGSAC = NcbGSAC * OneOverTauNRcbGSAC;
				CudaComplex RecombinationsGSAC = NonRadiativeRecombGSAC + SpontaneousEmissionGSAC + StimulatedEmissionGSAC;
				CudaComplex dNcbGSAC = CaptureFromES1AC2GSAC - RecombinationsGSAC - EscapeFromGSAC2ES1AC;
#endif // includeAC

				DataType RecombinationsGS =                                                     //[NumPops, NumSlices] Total recombination rates for GS[ns^-1]
					NonRadiativeRecombGS +                                             //[NumPops, NumSlices] Non radiative recombination rate for GS[ns^-1]
					SpontaneousEmissionGS +                                            //[NumPops, NumSlices] Spontaneous emission rate for GS without grating[ns^-1]
					StimulatedEmissionGS;                                               //[NumPops, NumSlices] Stimulated emission rate for GS[ns^-1]
				DataType dNcbGS =                                                               //[NumPops, NumSlices] Variation rate for GS[ns^-1]
					CaptureFromES12GS                                                //[NumPops, NumSlices] Capture rate from ES1 to GS[ns^-1]
					- RecombinationsGS                                                //[NumPops, NumSlices] Total recombination rates for GS[ns^-1]
					- EscapeFromGS2ES1;                                                  //[NumPops, NumSlices] Escape rate from GS to ES1[ns^-1]
					//update all the carriers numbers
				NcbWL[j] = NcbWL[j] + dt * dNcbWL;
				if (NcbWL[j] < 0)
					NcbWL[j] = 0;                            //[1, NumSlices]   TO DO
				NcbES1[j] = NcbES1[j] + dt * dNcbES1;                                           //[NumPops, NumSlices]  TO DO
				RhocbES1[j] = NcbES1[j] / DotNumberByDegeneracyES1;                              //[NumPops, NumSlices] Occupation probability for ES1[>= 0, <= 1]

				NcbGS[j] = NcbGS[j] + dt * dNcbGS;                                              //[NumPops, NumSlices]      TO DO
				RhocbGS[j] = NcbGS[j] / DotNumberByDegeneracyGS;                                 //[NumPops, NumSlices] Occupation probability for GS[>= 0, <= 1]
				//float myrandf = curand_uniform(&randstate);

#ifdef includeAC
				NcbWLAC = NcbWLAC + dt * dNcbWLAC;
				NcbES1AC = NcbES1AC + dt * dNcbES1AC;
				RhocbES1AC = NcbES1AC / DotNumberByDegeneracyES1;
				NcbGSAC = NcbGSAC + dt * dNcbGSAC;
				RhocbGSAC = NcbGSAC / DotNumberByDegeneracyGS;
#endif // includeAC


				//optical gain
				GainGS[j][pp] = 0.5*(MaterialGainGS*(2 * RhocbGS[j] - 1));

#ifndef includeAC
				GainGS[j][pp] = GainGS[j][pp] / (1 + GainCompressionGSCoeffRef * GammaYCorrection*CurrentPower_eV);
#endif // !includeAC

				ModalGainGSpp[j] = GainGS[j][pp] * GammaXY;
				ModalGainGS[j][pp] = ModalGainGSpp[j];

#ifdef includeAC
				GainGSAC[pp] = 0.5*(MaterialGainGS * 2 * RhocbGSAC);
				ModalGainGSACpp = GainGSAC[pp] * GammaXY;
				ModalGainGSAC[pp] = ModalGainGSACpp;
#endif // includeAC


				//electric field
				//CudaComplex sumEspProg;
				//CudaComplex sumEspRegr;
				if (current_time <= TimeStopEsp)
				{
					CudaComplex PhaseNoiseGSProgOld = PhaseNoiseGSProgNew[j];
					CudaComplex PhaseNoiseGSRegrOld = PhaseNoiseGSRegrNew[j];
					PhaseNoiseGSProgNew[j] = CudaComplex(cos(2.0*(idx * slicePerThread + j )*it / 10), sin(2.0*(idx * slicePerThread + j) / 10 * it));
					PhaseNoiseGSRegrNew[j] = CudaComplex(cos((idx * slicePerThread + j)*it / 10.0), sin((idx * slicePerThread + j) *it / 10.0));
#ifdef xor_rand
					//PhaseNoiseGSProgNew = rand_linear(&seed);
					//PhaseNoiseGSRegrNew = rand_linear(&seed);
#else
					//getRand(&PhaseNoiseGSProgNew, &state, 1);
					//getRand(&PhaseNoiseGSRegrNew, &state, 1);
#endif
					CudaComplex IspontGSProgOld = IspontGSProgNew[j];
					CudaComplex IspontGSRegrOld = IspontGSRegrNew[j];
					IspontGSProgNew[j] = IspontGSProgOld * FilterCoeffGS + RatioExpGS * (FilterCoeffGS * PhaseNoiseGSProgOld + PhaseNoiseGSProgNew[j]);
					IspontGSRegrNew[j] = IspontGSRegrOld * FilterCoeffGS + RatioExpGS * (FilterCoeffGS * PhaseNoiseGSRegrOld + PhaseNoiseGSRegrNew[j]);
					DataType RspGS = (NcbGS[j]*RhocbGS[j])*(EnergyGS*OneOverTauSpontGS);

					DataType sqrtRspGSOverHomoGS = sqrt(RspGS / HomoGS);
					DataType sqrtBetaOver2piCurrent = sqrtBetaOver2pi / sqrt(GammaYCorrection);
					sumEspProg = sqrtBetaOver2piCurrent * sqrtRspGSOverHomoGS * IspontGSProgNew[j];
					sumEspRegr = sqrtBetaOver2piCurrent * sqrtRspGSOverHomoGS * IspontGSRegrNew[j];
				}
				else
				{
					sumEspProg = CudaComplex(0., 0.);
					sumEspRegr = CudaComplex(0., 0.);
				}

				DataType MatLosses = 1 / Halfdz - MaterialLossedHalfed - SaturableLossesHalfed / (1 + SaturationCoefficient * (PowerProg[j] + PowerRegr[j]));
				TmpProg = Halfdz * (MatLosses*EProg[j][cc] + ModalGainGS[j][cc] * IfieldGSProg[j][cc]);
				TmpProgAC = Halfdz * ModalGainGSAC[j][cc] * IfieldGSRegr[j][cc];
				TmpRegr = Halfdz * (MatLosses*ERegr[j][cc] + ModalGainGS[j][cc] * IfieldGSRegr[j][cc]);
				TmpRegrAC = Halfdz * conj(ModalGainGSAC[j][cc])*IfieldGSProg[j][cc];
				Prog_Vec[idx * slicePerThread + j] = TmpProg;
				Regr_Vec[idx * slicePerThread + j] = TmpRegr;
				sumProg[idx * slicePerThread + j] = sumEspProg;
				sumRegr[idx * slicePerThread + j] = sumEspRegr;
				Prog_AC[idx * slicePerThread + j] = TmpProgAC;
				Regr_AC[idx * slicePerThread + j] = TmpRegrAC;

			}
#pragma omp barrier
			for (int j = 0; j < slicePerThread; j++)
			{
				DataType A = Halfdz * (ModalGainGSpp[j]*RatioExpGS - MaterialLossedHalfed);
				CudaComplex AAc = Halfdz * ModalGainGSAC[j][cc] * RatioExpGS;
				CudaComplex Bprog = Halfdz * (ModalGainGSpp[j]*FilterCoeffGS*(IfieldGSProg[j][cc] + RatioExpGS * EProg[j][cc]));
				CudaComplex BprogAC = Halfdz * (ModalGainGSACpp[j]*FilterCoeffGS*(IfieldGSRegr[j][cc] + RatioExpGS * ERegr[j][cc]));
				CudaComplex Bregr = Halfdz * (ModalGainGSpp[j]*FilterCoeffGS*(IfieldGSRegr[j][cc] + RatioExpGS * ERegr[j][cc]));
				CudaComplex BregrAC = Halfdz * (conj(ModalGainGSACpp[j])*FilterCoeffGS*(IfieldGSProg[j][cc] + RatioExpGS * EProg[j][cc]));;
				DataType A11 = 1 - A;
				DataType A22 = A11;
				CudaComplex A12 = -AAc;
				CudaComplex A21 = -conj(AAc);
				CudaComplex detA = A11 * A22 - A12 * A21;

				if (idx*slicePerThread +j == 0)
					C1 = Bprog + BprogAC + r0 * (Regr_Vec[0] + Regr_AC[0]);
				else
					C1 = Bprog + BprogAC + Prog_Vec[idx * slicePerThread + j - 1] + Prog_AC[idx * slicePerThread + j - 1];

				//printf("C1:%f",C1);
				if (idx * slicePerThread + j == ArraySize - 1)
					C2 = Bregr + BregrAC + rL * (Prog_Vec[ArraySize - 1] + Prog_AC[ArraySize - 1]);
				else
					C2 = Bregr + BregrAC + Regr_Vec[idx * slicePerThread + j + 1] + Regr_AC[idx * slicePerThread + j + 1];

				EProg[j][pp] = (C1*A22 - A12 * C2) / detA;
				ERegr[j][pp] = (C2*A11 - A21 * C1) / detA;

				if (idx * slicePerThread + j == 0)
					EProg[j][pp] = EProg[j][pp] + r0 * sumRegr[0];
				else
					EProg[j][pp] = EProg[j][pp] + sumProg[idx * slicePerThread + j - 1];
				if (idx * slicePerThread + j == ArraySize - 1)
					ERegr[j][pp] = ERegr[j][pp] + rL * sumProg[ArraySize-1];
				else
					ERegr[j][pp] = ERegr[j][pp] + sumRegr[idx * slicePerThread + j + 1];



				//SPLIT STEP

				PowerProg[j] = abs(EProg[j][pp])*abs(EProg[j][pp])*Photons2Power;
				PowerRegr[j] = abs(ERegr[j][pp])*abs(ERegr[j][pp])*Photons2Power;

				//FILTERED FIELDS
				CudaComplex EProgppExpanded = EProg[j][pp];
				CudaComplex ERegrppExpanded = ERegr[j][pp];
				IfieldGSProg[j][pp] = IfieldGSProg[j][cc] * FilterCoeffGS + RatioExpGS * (FilterCoeffGS*EProg[j][cc] + EProgppExpanded);
				IfieldGSRegr[j][pp] = IfieldGSRegr[j][cc] * FilterCoeffGS + RatioExpGS * (FilterCoeffGS*ERegr[j][cc] + ERegrppExpanded);
				EnergyGSProgProg[j] = EProgppExpanded * conj(IfieldGSProg[j][pp]);
				EnergyGSRegrRegr[j] = ERegrppExpanded * conj(IfieldGSRegr[j][pp]);
#ifdef includeAC
				EnergyGSProgRegr = EProgppExpanded * conj(IfieldGSRegr[pp]);
				EnergyGSRegrProg = ERegrppExpanded * conj(IfieldGSProg[pp]);
#endif // includeAC
				
			}
			//Shift vectors
			pp = 1 - pp;
			cc = 1 - cc;
#pragma omp barrier
		}

		
		

		free(PhaseNoiseGSProgNew);
		free(PhaseNoiseGSRegrNew);
		free(IspontGSProgNew);
		free(IspontGSRegrNew);
		free(EnergyGSProgProg);
		free(EnergyGSRegrRegr);
		free(EnergyGSProgRegr);
		free(EnergyGSRegrProg);
		free(PowerProg);
		free(PowerRegr);
		free(ModalGainGS);
		free(ModalGainGSACpp);
		
	}


	for (int i = 0; i < ArraySize; i++) {
		Regr_Vec[i].display();
	}
	finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("%f seconds\n", duration);
	
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
