double CATS3D::GetTotMomentum3D(const unsigned& WhichMomBin_out, const unsigned& WhichMomBin_side, const unsigned& WhichMomBin_long) const{
    if(NumMomBins<=WhichMomBin) return 0;
    return std::sqrt(std::pow(MomBinCenter_out[WhichMomBin_out],2) + 
                     std::pow(MomBinCenter_side[WhichMomBin_side],2) +
                     std::pow(MomBinCenter_long[WhichMomBin_long],2));
    //return 0.5*(MomBin[WhichMomBin]+MomBin[WhichMomBin+1]);
}

double CATS3D::CoulombPotential(const double& Radius) const{
    return Gamow?0:Q1Q2*AlphaFS/(fabs(Radius)+1e-64);
}

//the differential equation for the Schroedinger equation
void CATS3D::PropagatingFunction(double& Basic, double& Full,
                                 const double& Radius, const double& Momentum,
                                 const unsigned short& usPW, const unsigned short& usCh){
    //make sure that there is no division by zero by adding 1e-64
    //the Basic result is the Prop.Fun. WITHOUT a short range potential

    Basic = 2*RedMass*CoulombPotential(Radius) + double(usPW)*(double(usPW)+1)/(Radius*Radius+1e-64) - Momentum*Momentum;
    //the Full result is the Prop.Fun. WITH a short range potential
    double* Parameters = NULL;
    if(PotPar[usCh][usPW]){
        PotPar[usCh][usPW]->SetVariable(0,Radius*NuToFm,true);
        PotPar[usCh][usPW]->SetVariable(1,Momentum,true);
        Parameters = PotPar[usCh][usPW]->GetParameters();
    }
    //else if(PotParArray[usCh][usPW]){
    //    PotParArray[usCh][usPW][0] = Radius*NuToFm;
    //    PotParArray[usCh][usPW][1] = Momentum;
    //    Parameters = PotParArray[usCh][usPW];
    //}
    else{
        return;
    }
if(DEBUG==true)
{
printf("  Basic=%.15e\n",Basic);
printf("  ShortRangePotential[%u][%u]=%.15e\n",usCh,usPW,ShortRangePotential[usCh][usPW](Parameters));
printf("  Parameters[0]=%.15e\n",Parameters[0]);
printf("  Parameters[1]=%.15e\n",Parameters[1]);
printf("  Parameters[2]=%.15e\n",Parameters[2]);
printf("  Full=%.15e\n",Full);
}
    //! In principle this should be executed only if ShortRangePotential[usCh][usPW] is defined.
    //! Do note that this function is NEVER called in case this is not true! Make sure that this stays so!
    Full = Basic + 2*RedMass*ShortRangePotential[usCh][usPW](Parameters);
}

void CATS3D::ComputeWaveFunction(){

    //the sofia model does not require to set the WF, as we will work with a simple
    //phase shift which will be implemented later when we get the effective WF
    // if(use_sofia_model){
    //     ComputedWaveFunction = true;
    //     return;
    // }

    if(!MomBinConverged){
        MomBinConverged = new bool [NumMomBins_out][NumMomBins_side][NumMomBins_long];
    }
    for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
        for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
            for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                MomBinConverged[uMomBin_out][uMomBin_side][uMomBin_long] = true;
            }
        }
    }

    //here one finds the max num PWs, this is needed later on for
    //the correct mapping of all variables
    unsigned short MaxNumPW=0;
    //we also find the total number of steps to iterate over
    unsigned TotNumSteps=0;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        if(NumPW[usCh]>MaxNumPW){
            MaxNumPW = NumPW[usCh];
        }
        for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
            //if(!ShortRangePotential[usCh][usPW] && !ExternalWF[usCh][usPW]) continue;
            if(WfType[usCh][usPW]==wSchroedinger && !ShortRangePotential[usCh][usPW]) continue;
            else if(WfType[usCh][usPW]==wExternal && !ExternalWF[usCh][usPW]) continue;
            else if(WfType[usCh][usPW]==wSquareWell && Sw_Pars[usCh][usPW][0]*Sw_Pars[usCh][usPW][1]) continue;
            TotNumSteps++;
        }
    }
    // unsigned TotalNumberOfBins = NumMomBins*NumCh*MaxNumPW;
    // TotNumSteps *= NumMomBins;

    // unsigned uMomBin;
    // unsigned short usCh;
    // unsigned short usPW;

    // double Time;
    // int pTotal;
    // int pTotalOld=0;
    // double Progress;
    // char* cdummy = new char [512];
    // double TotalSteps = TotNumSteps;
    // long long CurrentStep=0;
    // DLM_Timer dlmTimer;

    //the problem with the omp is that the same PotPar are used, we need separate instance for each thread if we want it to work
    //this is however difficult as currently I only pass a single pointer as PotPar, and I do not know how many arguments there are
    //the same problem should occur for the source
    unsigned short NumThreads = omp_get_num_procs();
    if(NumThreads>MaxNumThreads) NumThreads = MaxNumThreads;
    //#pragma omp parallel for private(uMomBin,usCh,usPW) num_threads(NumThreads)
    for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
        for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
            for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
        //compute to which MomBin, Polarization and PW corresponds this MPP,
        //i.e. map uMomBin, usCh and usPW to uMPP
        // uMomBin = uMPP/(NumCh*MaxNumPW);
        usCh = 0; // for pion-pion interacton ...? TODO !!
        usPW = 0; // for pion-pion interacton ...? TODO !!
        //since uMPP is build under the assumption all NumPW are the same
        //one has to check if usPW has a meaningful value!
        if(usPW>=NumPW[usCh]) continue;
        //if the potential for this channel and PW is zero => continue
        //if(!ShortRangePotential[usCh][usPW] && !ExternalWF[usCh][usPW]) continue;
        if(WfType[usCh][usPW]==wSchroedinger && !ShortRangePotential[usCh][usPW]) continue;
        else if(WfType[usCh][usPW]==wExternal && !ExternalWF[usCh][usPW]) continue;
        else if(WfType[usCh][usPW]==wSquareWell && Sw_Pars[usCh][usPW][0]*Sw_Pars[usCh][usPW][1]) continue;
        //skip momentum bins which have obtained an error code
        if(!MomBinConverged[uMomBin] && ExcludeFailedConvergence) continue;

        SavedWaveFunBins[uMomBin_out][uMomBin_side][uMomBin_long][usCh][usPW]=0;

        //if s+l is odd, than this partial wave will cancel out during the
        //symmetrization for identical particles
        if( IdenticalParticles && (usPW+Spin[usCh])%2 ) continue;

        double Momentum;
        //the momentum is taken from the center of the bin
        Momentum = GetTotMomentum3D(uMomBin);
        //Momentum = 0.5*(MomBin[uMomBin]+MomBin[uMomBin+1]);

        //perform the numerical computation
        if(WfType[usCh][usPW]==wSchroedinger){
            double* BufferWaveFunction;
            double* BufferRad;

            int q1q2 = Gamow?0:Q1Q2;
            unsigned NumComputedPoints = 2;//counting the initial two starting points
            //index that keeps track of which value in the arrays correspond to the previous, current and next computing step.
            short kOld;
            short kCurrent;
            short kNew;
            double WaveFun[3];
            double PosRad[3];
            double Rho[3];
            //value of the propagating function
            double PropFunVal[3];
            //the value of the prop. function without strong interaction.
            //in case it is equal to the PropFunVal, than the algorithm must have converged
            double PropFunWithoutSI[3];
            //Step size
            double DeltaRad[3];
            //Step size^2
            double DeltaRad2[3];

            double MaxDeltaRad;
            double MinDeltaRad;

            PropagatingFunction(PropFunWithoutSI[0], PropFunVal[0], StartRad, Momentum, usPW, usCh);
            MinDeltaRad = sqrt(fabs(EpsilonProp/(PropFunVal[0]+1e-64)));
            MaxDeltaRad = sqrt(EpsilonProp/(Momentum*Momentum));

            if(MinDeltaRad>MaxDeltaRad) MinDeltaRad=MaxDeltaRad;

            kOld=0;
            kCurrent=1;
            kNew=2;

            //set the initial step in r
            //!This definition of DeltaRad is used in other functions as well.
            //!It is vital that in case there is a need to change the definition, it is changed
            //!in all other functions as well!
            //DeltaRad[kOld] = RhoStep/Momentum;
            //DeltaRad2 = DeltaRad*DeltaRad;

            PosRad[kOld] = 0;
            Rho[kOld] = Momentum*PosRad[kOld];
            DeltaRad[kOld] = MinDeltaRad;
            DeltaRad2[kOld] = DeltaRad[kOld]*DeltaRad[kOld];
            PropagatingFunction(PropFunWithoutSI[kOld], PropFunVal[kOld], PosRad[kOld], Momentum, usPW, usCh);

            PosRad[kCurrent] = PosRad[kOld]+DeltaRad[kOld];
            Rho[kCurrent] = Momentum*PosRad[kCurrent];
            DeltaRad[kCurrent] = MinDeltaRad;
            DeltaRad2[kCurrent] = DeltaRad[kCurrent]*DeltaRad[kCurrent];
            PropagatingFunction(PropFunWithoutSI[kCurrent], PropFunVal[kCurrent], PosRad[kCurrent], Momentum, usPW, usCh);

            //the initial values for the wave function are set based on the solution without the strong potential.
            //this will of course lead to a wrong normalization in the asymptotic region, but this will be corrected for later on,
            //at this stage it is only important that the algorithm gets a meaningful guess so that we do not encounter
            //overflow problems during the calculation.
            WaveFun[kOld] = ReferencePartialWave(PosRad[kOld], Momentum, usPW, q1q2);
            WaveFun[kCurrent] = ReferencePartialWave(PosRad[kCurrent], Momentum, usPW, q1q2);

            bool Convergence = false;
            bool Converged = false;
            //at which point the convergence criteria first occurred
            double ConvergenceRadius=0;
            //how much after the convergence was the WF propagated
            //(this is needed for the normalization)
            const double ConvIntervalLength=3.14;
            const double ConvIntervalRadLength=ConvIntervalLength/Momentum;
            unsigned ConvIntervalSteps=0;

            unsigned StepOfMaxConvergedNumWF=0;
            double MaxConvergedNumWF=0;
            double MaxConvRho=0;
            //the last radius at which the computation was performed
            double MaxConvRad=0;
            double DeltaRadAtMaxConvPrev=DeltaRad[kOld];
            double DeltaRadAtMaxConv=DeltaRad[kCurrent];
            double ConvPointOldWeight=-1;
            double ConvPointWeight=0;

            const unsigned short MinConvSteps = 32;

            //the next two variables are estimates of the max. required bin numbers in case one fails to converge
            //MaxNumRadSteps assumes the worst case scenario -> we compute always with the worst possible step length
            //the second one is the optimistic -> we always compute with the best possible step length
            //bins that converge have usually less than MinNumRadSteps entries, bins that fail to converge a bit above.
            //EstNumRadSteps is a fair small over-estimation of the actual bins needed. It is computed by assuming we have the worst-case
            //scenario up until rho=1.57 or r=4 fm (whichever occurs first at this Momentum) and assumes that the rest as best case scenario.
            //This number looks to be c.a. 10x bigger than the actual num. bins, so use it for memory allocation.

            unsigned EstNumRadSteps = ceil(( (MaxRad*Momentum)>MaxRho?MaxRad-1.57/Momentum:MaxRho/Momentum+ConvIntervalRadLength-1.57/Momentum)/MaxDeltaRad)+MinConvSteps+1 +
                                        ceil(( 4.*FmToNu<1.57/Momentum?4.*FmToNu:1.57/Momentum )/MinDeltaRad);

            BufferWaveFunction = new double [EstNumRadSteps];
            BufferRad = new double [EstNumRadSteps];

            BufferWaveFunction[0] = WaveFun[kOld];
            BufferWaveFunction[1] = WaveFun[kCurrent];

            BufferRad[0] = PosRad[kOld];
            BufferRad[1] = PosRad[kCurrent];

            //!this is the main loop, which computes each next WF point until the result converges or
            //!a maximum value of rho is reached.
            //N.B. if the result is currently converging, the numerical method should not be interrupted even
            //if the MaxRad condition is met!
            while( (!Converged||(PosRad[kCurrent]<ConvergenceRadius+ConvIntervalRadLength))
                  && (PosRad[kOld]<MaxRad || Rho[kOld]<MaxRho || Converged) ){

                PosRad[kNew] = PosRad[kCurrent] + DeltaRad[kCurrent];
                Rho[kNew] = PosRad[kNew]*Momentum;

                WaveFun[kNew] = WaveFun[kCurrent]*(1.+DeltaRad[kCurrent]/DeltaRad[kOld]) -
                                WaveFun[kOld]*DeltaRad[kCurrent]/DeltaRad[kOld] +
                                PropFunVal[kCurrent]*WaveFun[kCurrent]*DeltaRad2[kCurrent];

                //if we run out of memory...
                if(NumComputedPoints>=EstNumRadSteps){
                    EstNumRadSteps*=2;

                    double* BufferTempWF = new double [EstNumRadSteps];
                    for(unsigned uTmp=0; uTmp<NumComputedPoints-1; uTmp++){
                        BufferTempWF[uTmp] = BufferWaveFunction[uTmp];
                    }
                    delete [] BufferWaveFunction;
                    BufferWaveFunction = BufferTempWF;

                    double* BufferTempRad = new double [EstNumRadSteps];
                    for(unsigned uTmp=0; uTmp<NumComputedPoints-1; uTmp++){
                        BufferTempRad[uTmp] = BufferRad[uTmp];
                    }
                    delete [] BufferRad;
                    BufferRad = BufferTempRad;
                }

                BufferWaveFunction[NumComputedPoints] = WaveFun[kNew];
if(NumComputedPoints<32 && NumComputedPoints%1==0 && uMomBin==0 && uMPP==0 && false){
printf("Momentum=%f\n",Momentum);
DEBUG=omp_get_thread_num();
}
else{
    //#pragma omp critical
    {
    DEBUG=-1;
    }
}
                PropagatingFunction(PropFunWithoutSI[kNew], PropFunVal[kNew], PosRad[kNew], Momentum, usPW, usCh);

                DeltaRad2[kNew] = EpsilonProp/(fabs(PropFunVal[kNew])+1e-64);
                DeltaRad[kNew] = sqrt(DeltaRad2[kNew]);
                if(DeltaRad[kNew]<MinDeltaRad){
                    DeltaRad[kNew]=MinDeltaRad;
                    DeltaRad2[kNew]=DeltaRad[kNew]*DeltaRad[kNew];
                }
                else if(DeltaRad[kNew]>MaxDeltaRad){
                    DeltaRad[kNew]=MaxDeltaRad;
                    DeltaRad2[kNew]=DeltaRad[kNew]*DeltaRad[kNew];
                }

                BufferRad[NumComputedPoints] = PosRad[kNew];
if(DEBUG==omp_get_thread_num()){
printf("Stuff: \n");
printf(" PropFunWithoutSI[%u]=%.15e\n",kNew,PropFunWithoutSI[kNew]);
printf(" PropFunVal[%u]=%.15e\n",kNew,PropFunVal[kNew]);
printf(" BufferRad[%u]=%.15e\n",NumComputedPoints,BufferRad[NumComputedPoints]);
DEBUG=-1;
}

                Convergence =   fabs( (PropFunWithoutSI[kOld]-PropFunVal[kOld])/(fabs(PropFunWithoutSI[kOld]+PropFunVal[kOld])+1e-64) ) < EpsilonConv &&
                                fabs( (PropFunWithoutSI[kCurrent]-PropFunVal[kCurrent])/(fabs(PropFunWithoutSI[kCurrent]+PropFunVal[kCurrent])+1e-64) ) < EpsilonConv &&
                                fabs( (PropFunWithoutSI[kNew]-PropFunVal[kNew])/(fabs(PropFunWithoutSI[kNew]+PropFunVal[kNew])+1e-64) ) < EpsilonConv;

                //in case we have detected a Convergence
                //1) make sure the convergence is real and not a local artifact of the potential
                //2) reset all variables used to characterize the convergence region
                if( Convergence && !Converged){
                    Converged = true;
                    ConvergenceRadius = PosRad[kNew];

                    //2):
                    MaxConvergedNumWF = 0;
                    MaxConvRho = 0;
                    MaxConvRad = 0;
                    StepOfMaxConvergedNumWF = 0;
                    DeltaRadAtMaxConvPrev = MinDeltaRad;
                    DeltaRadAtMaxConv = MinDeltaRad;
                }
                //this condition makes sure that 1) is fulfilled.
                else{
                    Converged = Convergence;
                }

                //this will reset the ConvIntervalSteps in case the convergence is reset
                ConvIntervalSteps *= Converged;
                //this will count the number of steps computed after convergence
                ConvIntervalSteps += Converged;

                ConvPointWeight = (PosRad[kNew]-ConvergenceRadius)/ConvIntervalRadLength;

                if(fabs(MaxConvergedNumWF)*ConvPointOldWeight<=fabs(WaveFun[kNew])*ConvPointWeight){
                    MaxConvergedNumWF = WaveFun[kNew];
                    MaxConvRho = Rho[kNew];
                    MaxConvRad = PosRad[kNew];
                    StepOfMaxConvergedNumWF = NumComputedPoints;
                    DeltaRadAtMaxConvPrev = DeltaRad[kCurrent];
                    DeltaRadAtMaxConv = DeltaRad[kNew];
                    ConvPointOldWeight = ConvPointWeight;
                }

                NumComputedPoints++;

                kNew++;
                kNew=kNew%3;
                kCurrent++;
                kCurrent=kCurrent%3;
                kOld++;
                kOld=kOld%3;
            }//while(!Converged && PosRad[kOld]<MaxRad)

            //this is not desired for the later computation
            if(StepOfMaxConvergedNumWF==NumComputedPoints-1){
                StepOfMaxConvergedNumWF--;
                MaxConvergedNumWF = BufferWaveFunction[StepOfMaxConvergedNumWF];
                MaxConvRho -= DeltaRadAtMaxConvPrev*Momentum;
                MaxConvRad -= DeltaRadAtMaxConvPrev;
                DeltaRadAtMaxConv = DeltaRadAtMaxConvPrev;
            }
    //! SOMETIMES AT HIGH MOMENTA THE RESULT SIMPLY DOES NOT CONVERGE TO THE REQUIRED VALUE
    //за това дай на юзера опцията да реши какво да прави с неконвергирали резултати
    //1) изхвърли ги
    //2) запаза ги, като нормировката я направи на база максимума в MaxRho-3.14 (btw. make 3.14 the min. value for MaxRho!!!)


            //if the maximum wave-function is zero the computation will fail!
            //By design this should not really happen.
            if(!MaxConvergedNumWF && Notifications>=nWarning){
                printf("\033[1;33mWARNING:\033[0m MaxConvergedNumWF is zero, which is not allowed and points to a bug in the code!\n");
                printf("         Please contact the developers and do not trust your current results!\n");
            }

            bool NR_Status;

            //!Now follows the normalization of the numerical wave function to the asymptotic solution
            //the individual steps are explained in detail in the official CATS documentation
            if(MomBinConverged[uMomBin] || !ExcludeFailedConvergence){
                double ShiftRadStepLen = 0.1/Momentum;
                const unsigned MaximumShiftIter = 121;

                double ShiftRad;

                double DownShift=0;
                double UpShift=ShiftRadStepLen;
                if((!BufferWaveFunction[StepOfMaxConvergedNumWF] || !BufferWaveFunction[StepOfMaxConvergedNumWF+1])&&Notifications>=nWarning){
                    printf("\033[1;33mWARNING:\033[0m BufferWaveFunction is zero, which is not allowed and points to a bug in the code!\n");
                    printf("         Please contact the developers and do not trust your current results!\n");
                }
                double NumRatio = BufferWaveFunction[StepOfMaxConvergedNumWF+1]/BufferWaveFunction[StepOfMaxConvergedNumWF];
                double DownValue;
                double UpValue;
                double SignProduct;
                double SignRefAtLimits;

                //CHECK IF ZERO IS A VIABLE OPTION!
                DownShift = -0.5*ShiftRadStepLen;
                if(fabs(DownShift)>MaxConvRad) DownShift = -MaxConvRad;
                UpShift = 0.5*ShiftRadStepLen;

                //a potential solution should be locked in a region where the ShiftedReferenceWave*NumericalSol have the same sign.
                //furthermore the ShiftedReferenceWave cannot possibly change sign in this region
                for(unsigned uShift=0; uShift<MaximumShiftIter; uShift++){
                    //positive side
                    CurrentRhoStep = DeltaRadAtMaxConv*Momentum;
                    DownValue = AsymptoticRatio(MaxConvRad+DownShift, Momentum, usPW, q1q2) - NumRatio;
                    UpValue = AsymptoticRatio(MaxConvRad+UpShift, Momentum, usPW, q1q2) - NumRatio;

                    if(DownValue*UpValue<0){
                        ShiftRad = 0.5*(DownShift+UpShift);
                        SignProduct = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW, q1q2)*MaxConvergedNumWF;
                        SignRefAtLimits = ReferencePartialWave(MaxConvRad+DownShift, Momentum, usPW, q1q2)*ReferencePartialWave(MaxConvRad+UpShift, Momentum, usPW, q1q2);
                        if(SignProduct>0 && SignRefAtLimits>0) break;
                    }

                    //negative side
                    if(UpShift<=MaxConvRad){
                        DownValue = AsymptoticRatio(MaxConvRad-DownShift, Momentum, usPW, q1q2) - NumRatio;
                        UpValue = AsymptoticRatio(MaxConvRad-UpShift, Momentum, usPW, q1q2) - NumRatio;
                        if(DownValue*UpValue<0){
                            ShiftRad = -0.5*(DownShift+UpShift);
                            SignProduct = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW, q1q2)*MaxConvergedNumWF;
                            SignRefAtLimits = ReferencePartialWave(MaxConvRad-DownShift, Momentum, usPW, q1q2)*ReferencePartialWave(MaxConvRad-UpShift, Momentum, usPW, q1q2);
                            if(SignProduct>0 && SignRefAtLimits>0){
                                double Temp = DownShift;
                                DownShift = -UpShift;
                                UpShift = -Temp;
                                break;
                            }
                        }
                    }

                    //the factor 0.95 is there just to make sure that we don't by accident
                    //miss a zero just around the limiting values
                    UpShift += 0.95*ShiftRadStepLen;
                    DownShift = UpShift-ShiftRadStepLen;
                }

                ShiftRad = NewtonRapson(&CATS::AsymptoticRatio,
                        DeltaRadAtMaxConv, usPW, Momentum, q1q2, MaxConvRad+DownShift, MaxConvRad+UpShift, NumRatio,NR_Status) - MaxConvRad;
                //if(!NR_Status&&Notifications>=nWarning) {printf("The above warning triggered for bin nr. %u (k=%.2f)\n",uMomBin,GetMomentum(uMomBin));}
                if(!NR_Status) {MomBinConverged[uMomBin]=false;}

                double ShiftRho = ShiftRad*Momentum;
                double Norm = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW, q1q2)/MaxConvergedNumWF;

                PhaseShift[uMomBin][usCh][usPW] = ShiftRho;
                PhaseShiftF[usCh][usPW][uMomBin] = ShiftRho;

                //we save all entries up to the point in which the normalization was performed. For higher values
                //we will later on use the asymptotic form. N.B. in principle one could save a bit of space and
                //save the result at the first moment of detected convergence, however than in case the convergence
                //was not "perfect" there might be some inaccuracy in the interval up to the normalization point.
                //btw. the step size is set to be the MaxDeltaRad

                SavedWaveFunBins[uMomBin_out][uMomBin_side][uMomBin_long][usCh][usPW] = StepOfMaxConvergedNumWF+1;
                unsigned& SWFB = SavedWaveFunBins[uMomBin_out][uMomBin_side][uMomBin_long][usCh][usPW];

                if(WaveFunRad[uMomBin][usCh][usPW]) delete [] WaveFunRad[uMomBin][usCh][usPW];
                WaveFunRad[uMomBin][usCh][usPW] = new double [SWFB+1];

                if(WaveFunctionU[uMomBin][usCh][usPW]) delete [] WaveFunctionU[uMomBin][usCh][usPW];
                WaveFunctionU[uMomBin][usCh][usPW] = new complex<double> [SWFB];

                for(unsigned uPoint=0; uPoint<SWFB; uPoint++){
                    WaveFunctionU[uMomBin][usCh][usPW][uPoint] = Norm*CPF[uMomBin]*BufferWaveFunction[uPoint];
                }

                //we set up those as bin-range (i.e. the middle of the bin should be buffer rad)
                //we start from r=0, than proceed by adding the limit between different bins as the mid-way to the next point
                WaveFunRad[uMomBin][usCh][usPW][0] = 0;
                for(unsigned uPoint=1; uPoint<SWFB; uPoint++){
                    WaveFunRad[uMomBin][usCh][usPW][uPoint] = (BufferRad[uPoint]+BufferRad[uPoint-1])*0.5;
                }
                //the very last point we simply set as the double distance between the last bin-limit and the last radius value
                WaveFunRad[uMomBin][usCh][usPW][SWFB] = 2*BufferRad[SWFB-1]-WaveFunRad[uMomBin][usCh][usPW][SWFB-1];

            }//if(MomBinConverged[uMomBin] || !ExcludeFailedConvergence)

            //in case the method failed to converge, the whole momentum bin is marked as unreliable and no
            //further computations are done. This point will be excluded from the final result
            if(!Converged){
                MomBinConverged[uMomBin] = false;
                if(Notifications>=nWarning)
                    printf("          \033[1;33mWARNING:\033[0m A momentum bin at %.2f MeV failed to converge!\n", Momentum);
                if(ExcludeFailedConvergence && Notifications>=nWarning){
                    printf("                   It will be excluded from the final result!\n");
                }
            }
            if(!NR_Status){
                MomBinConverged[uMomBin] = false;
                if(Notifications>=nWarning)
                    printf("          \033[1;33mWARNING:\033[0m The failed NewtonRapson occurred at %.2f MeV!\n", Momentum);
                if(ExcludeFailedConvergence && Notifications>=nWarning){
                    printf("                   The corresponding bin will be excluded from the final result!\n");
                }
            }

            delete [] BufferRad;
            delete [] BufferWaveFunction;
        }//end of the numerical computation
        //the case with external wave function
        else if(WfType[usCh][usPW]==wExternal){
            SavedWaveFunBins[uMomBin_out][uMomBin_side][uMomBin_long][usCh][usPW] = ExternalWF[usCh][usPW]->GetNbins(1);
            unsigned& SWFB = SavedWaveFunBins[uMomBin_out][uMomBin_side][uMomBin_long][usCh][usPW];

            if(WaveFunRad[uMomBin][usCh][usPW]) delete [] WaveFunRad[uMomBin][usCh][usPW];
            WaveFunRad[uMomBin][usCh][usPW] = new double [SWFB+1];

            if(WaveFunctionU[uMomBin][usCh][usPW]) delete [] WaveFunctionU[uMomBin][usCh][usPW];
            WaveFunctionU[uMomBin][usCh][usPW] = new complex<double> [SWFB];

            for(unsigned uPoint=0; uPoint<=SWFB; uPoint++){
                //the input from outside is supposed to be in fermi, hence to conversion
                WaveFunRad[uMomBin][usCh][usPW][uPoint] = ExternalWF[usCh][usPW]->GetBinLowEdge(1,uPoint)*FmToNu;
            }

            double EvalPoint[2];
            EvalPoint[0] = GetMomentum(uMomBin);
            for(unsigned uPoint=0; uPoint<SWFB; uPoint++){
                EvalPoint[1] = ExternalWF[usCh][usPW]->GetBinCenter(1,uPoint);
                WaveFunctionU[uMomBin][usCh][usPW][uPoint] = CPF[uMomBin]*ExternalWF[usCh][usPW]->Eval(EvalPoint)*FmToNu;
            }
        }
        else if(WfType[usCh][usPW]==wSquareWell){
          double& sw_depth = Sw_Pars[usCh][usPW][0];
          double& sw_width = Sw_Pars[usCh][usPW][1];
          double DeltaRad = EpsilonProp*50.;
          SavedWaveFunBins[uMomBin_out][uMomBin_side][uMomBin_long][usCh][usPW] = ceil(sw_width/DeltaRad);
          unsigned& SWFB = SavedWaveFunBins[uMomBin_out][uMomBin_side][uMomBin_long][usCh][usPW];

          if(WaveFunRad[uMomBin][usCh][usPW]) delete [] WaveFunRad[uMomBin][usCh][usPW];
          WaveFunRad[uMomBin][usCh][usPW] = new double [SWFB+1];

          if(WaveFunctionU[uMomBin][usCh][usPW]) delete [] WaveFunctionU[uMomBin][usCh][usPW];
          WaveFunctionU[uMomBin][usCh][usPW] = new complex<double> [SWFB];

          for(unsigned uPoint=0; uPoint<=SWFB; uPoint++){
              WaveFunRad[uMomBin][usCh][usPW][uPoint] = uPoint*DeltaRad;
          }

          double EvalPoint[2];
          EvalPoint[0] = GetMomentum(uMomBin);
          for(unsigned uPoint=0; uPoint<SWFB; uPoint++){
//WORK IN PROGRESS !! ??
              WaveFunctionU[uMomBin][usCh][usPW][uPoint] = CPF[uMomBin];
          }
        }
        else{
          if(Notifications>=nError)
              printf("          \033[1;31mERROR:\033[0m Major bug in ComputeWaveFunction \033[0m contact the developers!\n");
          return;
        }

        //#pragma omp atomic
        CurrentStep++;
        //#pragma omp critical
        {
            Progress = double(CurrentStep)/TotalSteps;
            pTotal = int(Progress*100);
            if(pTotal!=pTotalOld){
                Time = double(dlmTimer.Stop())/1e6;
                Time = round((1./Progress-1.)*Time);
                ShowTime((long long)(Time), cdummy, 2, true, 5);
                if(Notifications>=nAll) printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
                cout << flush;
                pTotalOld = pTotal;
            }
        }

                }
        }
    }//for(unsigned uMPP=0; uMPP<TotalNumberOfBins; uMPP++)
    ComputedWaveFunction = true;
    if(Notifications>=nAll) printf("\r\033[K");
    delete [] cdummy;
}

void CATS3D::ComputeTotWaveFunction(const bool& ReallocateTotWaveFun){
    if(!BaseSourceGrid) return;

    if(WaveFunction2 && ReallocateTotWaveFun){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                delete [] WaveFunction2[uMomBin][uGrid];
            }
            delete [] WaveFunction2[uMomBin];
        }
        delete [] WaveFunction2; WaveFunction2 = NULL;
    }

    NumGridPts = BaseSourceGrid->GetNumEndNodes();

    if(!WaveFunction2){
        WaveFunction2 = new double** [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            WaveFunction2[uMomBin] = new double* [NumGridPts];
            for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                WaveFunction2[uMomBin][uGrid] = new double [NumCh];
            }
        }
    }

    //double Momentum;
    double Radius;
    double CosTheta;

    double Time;
    int pTotal;
    int pTotalOld=0;
    double Progress;
    char* cdummy = new char [512];
    double TotalSteps = double(NumMomBins)*double(NumGridPts);
    long long CurrentStep=0;
    DLM_Timer dlmTimer;
    unsigned short NumThreads = omp_get_num_procs();
    if(NumThreads>MaxNumThreads) NumThreads = MaxNumThreads;
    //#pragma omp parallel for private(Radius,CosTheta) num_threads(NumThreads)
    for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
        for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
            for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                //Momentum = GetMomentum(uMomBin);
                for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                    // if(PoorManRenorm!=1 && UseAnalyticSource==false) Radius = BaseSourceGrid->GetParValue(uGrid, 0)*FmToNu*PoorManRenorm;
                    // else 

                    // Here is the most important part in 3-D femto CATS !!! TODO!! 
                    Radius = std::sqrt(std::pow(BaseSourceGrid->GetParValue(uGrid, 0),2)+ 
                                       std::pow(BaseSourceGrid->GetParValue(uGrid, 1),2)+
                                       std::pow(BaseSourceGrid->GetParValue(uGrid, 2),2))*FmToNu;         
                    CosTheta = ( MomBinCenter_out[uMomBin_out]* BaseSourceGrid->GetParValue(uGrid, 0) + MomBinCenter_side[uMomBin_side]* BaseSourceGrid->GetParValue(uGrid, 1) + MomBinCenter_long[uMomBin_long]* BaseSourceGrid->GetParValue(uGrid, 2) ) / 
                                (std::sqrt(std::pow(BaseSourceGrid->GetParValue(uGrid, 0),2)+ 
                                       std::pow(BaseSourceGrid->GetParValue(uGrid, 1),2)+
                                       std::pow(BaseSourceGrid->GetParValue(uGrid, 2),2)) + 
                                 std::sqrt(std::pow(MomBinCenter_out[uMomBin_out],2)+ 
                                       std::pow(MomBinCenter_side[uMomBin_side],2)+
                                       std::pow(MomBinCenter_long[uMomBin_long],2)));
        
        
                    for(unsigned usCh=0; usCh<NumCh; usCh++){
                        WaveFunction2[uMomBin][uGrid][usCh] =   ThetaDependentSource?
                                                                EffectiveFunctionTheta(uMomBin, Radius, CosTheta, usCh):
                                                                EffectiveFunction(uMomBin, Radius, usCh);
                    }
        
                    //#pragma omp atomic
                    CurrentStep++;
                    //#pragma omp critical
                    {
                        Progress = double(CurrentStep)/TotalSteps;
                        pTotal = int(Progress*100);
                        if(pTotal!=pTotalOld){
                            Time = double(dlmTimer.Stop())/1e6;
                            Time = round((1./Progress-1.)*Time);
                            ShowTime((long long)(Time), cdummy, 2, true, 5);
                            if(Notifications>=nAll) printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
                            cout << flush;
                            pTotalOld = pTotal;
                        }
                    }
        
                }
            }
        }
    }
    if(Notifications>=nAll) printf("\r\033[K");
    delete [] cdummy;
}

void CATS3D::SetUpSourceGrid(){
    double LIMIT = GridEpsilon;
    if(!LIMIT){
        LIMIT=1./1024.;
    }
    short MAXDEPTH = GridMaxDepth;
    if(!MAXDEPTH){
        MAXDEPTH=14;
    }
    short DIM = 3;
    unsigned MAXGRIDPTS = uipow(2,MAXDEPTH*DIM);

    double* MEAN = new double[DIM];
    double* LENGTH = new double[DIM];

    MEAN[0] = SourceMaxRad*0.5*NuToFm; // out / side / long to clarify ....
    LENGTH[0] = SourceMaxRad*NuToFm;

    MEAN[1] = SourceMaxRad*0.5*NuToFm;
    LENGTH[1] = SourceMaxRad*NuToFm;

    MEAN[2] = SourceMaxRad*0.5*NuToFm;
    LENGTH[2] = SourceMaxRad*NuToFm;
    
    if(BaseSourceGrid){
        delete BaseSourceGrid; BaseSourceGrid=NULL;
    }
    if(UseAnalyticSource){
        //for setting up the grid we take the "mean" value of q_out, q_side, q_long
        AnaSourcePar->SetVariable(0,(MomBin_out[0]+MomBin_out[NumMomBins_out])*0.5,false);
        AnaSourcePar->SetVariable(1,(MomBin_side[0]+MomBin_side[NumMomBins_side])*0.5,false);
        AnaSourcePar->SetVariable(2,(MomBin_long[0]+MomBin_long[NumMomBins_long])*0.5,false);

//printf(" AnaSourcePar[3]=%e\n",AnaSourcePar->GetParameter(0));
        BaseSourceGrid = new CATSelder(DIM, GridMinDepth, MAXDEPTH, LIMIT, MEAN, LENGTH,
                         this, AnaSourcePar, NULL, 0, AutoNormSource);
    }
    else if(NumPairs){
        //sorts the whole data according to the GridBoxId (all bins!)
        //this will setup the base grid
        SortAllData();
        BaseSourceGrid = new CATSelder(DIM, GridMinDepth, MAXDEPTH, LIMIT, MEAN, LENGTH,
                                        NULL, NULL, GridBoxId, NumPairs, AutoNormSource);
    }
    else{
        BaseSourceGrid = NULL;
    }

    int64_t ArrayPosition=0;

    /// To be implemented!!!!! TODO !!! To understand and to do !!
    // //sorts the data in momentum blocks and sorts them according to the GridBoxId
    // for(unsigned uPair=0; uPair<NumPairs; uPair++){
    //     GridBoxId[uPair] += int64_t(MAXGRIDPTS)*int64_t(PairMomBin[uPair]);
    // }
    // SortAllData();
    // for(unsigned uPair=0; uPair<NumPairs; uPair++){
    //     GridBoxId[uPair] -= (int64_t(GridBoxId[uPair])/int64_t(MAXGRIDPTS))*int64_t(MAXGRIDPTS);
    // } 

    //Resets kSourceGrid3D
    if(kSourceGrid3D){
        for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
            for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
                for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                    if(kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]) delete kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long];
                }
            }
        }
        delete [] kSourceGrid3D;
    }

    kSourceGrid3D = new CATSelder* [NumMomBins_out][NumMomBins_side][NumMomBins_long];
    for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
        for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
            for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]=NULL;
            }
        }
    }//sets up kSourceGrid3D

    for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
        for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
            for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                
                // if(!MomDepSource&&uMomBin>0){
                //     break;
                // } // Aim ? to be implement ...  TODO !!

                AnaSourcePar->SetVariable(0,MomBinCenter_out[uMomBin_out],false);
                AnaSourcePar->SetVariable(1,MomBinCenter_side[uMomBin_side],false);
                AnaSourcePar->SetVariable(2,MomBinCenter_long[uMomBin_long],false);            
                
                if(UseAnalyticSource){
                    kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long] = new CATSelder(BaseSourceGrid, this, AnaSourcePar, NULL, 0, AutoNormSource);
                }
                // else if(LoadedPairsPerMomBin[uMomBin_out][uMomBin_side][uMomBin_long]){
                //     kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long] = new CATSelder(BaseSourceGrid, NULL, NULL,
                //                                         &GridBoxId[ArrayPosition], LoadedPairsPerMomBin[uMomBin_out][uMomBin_side][uMomBin_long], AutoNormSource);
                //     ArrayPosition += LoadedPairsPerMomBin[uMomBin_out][uMomBin_side][uMomBin_long];
                // }
                else{
                    kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long] = NULL;
                }
            }
        }
    }

    if(kbSourceGrid){
        for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
            for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
                for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){         if(!kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long]) break;
                    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                        if(kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]) delete kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin];
                    }
                    delete [] kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long];
                }
            }
        }
        delete [] kbSourceGrid; kbSourceGrid=NULL;
    }

    //sorts the data in k,b blocks and sorts them according to the GridBoxId
    //this is needed only if we use a data source with and more than 1 b-bin
    if(!UseAnalyticSource && NumIpBins>1){
        ArrayPosition = 0;

        //sorts the data in momentum-b blocks and sorts them according to the GridBoxId
        // TODO !!!
        // for(unsigned uPair=0; uPair<NumPairs; uPair++){
        //     GridBoxId[uPair] += int64_t(MAXGRIDPTS)*int64_t(PairIpBin[uPair])+
        //                         int64_t(MAXGRIDPTS)*int64_t(NumIpBins)*int64_t(PairMomBin[uPair]);
        // }
        // SortAllData();
        // for(unsigned uPair=0; uPair<NumPairs; uPair++){
        //     GridBoxId[uPair] -= (int64_t(GridBoxId[uPair])/int64_t(MAXGRIDPTS))*int64_t(MAXGRIDPTS);
        // }

        kbSourceGrid3D = new CATSelder* [NumMomBins_out][NumMomBins_side][NumMomBins_long];
        for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
            for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
                for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                    kbSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]=NULL;
                }
            }
        }

        for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
            for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
                for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                    // if(!MomDepSource&&uMomBin>0){
                    //     break;
                    // }  // TODO !!
                    kbSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long] = new CATSelder* [NumIpBins];
                    AnaSourcePar->SetVariable(0,GetMomentum(uMomBin),false);
                    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                        if(UseAnalyticSource){
                            kbSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = new CATSelder(BaseSourceGrid, this, AnaSourcePar, NULL, 0, AutoNormSource);
                        }
                        // else if(LoadedPairsPerBin[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]){
                        //     kbSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = new CATSelder(BaseSourceGrid, NULL, NULL,
                        //                                                 &GridBoxId[ArrayPosition], LoadedPairsPerBin[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin], AutoNormSource);
                        //     ArrayPosition += LoadedPairsPerBin[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin];
                        // }
                        else{
                            kbSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = NULL;
                        }
                    }
                }
            }
        }
    }
    SourceGridReady = true;
    SourceUpdated = true;

    delete [] MEAN;
    delete [] LENGTH;
}

void CATS3D::UpdateSourceGrid(){
    if(BaseSourceGrid) BaseSourceGrid->Update();
    for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
        for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
            for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                if(kSourceGrid3D && kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]) kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->Update();
            }
        }
    }

    if(kbSourceGrid3D){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                if(kbSourceGrid3D && kbSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long] && kbSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin])
                    kbSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]->Update();
            }
        }
    }
    SourceUpdated = true;
}

void CATS3D::FoldSourceAndWF(){
    if(!kbCorrFun){
        kbCorrFun = new double* [NumMomBins_out][NumMomBins_side][NumMomBins_long];
        kbCorrFunErr = new double* [NumMomBins_out][NumMomBins_side][NumMomBins_long];
        kCorrFun = new double [NumMomBins_out][NumMomBins_side][NumMomBins_long];
        kCorrFunErr = new double [NumMomBins_out][NumMomBins_side][NumMomBins_long];
        for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
            for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
                for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                    kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] = new double [NumIpBins];
                    kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = new double [NumIpBins];
                    kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] = 0;
                    kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = 0;
                    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                        kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = 0;
                        kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = 0;
                    }
                }
            }
        }
    }
    else if(!ComputedCorrFunction){
        for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
            for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
                for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){
                    kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] = 0;
                    kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = 0;
                    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                        kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = 0;
                        kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = 0;
                    }
                }
            }
        }
    }
    else return;

    //this can happen if we change something about the source, but force CATS to use the same grid.
    //in such a case before folding the date one needs to update the values for the source!
    if(!SourceUpdated) UpdateSourceGrid();

    unsigned NumGridPts;
    double SourceVal;
    double GridPointSize;
    double WaveFunVal;
    double Integrand;
    double SourceInt;
    double SourceIntCut;
    double LostInSofia;
    unsigned short NumThreads = omp_get_num_procs();
    if(NumThreads>MaxNumThreads) NumThreads = MaxNumThreads;
    unsigned uMomBinSource;
    //#pragma omp parallel for private(NumGridPts,SourceVal,WaveFunVal,Integrand,SourceInt,SourceIntCut) num_threads(NumThreads)
    for(unsigned uMomBin_out=0; uMomBin_out<NumMomBins_out; uMomBin_out++){
        for(unsigned uMomBin_side=0; uMomBin_side<NumMomBins_side; uMomBin_side++){
            for(unsigned uMomBin_long=0; uMomBin_long<NumMomBins_long; uMomBin_long++){        
                uMomBinSource = MomDepSource?uMomBin:0;
                //printf("uMomBinSource = %u\n",uMomBinSource);
                kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] = 0;
                kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = 0;
                //N.B. this whole for-looped was stopped to be maintained, there could be issues with 
                //the normalization of the source and the correlation
                //either way, this does NOT execute at all for Analytic sources
                for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                    //perform the k,b analysis only for a data-source which has at least 2 b-bins
                    if(UseAnalyticSource || NumIpBins<=1) continue;

                    if(!kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]){
                        kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = 0;
                        kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = 1e6;
                        continue;
                    }
                    kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = 0;
                    kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = 0;

                    NumGridPts = kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]->GetNumEndNodes();
                    SourceInt = 0;
                    SourceIntCut = 0;
                    uMomBinSource = MomDepSource?uMomBin:0;
                    for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                        double Radius = kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]->GetParValue(uGrid, 0)*FmToNu;
                        //this should return zero in case of !ThetaDependentSource,
                        //maybe worth doing a QA to make sure
                        //CosTheta = kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]->GetParValue(uGrid, 1);
                        SourceVal = kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]->GetGridValue(uGrid);
                        if(!SourceVal) continue;
                        SourceInt += SourceVal;
                        if(Radius<SourceMinRad || Radius>SourceMaxRad) {continue;}
                        SourceIntCut += SourceVal;
                        WaveFunVal=0;
                        int uTotMomBin = ???// TODO !!! q_total = sqrt(q_o^2 + q_s^2 + q_l^2)
                        for(unsigned usCh=0; usCh<NumCh; usCh++) WaveFunVal+=(WaveFunction2[uTotMomBin][uGrid][usCh])*ChannelWeight[usCh];
                        Integrand = SourceVal*WaveFunVal;
                        kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] += Integrand;
                        kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] += pow(Integrand*kbSourceGrid[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]->GetGridError(uGrid),2);
                    }//uGrid
            
                    //ideally this should not happen, but in case it does (i.e. the source is normalized to value above 1)
                    //here this is corrected for
                    if(SourceInt>1+1e-4){
                        if(Notifications>=nWarning)
                            printf("\033[1;33mWARNING:\033[0m The source had to be renormalized (Integrated value=%.4e), please check the validity of your input function!\n",SourceInt);
                        //if(AnaSourcePar)
                            //printf("AnaSourcePar: p0=%e; p1=%e",AnaSourcePar->GetParameter(0),AnaSourcePar->GetParameter(1));
                        //if(ForwardedSourcePar)
                            //printf("ForwardedSourcePar: p0=%e; p1=%e",ForwardedSourcePar->GetParameter(0),ForwardedSourcePar->GetParameter(1));
                        kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] /= SourceInt;
                        kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = sqrt(kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin])/SourceInt;
                        SourceIntCut/=SourceInt;
                        SourceInt=1.;
                    }
                    //printf("\033[1;36mDEBUG (1):\033[0m Integrated value=%.4e\n",SourceInt);
                    //if the source is assumed normalized, any deviation below unity is considered to be an effect of a
                    //long-range tail that is unaccounted for due to the numerical precision. This tail should result in
                    //a flat residual that is added to the correlation
                    if(NormalizedSource && SourceIntCut<1.){
                        //this is an additional scale factor, evaluated based on the weight of channels with non-zero asymptotic solution
                        double AsympWeight = 0;
                        for(unsigned short usCh=0; usCh<NumCh; usCh++){
                          AsympWeight += ChannelWeight[usCh]*(!OnlyNumPw[usCh]);
                        }
                        kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] += (1.-double(SourceIntCut))*AsympWeight;
                        kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = sqrt(kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin]);
                    }
                    //else simply re-normalize the source to the total value of the integral (even if below 1)
                    //I am not sure if this even makes sense, consider in the future to use only the NormalizedSource case
                    else{
                        kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] *= double(SourceInt)/double(SourceIntCut);
                        kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin] = sqrt(kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin])*double(SourceInt)/double(SourceIntCut);
                    }
        
                    //in case we use event mixing (data-source have more than a single b-bin), this is how we proceed
                    if(MixingDepth>1){
                        kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] += WeightIp[uIpBin]*kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin];
                        kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] += pow(WeightIpError[uIpBin]*kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin],2) +
                                                            pow(WeightIpError[uIpBin]*kbCorrFun[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin],2) +
                                                            pow(WeightIp[uIpBin]*kbCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long][uIpBin],2);
                    }
                    //printf("SourceInt = %f\n",SourceInt);
                }


                //the thing of interest
                if( UseAnalyticSource || NumIpBins<=1 || MixingDepth==1 ){
                    if(!kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]){
                        kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = 1e6;
                    }
                    //this is where I will always end up with an analytic source
                    else{
                        NumGridPts = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetNumEndNodes();
                        SourceInt = 0;
                        SourceIntCut = 0;
                        LostInSofia = 0;
                        //printf("\033[1;36m DEBUG (2):\033[0m NumGridPts=%u\n",NumGridPts);
                        for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                            double Radius = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetParValue(uGrid, 0)*FmToNu;
                            double NextRadius;
                            if(uGrid<NumGridPts-1)
                                NextRadius = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetParValue(uGrid+1, 0)*FmToNu;
                            else NextRadius = 1e16;
                            //this should return zero in case of !ThetaDependentSource,
                            //maybe worth doing a QA to make sure
                            //if(uMomBinSource==0&&uGrid==0){
                            //double rad = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetParValue(uGrid, 0);
                            //double cth = ThetaDependentSource?kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetParValue(uGrid, 1):-100;
                            //printf("%u r = %.3f; cth = %.5f\n",uGrid,rad,cth);
                            //usleep(100e3);
                            //}
                            SourceVal = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetGridValue(uGrid);
                            GridPointSize = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetGridSize(uGrid);
                            if(use_sofia_model && Radius<Sofia_Cutoff*FmToNu){
                                LostInSofia += SourceVal;
                                SourceVal = 0;
                            }
                            
                            if(compute_true_source){
                                double source_box_rad_0 = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetParValue(uGrid, 0);
                                double source_box_rad_1 = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetParValue(uGrid+1, 0);
                                double source_box_rad = (source_box_rad_1+source_box_rad_0)*0.5;
                                double source_box_cth_0 = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetParValue(uGrid, 1);
                                double source_box_cth_1 = kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetParValue(uGrid+1, 1);  
                                double source_box_cth = (source_box_cth_1+source_box_cth_0)*0.5;
                                //if(uGrid==NumGridPts-1)              
                                //printf("%f %f %f %f\n",source_box_rad_0,source_box_rad_1,source_box_cth_0,source_box_cth_1);
                                //usleep(100e3);
                                unsigned which_ts_bin_rad = TheTrueSource->GetBin(1, source_box_rad);
                                unsigned which_ts_bin_cth = TheTrueSource->GetBin(2, source_box_cth);
                                unsigned which_bin[3];
                                which_bin[0] = uMomBin;
                                which_bin[1] = which_ts_bin_rad;
                                which_bin[2] = which_ts_bin_cth;
                                
                                //printf("%u(%.2f) %u(%.2f) %u(%.2f) %f\n",uMomBin,MomBinCenter[uMomBin_out][uMomBin_side][uMomBin_long],which_ts_bin_rad,source_box_rad,which_ts_bin_cth,source_box_cth,SourceVal);
                                //usleep(50e3);
        
                                if(         ThetaDependentSource &&
                                            source_box_rad_0<TheTrueSource->GetBinCenter(1,which_ts_bin_rad) && 
                                            source_box_rad_1>TheTrueSource->GetBinCenter(1,which_ts_bin_rad) &&
                                            source_box_cth_0<TheTrueSource->GetBinCenter(2,which_ts_bin_rad) &&
                                            source_box_cth_1<TheTrueSource->GetBinCenter(2,which_ts_bin_rad)){
                                        TheTrueSource->SetBinContent(which_bin, SourceVal/GridPointSize);
                                }
                                else if(    !ThetaDependentSource &&
                                            source_box_rad_0<TheTrueSource->GetBinCenter(1,which_ts_bin_rad) && 
                                            source_box_rad_1>TheTrueSource->GetBinCenter(1,which_ts_bin_rad)){
                                    for(unsigned uCth=0; uCth<TheTrueSource->GetNbins(2); uCth++){
                                        which_bin[2] = uCth;
                                        TheTrueSource->SetBinContent(which_bin, SourceVal/GridPointSize);
                                    }
                                }
                                
                            }
                            //if(use_sofia_model){
                            //printf("\033[1;36m  DEBUG (2):\033[0m SourceVal=%.4e\n",SourceVal);
                            //}
                            if(!SourceVal) continue;
                            SourceInt += SourceVal;
                            if(Radius<SourceMinRad || Radius>SourceMaxRad) {continue;}
                            SourceIntCut += SourceVal;
                            //
                            WaveFunVal=0;
                            int uTotMomBin = ???// TODO !!! q_total = sqrt(q_o^2 + q_s^2 + q_l^2)
                            for(unsigned usCh=0; usCh<NumCh; usCh++) WaveFunVal+=(WaveFunction2[uTotMomBin][uGrid][usCh])*ChannelWeight[usCh];
        
                            //if(use_sofia_model){
                            //printf("\033[1;36m  DEBUG (2):\033[0m %u SourceVal*WaveFunVal = %.4f * %.4f\n",uMomBin,SourceVal,WaveFunVal);
                            //usleep(5e3);
                            //}
                            Integrand = SourceVal*WaveFunVal;
                            kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] += Integrand;
                            kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] += pow(Integrand*kSourceGrid3D[uMomBin_out][uMomBin_side][uMomBin_long]->GetGridError(uGrid),2);
                        }//uGrid
        
                        if(use_sofia_model){
                            kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] *= (SourceInt+LostInSofia)/SourceInt;
                            kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = sqrt(kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long])*(SourceInt+LostInSofia)/SourceInt;
                            SourceIntCut *= (SourceInt+LostInSofia)/SourceInt;
                            SourceInt += LostInSofia;
        
                            if(compute_true_source){
                                for(unsigned uRad=0; uRad<TheTrueSource->GetNbins(1); uRad++){
                                    for(unsigned uCth=0; uCth<TheTrueSource->GetNbins(2); uCth++){
                                        unsigned tot_bin = TheTrueSource->GetTotBin(uMomBin,uRad,uCth);
                                        TheTrueSource->ScaleBin(tot_bin, (SourceInt+LostInSofia)/SourceInt);
                                    }
                                }
                            }
        
                        }
        
                        //ideally this should not happen, but in case it does (i.e. the source is normalized to value above 1)
                        //here this is corrected for
                        if(SourceInt>1+1e-4){
                            if(Notifications>=nWarning)
                                printf("\033[1;33mWARNING:\033[0m The source had to be renormalized (Integrated value=%.4e), please check the validity of your input function!\n",SourceInt);
                            //printf(" SourceInt=%f\n",SourceInt);
                            //if(AnaSourcePar)
                            //printf("AnaSourcePar: p0=%e; p1=%e\n",AnaSourcePar->GetParameter(0),AnaSourcePar->GetParameter(1));
                            //if(ForwardedSourcePar)
                            //printf("ForwardedSourcePar: p0=%e; p1=%e\n",ForwardedSourcePar->GetParameter(0),ForwardedSourcePar->GetParameter(1));
                            kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] /= SourceInt;
                            kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = sqrt(kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long])/SourceInt;
                            SourceIntCut/=SourceInt;
                            SourceInt=1.;
        
                            if(compute_true_source){
                                for(unsigned uRad=0; uRad<TheTrueSource->GetNbins(1); uRad++){
                                    for(unsigned uCth=0; uCth<TheTrueSource->GetNbins(2); uCth++){
                                        unsigned tot_bin = TheTrueSource->GetTotBin(uMomBin,uRad,uCth);
                                        TheTrueSource->ScaleBin(tot_bin, 1./SourceInt);
                                    }
                                }
                            }                    
                        }
        
                        //if the source is assumed normalized, any deviation below unity is considered to be an effect of a
                        //long-range tail that is unaccounted for due to the numerical precision. This tail should result in
                        //a flat residual that is added to the correlation
                        if(NormalizedSource && SourceIntCut<1.){
                            //this is an additional scale factor, evaluated based on the weight of channels with non-zero asymptotic solution
                            double AsympWeight = 0;
                            for(unsigned short usCh=0; usCh<NumCh; usCh++){
                              AsympWeight += ChannelWeight[usCh]*(!OnlyNumPw[usCh]);
                            }
                            kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] += (1.-double(SourceIntCut))*AsympWeight;
                            kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = sqrt(kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long]);
                            //printf("Added %f to C(%.0f)\n",1.-double(SourceIntCut),GetMomentum(uMomBin));
                        }
                        //else simply re-normalize the source to the total value of the integral (even if below 1)
                        //I am not sure if this even makes sense, consider in the future to use only the NormalizedSource case
                        else{
//printf        ("We are in a questionable state: %f at C(%.0f)\n",SourceIntCut,GetMomentum(uMomBin));
                            kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] *= double(SourceInt)/double(SourceIntCut);
                            kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = sqrt(kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long])*double(SourceInt)/double(SourceIntCut);
                        }
        
                        kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long] *= double(SourceInt)/double(SourceIntCut);
                        kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long] = sqrt(kCorrFunErr[uMomBin_out][uMomBin_side][uMomBin_long])*double(SourceInt)/double(SourceIntCut);
        
                        //printf("\033[1;36mDEBUG (2):\033[0m Integrated value=%.4e <> %.4e kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long]=%.4e\n",SourceInt, SourceIntCut, kCorrFun[uMomBin_out][uMomBin_side][uMomBin_long]);
                    }//if analytic source
                }
            }
        }
    }//uMom
    ComputedCorrFunction = true;
}