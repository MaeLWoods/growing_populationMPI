
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                                                    ///
/// Evolutionary model ///  Model 1: insertions + deletions Model 2: insertions + deletions + translocations                           ///
///                         Model 3: insertions + deletions + selection Model 4: insertions + deletions + translocations + selection   ///
/// ------------------------------------------------------------------------------------------------------------                       ///
///                                                                                                                                    ///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <random>
#include <sstream>
#include <ctime>
#include <mpi.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

using namespace std;

int read_params(const string& filename, int id);

int read_constnts( int id, const string& filename, int &NChr, int &Npop, int &ngen, double &lowB, double &chrlength, double &minchrlen, double &Curvemax, double &minchr, double &mu_b, double &minSV, double &p_birth, double &p_death, int &Nparticles, int &MaxCellNo , int &MaxMutNo, double &Lchr_loss, double &Uchr_loss, double &Lp_tran, double &Up_tran, double &Lmu_ki, double &Umu_ki, double &Lmu_kd, double &Umu_kd, double &Lfitness, double &Ufitness, double &Lmaxchr, double &Umaxchr, double &Lgnmdou, double &Ugnmdou );

double rtest(int id, gsl_rng* r, double a, double b);

double sgn(int id, double r);

double rtest(int beta, gsl_rng* r, double a, double b){
double a1=gsl_rng_uniform (r);
  double myrandom = a + (b-a)*a1;
 	while(myrandom==0){
double a2=gsl_rng_uniform (r);
 	myrandom = a + (b-a)*a2;
 	}
 	return myrandom;

}


double sgn(int id, double r1){
  
     if(id==0){
		double sgnr;

		if(r1<0){
		sgnr = -1;
		}
		else if(r1==0){
		sgnr=0;
		}
		else{
		sgnr=1;
		}

 	return sgnr;
 	}
}


//////////////////////////////////////////////////////////
///                                                    ///
///   MAIN                                             ///
///                                                    ///
//////////////////////////////////////////////////////////
int main (int argc, char *argv[]) {
int id;
int my_proc;

  gsl_rng * r;
  const gsl_rng_type * T;
 
    int NChr = 0;
    int Npop = 0;
    int ngen = 0; 
    int Nparticles = 0;
    int MaxCellNo = 0; 
    int MaxMutNo = 0;
    double lowB = 0;
    double chrlength = 0;
    double minchrlen = 0;
    double Curvemax = 0;
    double minchr = 0; 
    double mu_b = 0; 
    double minSV = 0;
    double p_birth = 0;
    double p_death = 0;
    double Lchr_loss = 0;
    double Uchr_loss = 0;
    double Lp_tran = 0;
    double Up_tran = 0;
    double Lmu_ki = 0;
    double Umu_ki = 0;
    double Lmu_kd = 0;
    double Umu_kd = 0;
    double Lfitness= 0;
    double Ufitness = 0;
    double Lmaxchr = 0;
    double Umaxchr = 0;
    double LSV_mean = 0;
    double USV_mean = 0;
    double Lgnmdou = 0;
    double Ugnmdou = 0;
    int set_seed = 0;
    double Final_ep = 1;
    int Npart = 1;
    int Job = 1;
    double alph = 1;
    int nparam = 1;
    int position = -1;
  
    int extinction;
    int Ncurr = 0;
	int Ndead = 0;
  
    string CFile = "ConstPars";
    unsigned int nConstPars = 0;
    
     // Initialise MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &my_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
 
    if(id==0){
  nConstPars = read_constnts(id,CFile, NChr, Npop, ngen, lowB, chrlength, minchrlen, Curvemax, minchr, mu_b, minSV,p_birth,p_death, Nparticles, MaxCellNo,MaxMutNo,Lchr_loss,Uchr_loss,Lp_tran,Up_tran,Lmu_ki,Umu_ki,Lmu_kd,Umu_kd,Lfitness, Ufitness, Lmaxchr, Umaxchr, Lgnmdou,Ugnmdou);
 }
 MPI_Barrier(MPI_COMM_WORLD);
 	MPI_Bcast(&ngen, 1, MPI_INT, 0, MPI_COMM_WORLD);
 	MPI_Bcast(&MaxCellNo, 1, MPI_INT, 0, MPI_COMM_WORLD);
 	MPI_Bcast(&MaxMutNo, 1, MPI_INT, 0, MPI_COMM_WORLD);
 	MPI_Barrier(MPI_COMM_WORLD);

 	double Mprev[MaxCellNo*NChr*NChr];
	double MCprev[MaxCellNo*NChr*NChr];
	double NTprev[MaxCellNo*NChr*NChr];
	double NIprev[MaxCellNo*NChr*NChr];
	double NDprev[MaxCellNo*NChr*NChr];
	double Cprev[MaxCellNo*NChr];
	double DivCprev[MaxCellNo*NChr];
	double CMixprev[MaxCellNo*NChr];
	double CSize[MaxCellNo*NChr];
	double M[MaxCellNo*NChr*NChr];
	double MC[MaxCellNo*NChr*NChr];
	double C[MaxCellNo*NChr];
	double DivC[MaxCellNo*NChr];
	double CMix[MaxCellNo*NChr];
	double NT[MaxCellNo*NChr*NChr];
	double NI[MaxCellNo*NChr*NChr];
	double ND[MaxCellNo*NChr*NChr];
	int n_chr[MaxCellNo];
	int Iborn[MaxCellNo];
	int Idead[MaxCellNo];
	int n_chrprev[MaxCellNo];
	double CmaxL[MaxCellNo];
	double CminL[MaxCellNo];
	int seed;
	double qq;
	///Rates
	double rdel[MaxCellNo*NChr];
	double rins[MaxCellNo*NChr];
	double rtrans[MaxCellNo*NChr];
	double rdelprev[MaxCellNo*NChr];
	double rinsprev[MaxCellNo*NChr];
	double rtransprev[MaxCellNo*NChr];
	
	double parameters[nparam];
	
	double Print_rtrans[ngen];
	double Print_rdel[ngen];
	double Print_rins[ngen];
	double Print_Hamming[ngen];
	int Print_Ncells[ngen];
	int Print_Ndead[ngen];
	int Print_Nborn[ngen];
	

    int Iremaining[MaxCellNo];
    int Nborn;
    int Ncurrprev;
    int Nnextlength;
    int Hsizeprev;
    int Hmsize;
    int Hmutval;
    
	
	int Nmut[MaxCellNo];
	int Hmut[MaxCellNo*MaxMutNo];
	int Hms[MaxCellNo];
	
	int Ngensample;	
    //parameters inferred
    double chr_loss = 0;
    double p_tran = 0;     
    double mu_ki = 0;
    double mu_kd = 0;
    double fitness = 0;
    double maxchr = 0; 
    double rand_forbirth = 0; 
    double rand_fordeath = 0; 

  MPI_Barrier(MPI_COMM_WORLD);
  if (id==0)
    {
      
	

  
  ////////////////////////////////////////////////////////////
  ///Generating parameters and random seed from bash script///
  ////////////////////////////////////////////////////////////
	
  for (int CommandLineCounter=1; CommandLineCounter < argc; CommandLineCounter++) {
  

				
    for (int StringPosition=0; StringPosition < strlen(argv[CommandLineCounter]); StringPosition++) {
						
      int EqualPosition = -1;
			
      if (argv[CommandLineCounter][StringPosition]=='=') {

	EqualPosition = StringPosition;				
				
      }
			
      if (EqualPosition > 1) {
	
	int VariableNameLength = EqualPosition;
	
	char VariableNameChar[VariableNameLength];
				
	for (int i=0; i<EqualPosition; i++) {
					
	  VariableNameChar[i] = argv[CommandLineCounter][i];

	}
				
	int VariableStringLength = strlen(argv[CommandLineCounter])-EqualPosition;

	char VariableChar[VariableStringLength];
				
	for (int i=0; i<VariableStringLength; i++) {
					
	  VariableChar[i] = argv[CommandLineCounter][i+1+EqualPosition];
					
	}
								
	double t1,t2; 
	char Bunk[2];
	Bunk[0] = VariableNameChar[0];
	Bunk[1] = VariableNameChar[1];
	
	if(Bunk[0]=='S'){
	  
	  t1 = strtod(VariableChar,NULL);
	  set_seed = int(t1);
	  					
	}else if(Bunk[0]=='P'){
	  
	  Npart = atoi(VariableChar);

	}else if(Bunk[0]=='F'){
	  
	  Final_ep = atof(VariableChar);
	 				
	}else if(Bunk[0]=='N'){
	  
	  nparam = atoi(VariableChar);
	 				
	}else if(Bunk[0]=='J'){

	  Job = atoi(VariableChar);

	}else if(Bunk[0]=='A'){

	  alph = atof(VariableChar);

	}else if(Bunk[0]=='T'){

	  position = atoi(VariableChar);

	}
      }
    }
  }

  cout << "Read input arguments" << endl;
  cout << "\tNparticles : " << Npart << endl;
  cout << "\tset_seed  : " << set_seed << endl;
  cout << "\tFinal_ep      : " << Final_ep << endl;
  cout << "\talpha       : " << alph << endl;
  cout << "\tposition       : " << position << endl;
	
   ////////////////////////////////////////////////////////////
   ///   Loop over the particles                            ///
   ////////////////////////////////////////////////////////////


seed = Job;
  extinction = 0;
  }
  
  MPI_Bcast(&Job, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Barrier(MPI_COMM_WORLD);
 

  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
    
 MPI_Barrier(MPI_COMM_WORLD);
  ////////////////////////////////////////////////////////////
  ///   Containers                                         ///
  ////////////////////////////////////////////////////////////
  gsl_rng_set(r, set_seed );
    
  	if(id==0){


	cout << "Begin initialization: " << endl;
	
	std::default_random_engine generatorinitial;
    normal_distribution<double> distributioninitial(0.0,1);
	cout << "after normal distribution: " << endl;
    for (int i = 0; i < MaxCellNo; i++) {
    if(i<Npop){
    
	  Hmut[i*MaxMutNo + 0] =  i;

	  Hms[i] = 1;
	  
	  for(int s=1; s<MaxMutNo; s++){
	  
	  Hmut[i*MaxMutNo + s] = 0;
	  
	  
	  }

     
	  Nmut[i] = 0;
	  n_chr[i] = NChr;
      
      //step containers
      n_chrprev[i] = NChr;
      double Max_ChrSize[NChr];
	
      for(int j = 0; j < NChr; j++){
			  double rlen = ( rtest(id,r,0,1) );
			  C[i*NChr+j] = minchrlen + (chrlength-minchrlen)*rlen;
			  CMix[i*NChr+j] = minchrlen + (chrlength-minchrlen)*rlen;
			  CSize[i*NChr+j] = minchrlen + (chrlength-minchrlen)*rlen;
	         
			  Cprev[i*NChr+j]=C[i*NChr+j];
			  CMixprev[i*NChr+j]=CMix[i*NChr+j];
			  Max_ChrSize[j] = C[i*NChr+j];		
	
	DivC[i*NChr+j]=0;
	DivCprev[i*NChr+j]=DivC[i*NChr+j];
    
      }
      CmaxL[i] = *max_element(Max_ChrSize,Max_ChrSize+NChr);
      CminL[i] = *min_element(Max_ChrSize,Max_ChrSize+NChr);
      for(int j = 0; j < NChr; j++){
      for(int k = 0; k < NChr; k++){
	
		NT[i*NChr*NChr+j*NChr+k]=0;
		NI[i*NChr*NChr+j*NChr+k]=0;
		ND[i*NChr*NChr+j*NChr+k]=0;
		NTprev[i*NChr*NChr+j*NChr+k]=NT[i*NChr*NChr+j*NChr+k];
	    NIprev[i*NChr*NChr+j*NChr+k]=NI[i*NChr*NChr+j*NChr+k];
	    NDprev[i*NChr*NChr+j*NChr+k]=ND[i*NChr*NChr+j*NChr+k];
	    
	  if(j==k){
	    M[i*NChr*NChr+j*NChr+k]=1;
	    Mprev[i*NChr*NChr+j*NChr+k]=M[i*NChr*NChr+j*NChr+k];
	    MC[i*NChr*NChr+j*NChr+k]=C[i*NChr+k];
	    MCprev[i*NChr*NChr+j*NChr+k]=MC[i*NChr*NChr+j*NChr+k];
	  }
	  else{
	      M[i*NChr*NChr+j*NChr+k]=0;
	      MC[i*NChr*NChr+j*NChr+k]=0;
		  Mprev[i*NChr*NChr+j*NChr+k]=0;
		  MCprev[i*NChr*NChr+j*NChr+k]=0;
	  		}
		 }
		 
		 }
    
    }
    else{
    
      Hmut[i*MaxMutNo + 0] =  0;
	  Hmut[i*MaxMutNo + 1] =  0;
	  Hmut[i*MaxMutNo + 2] =  0;
	  
	  Hms[i] = 1;
	  for(int s=1; s<MaxMutNo; s++){
	  
	  Hmut[i*MaxMutNo + s] = 0;
	  
	  
	  }
    
	  Nmut[i] = 0;
	  n_chr[i] = 0;
      
      //step containers
      n_chrprev[i] = 0;
      double Max_ChrSize[NChr];
	
      for(int j = 0; j < NChr; j++){
	
			  double rlen = ( rtest(id,r,0,1) );
			  C[i*NChr+j] = 0;
			  CMix[i*NChr+j] = 0;
			  CSize[i*NChr+j] = 0;
	         
			  Cprev[i*NChr+j]= 0;
			  CMixprev[i*NChr+j]= 0;
			  Max_ChrSize[j] = 0;		
	
	DivC[i*NChr+j]=0;
	DivCprev[i*NChr+j]= 0;
    
	
      }
      
      CmaxL[i] = 0;
      CminL[i] = 0;
            
      for(int j = 0; j < NChr; j++){
      for(int k = 0; k < NChr; k++){
	
		NT[i*NChr*NChr+j*NChr+k]=0;
		NI[i*NChr*NChr+j*NChr+k]=0;
		ND[i*NChr*NChr+j*NChr+k]=0;
		NTprev[i*NChr*NChr+j*NChr+k]= 0;
	    NIprev[i*NChr*NChr+j*NChr+k]= 0;
	    NDprev[i*NChr*NChr+j*NChr+k]= 0;
	    
	  if(j==k){
	    M[i*NChr*NChr+j*NChr+k]= 0;
	    Mprev[i*NChr*NChr+j*NChr+k]= 0;
	    MC[i*NChr*NChr+j*NChr+k]= 0;
	    MCprev[i*NChr*NChr+j*NChr+k]= 0;
	  }
	  else{
	      M[i*NChr*NChr+j*NChr+k]=0;
	      MC[i*NChr*NChr+j*NChr+k]=0;
		  Mprev[i*NChr*NChr+j*NChr+k]=0;
		  MCprev[i*NChr*NChr+j*NChr+k]=0;
	  		}
		 }
		 
		 }
    
    }
    
    }
    

	/////////////////////////
	/// Read in the files ///
    /////////////////////////
    
	 char DatF[100];
      int DnF;
      DnF=sprintf(DatF,"param_folder/param-M%d.dat",position);
   

        ifstream infile (DatF);

  int counter = 0; 
  if (infile.is_open()){

    string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

		vector<string> split;
      string buf; 
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

	
	
      for(int k=0; k<nparam; k++){
     parameters[k] = atof( split[k].c_str() ) ;
     
     }
     
    }
      
  }
  else{
    cerr << "Error: open of constant parameters file unsuccessful: " << endl;
  }
    

	
	////////////////////////////////////////////////////////////
    /// Initialize mutation rates per base per cell division ///
    ////////////////////////////////////////////////////////////
    
    
    double trans_rate = 0;
    double var = 0;
    
      trans_rate = parameters[1];
    cout << "trans_rate: " << trans_rate << endl;
      var = parameters[2];
	cout << "var: " << var << endl;
      
      for(int i=0; i<MaxCellNo; i++){
      if(i<Npop){
      
      /////////////////////////////////////////////////
      ///Set cell independent rates of translocation///
      /////////////////////////////////////////////////
      
		double qq_t = trans_rate;
		double p_t = mu_b*pow(10,qq_t);		
		
	  /////////////////////////////////////////////////
      ///Set chromosome independent rates of indels ///
      /////////////////////////////////////////////////
      for(int j = 0; j < NChr; j++){

		int realrates = 0;
		double p_ins = 0;
		double p_del = 0;
		while(realrates==0){

			double rins1 = 0.2;//HERE( runiform(id,r,stream,-0.499,0.499 ) );
			double rins2 = 0.18;//HERE( runiform(id,r,stream,-0.499,0.499 ) );
			double sgnr1 = sgn(id,rins1);
			double sgnr2 = sgn(id,rins2);
			double qqi = -1*0.75*sgnr1*log(1-2*abs(rins1));
			double qqd = -1*0.75*sgnr2*log(1-2*abs(rins2));
		
			p_ins = 4-qqi;
			p_del = 4-qqd;
		
			if((p_ins>0)&&(p_ins<5.5)&&(p_del>0)&&(p_del<5.5)){
			 realrates = 1;
			}
		}
		
        rins[i*NChr+j] = mu_b*pow(10,p_ins);
		rdel[i*NChr+j] = mu_b*pow(10,p_del);	
		rdelprev[i*NChr+j] = rdel[i*NChr+j];
		rinsprev[i*NChr+j] = rins[i*NChr+j];
		rtrans[i*NChr+j]= p_t;
		rtransprev[i*NChr+j] = rtrans[i*NChr+j];
	  
      }
      
      }
      else{
      for(int j=0; j<NChr; j++){
       rins[i*NChr+j] = 0;
		rdel[i*NChr+j] = 0;	
		rdelprev[i*NChr+j] = 0;
		rinsprev[i*NChr+j] = 0;
		rtrans[i*NChr+j]= 0;
		rtransprev[i*NChr+j] = 0;
      }
      }
      
    }
    
    int TimeUntilDeath = 0;    
     	
    ////////////////////////////////////////////////////////////
    ///   Loop over the generations                          ///
    ////////////////////////////////////////////////////////////
    //cout << "Looping over generations" << endl;

    
    cout << "start of generations: " << endl;
	Ncurr = Npop;
	Ndead = 0;
	Hsizeprev =3;
	Hmsize =3;
	
	
	}


   
    //gsl_rng_set(r, set_seed );
 
   //qq=rtest(id, r, 0, 1);


	
	 MPI_Barrier(MPI_COMM_WORLD);
	
      for(int gen_clock=1; gen_clock<ngen; gen_clock++){

      
      //////////////////////////////////////////////////////
      /////////Clear Hamming distance variables/////////////
      //////////////////////////////////////////////////////
     MPI_Barrier(MPI_COMM_WORLD); 
     if (id==0){
      cout << "start of generations: " << gen_clock << endl;
      Nborn = 0;
      Ncurrprev = Ncurr;
      Nnextlength = Ncurr;
      
      int ThisNDead = 0;
      
      for(int q=0; q<MaxCellNo;q++){
       Iborn[q] = 0;
       Idead[q] = 0;
      }
      
      for(int q=0; q<MaxCellNo;q++){
       Iremaining[q] = 0;
       if(q<Ncurr){
       Iremaining[q] = q;
       }
      }
      
      
      //////////////////////////////////////////////////////
      /////////Clear Hamming distance variables/////////////
      //////////////////////////////////////////////////////
      
      
      //cout << "gen_clock: " << gen_clock << endl;
      Ngensample = Ncurr;
      //cout << "Ncurr: " << Ncurr << endl;
      //cout << "Hsizeprev: " << Hsizeprev << endl;
      //////////////////////////////////////////////////////
      /////////Mutate: insertions and deletions/////////////
      //////////////////////////////////////////////////////

		
				int checkindelnum = 0;
				for(int i=0; i<Ngensample; i++){
				for(int g_d=0; g_d<n_chrprev[i]; g_d++){
				for(int k=0; k<NChr; k++){
		
				if(NIprev[i*NChr*NChr+ g_d*n_chrprev[i]+k]>30000){
				checkindelnum = 1;
				gen_clock=ngen;
				
			}}}}
		
		

    	for(int ng_count=0; ng_count<Ngensample; ng_count++){	
    	
    	//cout << "Ngensample: " << Ngensample << endl;
    	
    	if((Ncurr<3)||(Ncurr>(MaxCellNo-100))){
      extinction = 1;
      //ng_count = Ngensample;
      //gen_clock = ngen;
      
      }
    	
    if(extinction==0){
    	//cout << " before rand birth gen_clock: " << gen_clock << endl;

        
    	rand_forbirth = (   rtest(id,r,0,1) ); 
    	rand_fordeath = (   rtest(id,r,0,1) ); 
    	
    	if(rand_forbirth<p_birth){
    	 		
    	 		//cout << " goes into birth: " << gen_clock << endl;		
    	
    	
								double pi_q = (   rtest(id,r,0,1) ); 
								int i = floor(pi_q*Ngensample);
								
	
								while(i>(Ncurr-1)){
								double pi_qq = (   rtest(id,r,0,1) ); 
								i = floor(pi_qq*Ngensample);
	
								}
								
								
								default_random_engine generatorindels;
								normal_distribution<double> distributionindels(0.0,1.5);
								//////////////////////////////////////////////////////
								/////////Mutate: insertions and deletions/////////////
								//////////////////////////////////////////////////////
								int anymutation = 0;
								int number_ins = 0;
								int number_dels = 0;
								int size_of_Hmut = 0;
								
								//cout << " goes into birth before indels: " << gen_clock << endl;	
		
								for(int g_d=0; g_d<n_chrprev[i]; g_d++){
								for(int k=0; k<NChr; k++){
														
										 int N_i = gsl_ran_poisson (r, MCprev[i*NChr*NChr+g_d*NChr+k]*rinsprev[i*NChr+k]);
										 int N_d = gsl_ran_poisson (r, MCprev[i*NChr*NChr+g_d*NChr+k]*rdelprev[i*NChr+k]);	
				
											  for(int tot=0; tot<N_i; tot++){
											
												double qq = distributionindels(generatorindels);		
												double sz = exp(15 + qq);
												if(sz<minSV){}
												else{
												DivCprev[i*NChr+k] += sz;
												Cprev[i*NChr+k] = Cprev[i*NChr+k] + sz;
												MCprev[i*NChr*NChr+g_d*NChr+k] = MCprev[i*NChr*NChr+g_d*NChr+k] + sz;
												NIprev[i*NChr*NChr+g_d*NChr+k] = NIprev[i*NChr*NChr+g_d*NChr+k] + 1;
						
												anymutation = 1;
												size_of_Hmut += 1;
												number_ins += 1;
												}
						
												double Cmax = 1.1*CmaxL[i];
												double Cmin = 0.9*CminL[i];
												double checklen = 0;
												for(int p=0; p<NChr; p++){
						
													checklen += MCprev[i*NChr*NChr+g_d*NChr+p];
						
												}
						
												if(checklen>10*Cmax){
						
												tot = N_i;
						
												}
		  
											  }
											  for(int tot=0; tot<N_d; tot++){					  	
					  
												double qq = distributionindels(generatorindels);		
												double sz = exp(15 + qq);
						
												if((MCprev[i*NChr*NChr+g_d*NChr+k]>=0)&&(sz>MCprev[i*NChr*NChr+g_d*NChr+k])){
													sz = MCprev[i*NChr*NChr+g_d*NChr+k];
													tot = N_d;
														}
						
												if(sz<minSV){}
												else{
												DivCprev[i*NChr+k] -= sz;
												Cprev[i*NChr+k] = Cprev[i*NChr+k] - sz;
												MCprev[i*NChr*NChr+g_d*NChr+k] = MCprev[i*NChr*NChr+g_d*NChr+k] - sz;
												NDprev[i*NChr*NChr+g_d*NChr+k] = NDprev[i*NChr*NChr+g_d*NChr+k] + 1;
						    					number_dels += 1;
												anymutation = 1;
												size_of_Hmut += 1;
												}
		
												}	     
											}	
						
										}
										
										//cout << " after indels: " << gen_clock << endl;
								//////////////////////////////////////////////////////
								/////////End indels loop//////////////////////////////
								//////////////////////////////////////////////////////
					
    	
								/////////////////////////////
								/////////UpdateM/////////////
								/////////////////////////////
								for(int g_d=0; g_d<n_chrprev[i]; g_d++){
		
			
									double sum_j = 0;
			
									for(int k=0; k<NChr; k++){
			
									
		
									double cprev_size = 0;
		
									cprev_size = Cprev[i*NChr + k];

		
									  //Avoid seg fault
									  if(cprev_size==0){
									  Mprev[i*NChr*NChr+g_d*NChr+k] = 0;
									  }
									  else{
									  Mprev[i*NChr*NChr+g_d*NChr+k] = MCprev[i*NChr*NChr+g_d*NChr+k]/cprev_size;
									  }
				
									  sum_j = Mprev[i*NChr*NChr+g_d*NChr+k]*Cprev[i*NChr+k] + sum_j;
							
									}
			
									CMixprev[i*NChr+g_d]=sum_j;
	  
								}
		
		                        // cout << "After update M: " << endl;
								//////////////////////////////////////////////////////
								//    Translocation
								//---------------------------------		
								

								int anytrans = 0;
								int number_trans = 0;
								
								
		 
							for(int j=0; j<n_chrprev[i]; j++){
							
							

 
								  for(int k=0; k<NChr; k++){
								  
								 
								  //default_random_engine generatortrans;
  									//	poisson_distribution<int> distributiontrans(MCprev[i*NChr*NChr+j*NChr+k]*rtransprev[i*NChr+k]);
								  //translocate between C1 and C2//
								 // int N_t = distributiontrans(generatortrans);
									int N_t = gsl_ran_poisson (r, MCprev[i*NChr*NChr+j*NChr+k]*rtransprev[i*NChr+k]);
									double f_t = rtest(id,r,0,1);
									
									 

									if(N_t>0){
									  ///Mutate: second chomsome///
									  for(int j1=j; j1<n_chrprev[i]; j1++){
		  
					
												if(j!=j1){
			
												  for(int k1=0; k1<NChr; k1++){
						  
															if(k1!=k){
								  
								  						//		default_random_engine generatortrans1;
  										//poisson_distribution<int> distributiontrans1(MCprev[i*NChr*NChr+j1*NChr+k1]*rtransprev[i*NChr+k1]);
  										//int N_t1 = distributiontrans1(generatortrans1);
															  int N_t1 = gsl_ran_poisson (r, MCprev[i*NChr*NChr+j1*NChr+k1]*rtransprev[i*NChr+k1]);
															  double f_tt = rtest(id,r,0,1);
															  if(N_t1>0){		
																		double sum_j0 = 0;
																		double sum_j1 = 0;
			
																		for(int ka=0; ka<NChr; ka++){
												

																				  double a1 = Mprev[i*NChr*NChr+j*NChr+ka];
																				  double b1 = Mprev[i*NChr*NChr+j1*NChr+ka];
				  
																				  double p = lowB + (1-2*lowB)*rtest(id,r,0,1);
				  
																				  Mprev[i*NChr*NChr+j*NChr+ka] = p*(a1+b1);
																				  Mprev[i*NChr*NChr+j1*NChr+ka] = (1-p)*(a1+b1);
														  
																				  sum_j0 += Mprev[i*NChr*NChr+j*NChr+ka]*Cprev[i*NChr+ka];
																				  sum_j1 += Mprev[i*NChr*NChr+j1*NChr+ka]*Cprev[i*NChr+ka];
																				}
			
																		CMixprev[i*NChr+j]=sum_j0;
																		CMixprev[i*NChr+j1]=sum_j1;
			
																		NTprev[i*NChr*NChr+j*NChr+k] = NTprev[i*NChr*NChr+j*NChr+k] + 1;
																		NTprev[i*NChr*NChr+j1*NChr+k1] = NTprev[i*NChr*NChr+j1*NChr+k1] + 1;
																		anytrans = 1;
																		number_trans += 1;
			
																	}
											
															}
						
														}
					
													}
			
										//end if not the same dynamic chr
			
										 }
				 
										 //loop over dynamic chr
									}
									//end mutate second chr
								  }
			
								}
 
 								
	
								if((anytrans == 1)||(anymutation == 1)){
													Hmsize = Hsizeprev + 1;
													Hmutval = 1;							
												}		
									
	  
								//////////////////////////////////////////////////////
								/////////Update M & MC ///////////////////////////////
								//////////////////////////////////////////////////////
		
								for(int g_d=0; g_d<n_chrprev[i]; g_d++){
		
			
									double sum_j = 0;
			
									for(int k=0; k<NChr; k++){
				
									  MCprev[i*NChr*NChr+g_d*NChr+k] = Mprev[i*NChr*NChr+g_d*NChr+k]*Cprev[i*NChr+k];
				
									  sum_j = Mprev[i*NChr*NChr+g_d*NChr+k]*Cprev[i*NChr+k] + sum_j;
							
									}
			
									CMixprev[i*NChr+g_d]=sum_j;
		
	  
								}
								
								//cout << " Checkpoint 3 " << endl;
								
								//////////////////////////////////////////////////////
								/////////Selection                        ////////////
								//////////////////////////////////////////////////////
								int select_out = 0;
		
								for(int k=0; k<NChr; k++){
		
		
								if(CMixprev[i*NChr+k]!=Cprev[i*NChr+k]){
		
								}

								double Cmax = 1.1*CmaxL[i];
								double Cmin = 0.9*CminL[i];
		
								if((CMixprev[i*NChr + k]>(Cmax))||(CMixprev[i*NChr + k]<(Cmin))){
								select_out = 1;
								}
		
								}
								
								
		
								//cout << "After selection: " << endl;
								//////////////////////////////////////////////////////
								/////////Resample                        /////////////
								//////////////////////////////////////////////////////
				
								if(select_out==1){
																
				                 Ncurr -= 1;
				                 Ndead += 1;
				                 ThisNDead += 1;
				                 
				                 if(i<Nnextlength){
				                 Nnextlength -= 1;
				                 
				                 
				                 for(int q=0; q<Nnextlength; q++){
				                 if(q==i){
				                 Idead[ThisNDead-1] = Iremaining[q];
				                 Iremaining[q] = 0;
				                 for(int r=q; r<(Nnextlength); r++){
				                 Iremaining[r] = Iremaining[r+1];
				                 }
				                 }
				                 }
				                 
				                 Iremaining[Nnextlength] = 0;
				                 
				                 }
				                 
				                 int s = i+1; 
				                 int Ncurrprev = Ncurr+2;
				                 for(int s_tilde=i; s_tilde<Ncurrprev; s_tilde++){
				                 
				                 	if(s_tilde==(Ncurr+1)){
				                 		
				                 		n_chrprev[s_tilde] = 0;
				                     CmaxL[s_tilde] = 0;
				                     CminL[s_tilde] = 0;
				                     n_chr[s_tilde] = 0;
				                     Nmut[s_tilde] = 0;
				                     
				                     for(int p=0; p<Hms[s_tilde]; p++){
				                     
				                     	Hmut[s_tilde*MaxMutNo + p] = 0;
				                     
				                     }
				                     
				                     Hms[s_tilde]=1;
				                     
				                 	for(int a=0; a<NChr; a++){
				                 		
				                 			CSize[s_tilde*NChr + a] = 0;
				                 			Cprev[s_tilde*NChr + a] = 0;
				                 			DivCprev[s_tilde*NChr + a] = 0;
				                 			CMixprev[s_tilde*NChr + a] = 0;
				                 			
				                 			rdelprev[s_tilde*NChr + a] = 0;
				                 			rinsprev[s_tilde*NChr + a] = 0;
				                 			rtransprev[s_tilde*NChr + a] = 0;
				                 			
				                 			rdel[s_tilde*NChr + a] = 0;
				                 			rins[s_tilde*NChr + a] = 0;
				                 			rtrans[s_tilde*NChr + a] = 0;
				                 			
				                 			C[s_tilde*NChr + a] = 0;
				                 			DivC[s_tilde*NChr + a] = 0;
				                 			CMix[s_tilde*NChr + a] = 0;
				                 						                 		
				                 		for(int d=0; d<NChr; d++){
				                 		   
				                 		    NTprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    NIprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    NDprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    
				                 		    NT[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    NI[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    ND[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    
				                 		    Mprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    MCprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    
				                 		    M[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    MC[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		
				                 		}
				                 		
				                 			
				                 		
				                 		}
				                 		
				                 		
				                 	
				                 	}
				                 	else{
				                 				                 
				                     n_chrprev[s_tilde] = n_chrprev[s];
				                     CmaxL[s_tilde] = CmaxL[s];
				                     CminL[s_tilde] = CminL[s];
				                     n_chr[s_tilde] = n_chr[s];
				                     Nmut[s_tilde] = Nmut[s];
				                     
				                  
				                     
				                     for(int p=0; p<Hms[s_tilde]; p++){
				                     Hmut[s_tilde*MaxMutNo + p] = 0;
				                     }
				                     for(int p=0; p<Hms[s]; p++){
				                     	Hmut[s_tilde*MaxMutNo + p] = Hmut[s*MaxMutNo + p];
				                     
				                     }
				                     
				                     Hms[s_tilde] = Hms[s];
				                     
				                 	for(int a=0; a<NChr; a++){
				                 		
				                 			CSize[s_tilde*NChr + a] = CSize[s*NChr + a];
				                 			Cprev[s_tilde*NChr + a] = Cprev[s*NChr + a];
				                 			DivCprev[s_tilde*NChr + a] = DivCprev[s*NChr + a];
				                 			CMixprev[s_tilde*NChr + a] = CMixprev[s*NChr + a];
				                 			
				                 			rdelprev[s_tilde*NChr + a] = rdelprev[s*NChr + a];
				                 			rinsprev[s_tilde*NChr + a] = rinsprev[s*NChr + a];
				                 			rtransprev[s_tilde*NChr + a] = rtransprev[s*NChr + a];
				                 			
				                 			rdel[s_tilde*NChr + a] = rdel[s*NChr + a];
				                 			rins[s_tilde*NChr + a] = rins[s*NChr + a];
				                 			
				                 			rtrans[s_tilde*NChr + a] = rtrans[s*NChr + a];
				                 			
				                 			C[s_tilde*NChr + a] = C[s*NChr + a];
				                 			DivC[s_tilde*NChr + a] = DivC[s*NChr + a];
				                 			CMix[s_tilde*NChr + a] = CMix[s*NChr + a];
				                 						                 		
				                 		for(int d=0; d<NChr; d++){
				                 		   
				                 		    NTprev[s_tilde*NChr*NChr + a*NChr + d] = NTprev[s*NChr*NChr + a*NChr + d];
				                 		    NIprev[s_tilde*NChr*NChr + a*NChr + d] = NIprev[s*NChr*NChr + a*NChr + d];
				                 		    NDprev[s_tilde*NChr*NChr + a*NChr + d] = NDprev[s*NChr*NChr + a*NChr + d];
				                 		    
				                 		    NT[s_tilde*NChr*NChr + a*NChr + d] = NT[s*NChr*NChr + a*NChr + d];
				                 		    NI[s_tilde*NChr*NChr + a*NChr + d] = NI[s*NChr*NChr + a*NChr + d];
				                 		    ND[s_tilde*NChr*NChr + a*NChr + d] = ND[s*NChr*NChr + a*NChr + d];
				                 		    
				                 		    Mprev[s_tilde*NChr*NChr + a*NChr + d] = Mprev[s*NChr*NChr + a*NChr + d];
				                 		    MCprev[s_tilde*NChr*NChr + a*NChr + d] = MCprev[s*NChr*NChr + a*NChr + d];
				                 		    
				                 		    M[s_tilde*NChr*NChr + a*NChr + d] = M[s*NChr*NChr + a*NChr + d];
				                 		    MC[s_tilde*NChr*NChr + a*NChr + d] = MC[s*NChr*NChr + a*NChr + d];
				                 		
				                 		}
				                 		
				                 			
				                 		
				                 		}
				                       s += 1;
				                      }
				                 	
				                 }
				                 
				                 Hmsize = Hsizeprev;
									
						  //End death rate
			
								}
								else{
							
									
								//replace one		
								    Ncurr += 1;
									int position = i;
									Nborn += 1;
									n_chrprev[(Ncurr-1)] = n_chrprev[i];
				                     CmaxL[(Ncurr-1)] = CmaxL[i];
				                     CminL[(Ncurr-1)] = CminL[i];
				                     n_chr[(Ncurr-1)] = n_chr[i];
				                     Nmut[(Ncurr-1)] = Nmut[i];
				                     
				                    
				                     if(i<Nnextlength){				                 
				                 
				                 Iborn[(Nborn-1)]=Iremaining[i];
				                 
				                 }
				                 else{
				                 Iborn[(Nborn-1)]=i-Nnextlength;
				                     }
				                     
				                    
				                     
				                     for(int s=0; s<Ncurr; s++){
				                     
				                     if(s==(Ncurr-1)){
				                     
				                     	
												if(s==(Ncurr-1)){
												    //cout << "Hmut size: " << Hms[i] << endl;
												    for(int h=0; h<Hms[i]; h++){
												    Hmut[s*MaxMutNo + h] = Hmut[i*MaxMutNo + h];
												    }
												   // cout << "Hms[i]: " << Hms[i] << endl;
													Hmut[s*MaxMutNo + Hms[i]] = Hmutval;
													Hmut[s*MaxMutNo + Hms[i]+1] = Hmsize-1;
													
													Hms[s] = Hms[i] + 2;
													
												}
												
											}
											
				                     	}
				                     
				                     
				                     Hsizeprev = Hmsize;
				                     
				                     
				           
				                     
				                 	for(int a=0; a<NChr; a++){
				                 	
				                 		
				                 			CSize[(Ncurr-1)*NChr + a] = CSize[i*NChr + a];
				                 			Cprev[(Ncurr-1)*NChr + a] = Cprev[i*NChr + a];
				                 			DivCprev[(Ncurr-1)*NChr + a] = DivCprev[i*NChr + a];
				                 			CMixprev[(Ncurr-1)*NChr + a] = CMixprev[i*NChr + a];
				                 			
				                 			rdelprev[(Ncurr-1)*NChr + a] = rdelprev[i*NChr + a];
				                 			rinsprev[(Ncurr-1)*NChr + a] = rinsprev[i*NChr + a];
				                 			rtransprev[(Ncurr-1)*NChr + a] = rtransprev[i*NChr + a];
				                 			
				                 			rdel[(Ncurr-1)*NChr + a] = rdel[i*NChr + a];
				                 			rins[(Ncurr-1)*NChr + a] = rins[i*NChr + a];
				                 			rtrans[(Ncurr-1)*NChr + a] = rtrans[i*NChr + a];
				                 			
				                 			C[(Ncurr-1)*NChr + a] = C[i*NChr + a];
				                 			DivC[(Ncurr-1)*NChr + a] = DivC[i*NChr + a];
				                 			CMix[(Ncurr-1)*NChr + a] = CMix[i*NChr + a];
				                 						                 		
				                 		for(int d=0; d<NChr; d++){
				                 		   
				                 		    NTprev[(Ncurr-1)*NChr*NChr + a*NChr + d] = NTprev[i*NChr*NChr + a*NChr + d];
				                 		    NIprev[(Ncurr-1)*NChr*NChr + a*NChr + d] = NIprev[i*NChr*NChr + a*NChr + d];
				                 		    NDprev[(Ncurr-1)*NChr*NChr + a*NChr + d] = NDprev[i*NChr*NChr + a*NChr + d];
				                 		    
				                 		    NT[(Ncurr-1)*NChr*NChr + a*NChr + d] = NT[i*NChr*NChr + a*NChr + d];
				                 		    NI[(Ncurr-1)*NChr*NChr + a*NChr + d] = NI[i*NChr*NChr + a*NChr + d];
				                 		    ND[(Ncurr-1)*NChr*NChr + a*NChr + d] = ND[i*NChr*NChr + a*NChr + d];
				                 		    
				                 		    Mprev[(Ncurr-1)*NChr*NChr + a*NChr + d] = Mprev[i*NChr*NChr + a*NChr + d];
				                 		    MCprev[(Ncurr-1)*NChr*NChr + a*NChr + d] = MCprev[i*NChr*NChr + a*NChr + d];
				                 		    
				                 		    M[(Ncurr-1)*NChr*NChr + a*NChr + d] = M[i*NChr*NChr + a*NChr + d];
				                 		    MC[(Ncurr-1)*NChr*NChr + a*NChr + d] = MC[i*NChr*NChr + a*NChr + d];
				                 		
				                 		}
				                 		
				                 		
				                 		}
				                 		
				                 
									  
									  for(int j = 0; j < NChr; j++){	
	
											Cprev[i*NChr + j]=C[i*NChr + j];
			
											for(int g_d=0; g_d<n_chrprev[i]; g_d++){
											Mprev[i*NChr*NChr + NChr*g_d + j] = M[i*NChr*NChr + NChr*g_d + j];
											MCprev[i*NChr*NChr + NChr*g_d + j] = MC[i*NChr*NChr + NChr*g_d + j];
											NTprev[i*NChr*NChr + NChr*g_d + j] = NT[i*NChr*NChr + NChr*g_d + j];
											NIprev[i*NChr*NChr + NChr*g_d + j] = NI[i*NChr*NChr + NChr*g_d + j];
											NDprev[i*NChr*NChr + NChr*g_d + j] = ND[i*NChr*NChr + NChr*g_d + j];
			
											}
											DivCprev[i*NChr + j] = DivC[i*NChr + j];
											CMixprev[i*NChr + j] = CMix[i*NChr + j];
											rdelprev[i*NChr + j] = rdel[i*NChr + j];
											rinsprev[i*NChr + j] = rins[i*NChr + j];
											rtransprev[i*NChr + j] = rtrans[i*NChr + j];
	
										  }
								
										n_chrprev[i] = n_chr[i];
										
									
				                
									
								}
							 
    	}
    	
    
    	
    	if(rand_fordeath<p_death){
    	
    	//cout << " goes into death: " << gen_clock << endl;	
    	
    	
    	
			double pi_q = (   rtest(id,r,0,1) ); 
			int i = floor(pi_q*Ngensample);
    
    	while(i>(Ncurr-1)){
    	
			double pi_qq = (   rtest(id,r,0,1) ); 
			i = floor(pi_qq*Ngensample);
    
    	}
    	
    					
    					
				                 Ncurr -= 1;
				                 Ndead += 1;
				                 ThisNDead += 1;
				                 
				                 if(i<Nnextlength){
				                 Nnextlength -= 1;
				                 
				                 
				                 for(int q=0; q<Nnextlength; q++){
				                 if(q==i){
				                 Idead[ThisNDead-1] = Iremaining[q];
				                 Iremaining[q] = 0;
				                 for(int r=q; r<(Nnextlength); r++){
				                 Iremaining[r] = Iremaining[r+1];
				                 }
				                 }
				                 }
				                 
				                 Iremaining[Nnextlength] = 0;
				                 
				                 }
				                 
				                
				                 
				                 int s = i+1; 
				                 int Ncurrprev = Ncurr+2;
				            
				                 for(int s_tilde=i; s_tilde<Ncurrprev; s_tilde++){
				                 
				                
				                
				                 	if(s_tilde==(Ncurr+1)){
				                 		
				                 		n_chrprev[s_tilde] = 0;
				                     CmaxL[s_tilde] = 0;
				                     CminL[s_tilde] = 0;
				                     n_chr[s_tilde] = 0;
				                     Nmut[s_tilde] = 0;
				                     
				                     for(int p=0; p<Hms[s_tilde]; p++){
				                     
				                     	Hmut[s_tilde*MaxMutNo + p] = 0;
				                     
				                     }
				                     //this is wrong
				                     Hms[s_tilde]=1;
				                     
				                 	for(int a=0; a<NChr; a++){
				                 		
				                 			CSize[s_tilde*NChr + a] = 0;
				                 			Cprev[s_tilde*NChr + a] = 0;
				                 			DivCprev[s_tilde*NChr + a] = 0;
				                 			CMixprev[s_tilde*NChr + a] = 0;
				                 			
				                 			rdelprev[s_tilde*NChr + a] = 0;
				                 			rinsprev[s_tilde*NChr + a] = 0;
				                 			rtransprev[s_tilde*NChr + a] = 0;
				                 			
				                 			rdel[s_tilde*NChr + a] = 0;
				                 			rins[s_tilde*NChr + a] = 0;
				                 			rtrans[s_tilde*NChr + a] = 0;
				                 			
				                 			C[s_tilde*NChr + a] = 0;
				                 			DivC[s_tilde*NChr + a] = 0;
				                 			CMix[s_tilde*NChr + a] = 0;
				                 						                 		
				                 		for(int d=0; d<NChr; d++){
				                 		   
				                 		    NTprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    NIprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    NDprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    
				                 		    NT[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    NI[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    ND[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    
				                 		    Mprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    MCprev[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    
				                 		    M[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		    MC[s_tilde*NChr*NChr + a*NChr + d] = 0;
				                 		
				                 		}
				                 		
				                 			
				                 		
				                 		}
				                 	
				                 	}
				                 	else{
				                 				                 
				                     n_chrprev[s_tilde] = n_chrprev[s];
				                     
				                     
				                     CmaxL[s_tilde] = CmaxL[s];
				                     CminL[s_tilde] = CminL[s];
				                     n_chr[s_tilde] = n_chr[s];
				                     Nmut[s_tilde] = Nmut[s];
				                     
				                      for(int p=0; p<Hms[s_tilde]; p++){
				                     Hmut[s_tilde*MaxMutNo + p] = 0;
				                     }
				                     for(int p=0; p<Hms[s]; p++){
				                     	Hmut[s_tilde*MaxMutNo + p] = Hmut[s*MaxMutNo + p];
				                     
				                     }
				                     
				                     Hms[s_tilde] = Hms[s];
				                     
				                 	for(int a=0; a<NChr; a++){
				                 	
				                 	
				                 		
				                 			CSize[s_tilde*NChr + a] = CSize[s*NChr + a];
				                 			Cprev[s_tilde*NChr + a] = Cprev[s*NChr + a];
				                 			DivCprev[s_tilde*NChr + a] = DivCprev[s*NChr + a];
				                 			CMixprev[s_tilde*NChr + a] = CMixprev[s*NChr + a];
				                 			
				                 			rdelprev[s_tilde*NChr + a] = rdelprev[s*NChr + a];
				                 			rinsprev[s_tilde*NChr + a] = rinsprev[s*NChr + a];
				                 			rtransprev[s_tilde*NChr + a] = rtransprev[s*NChr + a];
				                 			
				                 			rdel[s_tilde*NChr + a] = rdel[s*NChr + a];
				                 			rins[s_tilde*NChr + a] = rins[s*NChr + a];
				                 			rtrans[s_tilde*NChr + a] = rtrans[s*NChr + a];
				                 			
				                 			C[s_tilde*NChr + a] = C[s*NChr + a];
				                 			DivC[s_tilde*NChr + a] = DivC[s*NChr + a];
				                 			CMix[s_tilde*NChr + a] = CMix[s*NChr + a];
				                 						                 		
				                 		for(int d=0; d<NChr; d++){
				                 		   
				                 		    NTprev[s_tilde*NChr*NChr + a*NChr + d] = NTprev[s*NChr*NChr + a*NChr + d];
				                 		    NIprev[s_tilde*NChr*NChr + a*NChr + d] = NIprev[s*NChr*NChr + a*NChr + d];
				                 		    NDprev[s_tilde*NChr*NChr + a*NChr + d] = NDprev[s*NChr*NChr + a*NChr + d];
				                 		    
				                 		    NT[s_tilde*NChr*NChr + a*NChr + d] = NT[s*NChr*NChr + a*NChr + d];
				                 		    NI[s_tilde*NChr*NChr + a*NChr + d] = NI[s*NChr*NChr + a*NChr + d];
				                 		    ND[s_tilde*NChr*NChr + a*NChr + d] = ND[s*NChr*NChr + a*NChr + d];
				                 		    
				                 		    Mprev[s_tilde*NChr*NChr + a*NChr + d] = Mprev[s*NChr*NChr + a*NChr + d];
				                 		    MCprev[s_tilde*NChr*NChr + a*NChr + d] = MCprev[s*NChr*NChr + a*NChr + d];
				                 		    
				                 		    M[s_tilde*NChr*NChr + a*NChr + d] = M[s*NChr*NChr + a*NChr + d];
				                 		    MC[s_tilde*NChr*NChr + a*NChr + d] = MC[s*NChr*NChr + a*NChr + d];
				                 		
				                 		}
				                 		
				                 			
				                 		
				                 		}
				                       s += 1;
				                      }
				                 	
				                 }
				                 
				            
		
      }
      

      }
      else{
      
      }
      }
      
      
}
      

//////////Compute Hamming distance//////////
////////////////////////////////////////////
//cout << "gets to hamming distance: " << gen_clock << std::endl;
MPI_Barrier(MPI_COMM_WORLD);
if(id==0){
cout << "gets to hamming distance: " << my_proc << std::endl;
}
//cout << "Longest_C: " << Longest_C << endl;
MPI_Bcast(&Ncurr, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&MaxCellNo, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&MaxMutNo, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&Hsizeprev, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&Hms, MaxCellNo, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&Hmut, MaxCellNo*MaxMutNo, MPI_INT, 0, MPI_COMM_WORLD);

float local_hamming = 0;
MPI_Barrier(MPI_COMM_WORLD);
if((Ncurr/my_proc)<1){

if(id<Ncurr){


	for(int w=0; w<Hsizeprev; w++){

	int any_v = 0;
	if(Hms[id]<2){
	
	}
	else{
	for(int v=2; v<Hms[id]; v+=2){
		   if((w<3)&&(id<2)){
		   //cout << "w: " << w << " v: " << v << " HMs: " << Hms[d] << " Hmut[d*MaxMutNo + v]: " << Hmut[d*MaxMutNo + v] << endl;
		   }
		   if(Hmut[id*MaxMutNo + v]==w){
		   any_v = 1;
		   
		   for(int z=0; z<Ncurr; z++){

		int anysame = 0;
		if(Hms[z]<2){
		
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(Hmut[id*MaxMutNo + v-1]-Hmut[z*MaxMutNo + m-1]);
				anysame = 1;
				}
			
			
			}
			
			if(anysame==0){
			
			local_hamming += abs(Hmut[id*MaxMutNo + v-1]-0);
			}
			}
			
			
			}
			
			
    }
    


		}
		}
		if(any_v==0){
		
		
		
		 for(int z=0; z<Ncurr; z++){

		int anysame2 = 0;
		if(Hms[z]<2){
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
			//cout << "m: " << m << endl;
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(0-Hmut[z*MaxMutNo + m-1]);
				
				anysame2 = 1;
				}
			
			
			}
			}
			if(anysame2==0){
			
			local_hamming += 0;
			}
		   
		   else{
		   
		  
		   
		   }
	
			
			}
			
			
		
		}
		else{
		
		}

}


}
}
else{

int spacing = Ncurr/60;

if(spacing<1){

if(id==0){

for(int d= 0; d<Ncurr; d++){

	for(int w=0; w<Hsizeprev; w++){

	int any_v = 0;
	if(Hms[d]<2){
	
	}
	else{
	for(int v=2; v<Hms[d]; v+=2){
		   if((w<3)&&(d<2)){
		   //cout << "w: " << w << " v: " << v << " HMs: " << Hms[d] << " Hmut[d*MaxMutNo + v]: " << Hmut[d*MaxMutNo + v] << endl;
		   }
		   if(Hmut[d*MaxMutNo + v]==w){
		   any_v = 1;
		   
		   for(int z=0; z<Ncurr; z++){

		int anysame = 0;
		if(Hms[z]<2){
		
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(Hmut[d*MaxMutNo + v-1]-Hmut[z*MaxMutNo + m-1]);
				anysame = 1;
				}
			
			
			}
			
			if(anysame==0){
			
			local_hamming += abs(Hmut[d*MaxMutNo + v-1]-0);
			}
			}
			
			
			}
			
			
    }
    


		}
		}
		if(any_v==0){
		
		
		
		 for(int z=0; z<Ncurr; z++){

		int anysame2 = 0;
		if(Hms[z]<2){
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
			//cout << "m: " << m << endl;
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(0-Hmut[z*MaxMutNo + m-1]);
				
				anysame2 = 1;
				}
			
			
			}
			}
			if(anysame2==0){
			
			local_hamming += 0;
			}
		   
		   else{
		   
		  
		   
		   }
	
			
			}
			
			
		
		}
		else{
		
		}

}

}

}

}
else{

if(id==0){

for(int d= 0; d<60; d++){

	for(int w=0; w<Hsizeprev; w++){

	int any_v = 0;
	if(Hms[d]<2){
	
	}
	else{
	for(int v=2; v<Hms[d]; v+=2){
		   if((w<3)&&(d<2)){
		   //cout << "w: " << w << " v: " << v << " HMs: " << Hms[d] << " Hmut[d*MaxMutNo + v]: " << Hmut[d*MaxMutNo + v] << endl;
		   }
		   if(Hmut[d*MaxMutNo + v]==w){
		   any_v = 1;
		   
		   for(int z=0; z<Ncurr; z++){

		int anysame = 0;
		if(Hms[z]<2){
		
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(Hmut[d*MaxMutNo + v-1]-Hmut[z*MaxMutNo + m-1]);
				anysame = 1;
				}
			
			
			}
			
			if(anysame==0){
			
			local_hamming += abs(Hmut[d*MaxMutNo + v-1]-0);
			}
			}
			
			
			}
			
			
    }
    


		}
		}
		if(any_v==0){
		
		
		
		 for(int z=0; z<Ncurr; z++){

		int anysame2 = 0;
		if(Hms[z]<2){
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
			//cout << "m: " << m << endl;
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(0-Hmut[z*MaxMutNo + m-1]);
				
				anysame2 = 1;
				}
			
			
			}
			}
			if(anysame2==0){
			
			local_hamming += 0;
			}
		   
		   else{
		   
		  
		   
		   }
	
			
			}
			
			
		
		}
		else{
		
		}

}

}

}
else if(id<floor(Ncurr/60)){
for(int d=60*id; d<60*(id+1); d++){

	for(int w=0; w<Hsizeprev; w++){

	int any_v = 0;
	if(Hms[d]<2){
	
	}
	else{
	for(int v=2; v<Hms[d]; v+=2){
		   if((w<3)&&(d<2)){
		   //cout << "w: " << w << " v: " << v << " HMs: " << Hms[d] << " Hmut[d*MaxMutNo + v]: " << Hmut[d*MaxMutNo + v] << endl;
		   }
		   if(Hmut[d*MaxMutNo + v]==w){
		   any_v = 1;
		   
		   for(int z=0; z<Ncurr; z++){

		int anysame = 0;
		if(Hms[z]<2){
		
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(Hmut[d*MaxMutNo + v-1]-Hmut[z*MaxMutNo + m-1]);
				anysame = 1;
				}
			
			
			}
			
			if(anysame==0){
			
			local_hamming += abs(Hmut[d*MaxMutNo + v-1]-0);
			}
			}
			
			
			}
			
			
    }
    


		}
		}
		if(any_v==0){
		
		
		
		 for(int z=0; z<Ncurr; z++){

		int anysame2 = 0;
		if(Hms[z]<2){
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
			//cout << "m: " << m << endl;
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(0-Hmut[z*MaxMutNo + m-1]);
				
				anysame2 = 1;
				}
			
			
			}
			}
			if(anysame2==0){
			
			local_hamming += 0;
			}
		   
		   else{
		   
		  
		   
		   }
	
			
			}
			
			
		
		}
		else{
		
		}

}

}
}
else if(id==(floor(Ncurr/60) +1)){
for(int d=60*(id-1); d<Ncurr; d++){

	for(int w=0; w<Hsizeprev; w++){

	int any_v = 0;
	if(Hms[d]<2){
	
	}
	else{
	for(int v=2; v<Hms[d]; v+=2){
		   if((w<3)&&(d<2)){
		   //cout << "w: " << w << " v: " << v << " HMs: " << Hms[d] << " Hmut[d*MaxMutNo + v]: " << Hmut[d*MaxMutNo + v] << endl;
		   }
		   if(Hmut[d*MaxMutNo + v]==w){
		   any_v = 1;
		   
		   for(int z=0; z<Ncurr; z++){

		int anysame = 0;
		if(Hms[z]<2){
		
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(Hmut[d*MaxMutNo + v-1]-Hmut[z*MaxMutNo + m-1]);
				anysame = 1;
				}
			
			
			}
			
			if(anysame==0){
			
			local_hamming += abs(Hmut[d*MaxMutNo + v-1]-0);
			}
			}
			
			
			}
			
			
    }
    


		}
		}
		if(any_v==0){
		
		
		
		 for(int z=0; z<Ncurr; z++){

		int anysame2 = 0;
		if(Hms[z]<2){
		}
		else{
			for(int m=2; m<Hms[z]; m+=2){
			
			//cout << "m: " << m << endl;
			
				if(Hmut[z*MaxMutNo + m]==w){
				
				local_hamming += abs(0-Hmut[z*MaxMutNo + m-1]);
				
				anysame2 = 1;
				}
			
			
			}
			}
			if(anysame2==0){
			
			local_hamming += 0;
			}
		   
		   else{
		   
		  
		   
		   }
	
			
			}
			
			
		
		}
		else{
		
		}

}

}
}
}
}

    
MPI_Barrier(MPI_COMM_WORLD);
 float global_hamming;
MPI_Reduce(&local_hamming, &global_hamming, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);        
        if (id==0){
      
      if(extinction==0){   
      
     Print_Hamming[gen_clock] = global_hamming/(Ncurr*(Ncurr-1));
     
     }
     
     else{
     
     Print_Hamming[gen_clock] = -1;
     }
      
      
      double rtrans_mean = 0;
      double rdel_mean = 0;
      double rins_mean = 0;
      //Ncurr
      //Ndead
      //Hdist   
      
      for (int i = 0; i <Ncurr; i++) {

      
	for(int j = 0; j < NChr; j++){

	rtrans_mean += rtrans[i*NChr + j];
	rdel_mean += rdel[i*NChr + j];
	rins_mean += rins[i*NChr + j];
	
	
	}
	}
	
	Print_rtrans[gen_clock] = rtrans_mean/(Ncurr*NChr);
	Print_rdel[gen_clock] = rdel_mean/(Ncurr*NChr);
	Print_rins[gen_clock] = rins_mean/(Ncurr*NChr);
	
	Print_Ndead[gen_clock] = Ndead;
	Print_Ncells[gen_clock] = Ncurr;
	Print_Nborn[gen_clock] = Nborn;
	
	
	
      
   
// char Dirr[100];
//  int Dr;
//       Dr=sprintf(Dirr,"mkdir -p %d/",gen_clock);
//       const int dir= system(Dirr);
// 	char Datvar[100];
//     int Dnvar;
//       Dnvar=sprintf(Datvar,"%d/hd%d.dat",gen_clock,seed);
//     fstream myfileHamming; 
//     myfileHamming.open(Datvar,ios::out);
//     if(myfileHamming.is_open()){
//                                    
//     for(int w=0; w<Ncurr; w++){
//     
//     	for(int p=0; p<Hsizeprev; p++){
// 				    
// 				    if(p==(Hsizeprev-1)){
// 				    
// 				    if(w==(Ncurr-1)){
// 				    myfileHamming <<  Hmut[w*Hsizeprev + p];
// 				    }else{
// 				    myfileHamming <<  Hmut[w*Hsizeprev + p] << '\n';
// 				    }
// 				    }
// 				    
// 				    else{
// 				    myfileHamming <<  Hmut[w*Hsizeprev + p] << '\t';
// 				    }
// 				                     
// 				                     }
// 				
//     }
//     
//     myfileHamming.close();
//     }
//     
//     if(Nborn>0){
//     
//     char Datvar1[100];
//     int Dnvar1;
//       Dnvar1=sprintf(Datvar1,"%d/iborn%d.dat",gen_clock,seed);
//     fstream myfileborn; 
//     myfileborn.open(Datvar1,ios::out);
//     if(myfileborn.is_open()){
//                           
//                           
//     for(int w=0; w<Nborn; w++){
//     
//     if(w==(Nborn-1)){
//     myfileborn << Iborn[w];
//     }
//     else{
//     myfileborn << Iborn[w] << '\n';
//     }
//     
//     }         
//     
//     
//     myfileborn.close();
//     }
//     
//     }
//     
//     if(ThisNDead>0){
//     
//     char Datvar2[100];
//     int Dnvar2;
//       Dnvar2=sprintf(Datvar2,"%d/idead%d.dat",gen_clock,seed);
//     fstream myfiledead; 
//     myfiledead.open(Datvar2,ios::out);
//     if(myfiledead.is_open()){
//                                    
//     for(int w=0; w<ThisNDead; w++){
//     
//     if(w==(ThisNDead-1)){
//     myfiledead << Idead[w];
//     }
//     else{
//     myfiledead << Idead[w] << '\n';
//     }
//     
//     }
//     
//     myfiledead.close();
//     }
// 	}
	}
	
    } // Loop over generations

if (id==0){
    //////////////////////////////////////////////////////
    /////////Print out simulation results ////////////////
    //////////////////////////////////////////////////////
    
    char Dirr[100];
 int Dr;
      Dr=sprintf(Dirr,"mkdir -p Results/%d/",position);
      const int dir= system(Dirr);
	char Datvar[100];
    int Dnvar;
      Dnvar=sprintf(Datvar,"Results/%d/hd%d.dat",position,seed);
    fstream myfileHamming; 
    myfileHamming.open(Datvar,ios::out);
    if(myfileHamming.is_open()){
                                   
    for(int w=1; w<ngen; w++){

				    if(w==(ngen-1)){
				    
				    myfileHamming <<  Print_Hamming[w];
				    }else{
				    myfileHamming <<  Print_Hamming[w] << '\n';
				    }               
				
    }
    
    myfileHamming.close();
    }
	
	char Dat_rins[100];
    int Dn_rins;
      Dn_rins=sprintf(Dat_rins,"Results/%d/rins%d.dat",position,seed);
    fstream myfile_rins; 
    myfile_rins.open(Dat_rins,ios::out);
    if(myfile_rins.is_open()){
                                   
    for(int w=1; w<ngen; w++){

				    if(w==(ngen-1)){
				    
				    myfile_rins <<  Print_rins[w];
				    }else{
				    myfile_rins <<  Print_rins[w] << '\n';
				    }               
				
    }
    
    myfile_rins.close();
    }
	
	char Dat_rdel[100];
    int Dn_rdel;
      Dn_rdel=sprintf(Dat_rdel,"Results/%d/rdel%d.dat",position,seed);
    fstream myfile_rdel; 
    myfile_rdel.open(Dat_rdel,ios::out);
    if(myfile_rdel.is_open()){
                                   
    for(int w=1; w<ngen; w++){

				    if(w==(ngen-1)){
				    
				    myfile_rdel <<  Print_rdel[w];
				    }else{
				    myfile_rdel <<  Print_rdel[w] << '\n';
				    }               
				
    }
    
    myfile_rdel.close();
    }
	
	char Dat_rtran[100];
    int Dn_rtran;
      Dn_rtran=sprintf(Dat_rtran,"Results/%d/rtran%d.dat",position,seed);
    fstream myfile_rtran; 
    myfile_rtran.open(Dat_rtran,ios::out);
    if(myfile_rtran.is_open()){
                                   
    for(int w=1; w<ngen; w++){

				    if(w==(ngen-1)){
				    
				    myfile_rtran <<  Print_rtrans[w];
				    }else{
				    myfile_rtran <<  Print_rtrans[w] << '\n';
				    }               
				
    }
    
    myfile_rtran.close();
    }
    
    char Dat_Ncell[100];
    int Dn_Ncell;
      Dn_Ncell=sprintf(Dat_Ncell,"Results/%d/Ncell%d.dat",position,seed);
    fstream myfile_Ncell; 
    myfile_Ncell.open(Dat_Ncell,ios::out);
    if(myfile_Ncell.is_open()){
                                   
    for(int w=1; w<ngen; w++){

				    if(w==(ngen-1)){
				    
				    myfile_Ncell <<  Print_Ncells[w];
				    }else{
				    myfile_Ncell <<  Print_Ncells[w] << '\n';
				    }               
				
    }
    
    myfile_Ncell.close();
    }
    
    char Dat_Ndead[100];
    int Dn_Ndead;
      Dn_Ndead=sprintf(Dat_Ndead,"Results/%d/Ndead%d.dat",position,seed);
    fstream myfile_Ndead; 
    myfile_Ndead.open(Dat_Ndead,ios::out);
    if(myfile_Ndead.is_open()){
                                   
    for(int w=1; w<ngen; w++){

				    if(w==(ngen-1)){
				    
				    myfile_Ndead <<  Print_Ndead[w];
				    }else{
				    myfile_Ndead <<  Print_Ndead[w] << '\n';
				    }               
				
    }
    
    myfile_Ndead.close();
    }
    
    char Dat_Nborn[100];
    int Dn_Nborn;
      Dn_Nborn=sprintf(Dat_Nborn,"Results/%d/Nborn%d.dat",position,seed);
    fstream myfile_Nborn; 
    myfile_Nborn.open(Dat_Nborn,ios::out);
    if(myfile_Nborn.is_open()){
                                   
    for(int w=1; w<ngen; w++){

				    if(w==(ngen-1)){
				    
				    myfile_Nborn <<  Print_Nborn[w];
				    }else{
				    myfile_Nborn <<  Print_Nborn[w] << '\n';
				    }               
				
    }
    
    myfile_Nborn.close();
    }
    
    
    ///////
    //Time series and cells left
    ///////


   cout << "Did it complete generations? " << endl;


    
	cout << "position: " << position << endl;
    
    ///////////////////////////
    // Clear containers ///////
    ///////////////////////////
    
    for (int i = 0; i < MaxCellNo; i++) {
    
	  
      Hmut[i] = 0;
	  Nmut[i] = 0;
	  n_chr[i] = 0;
      Iborn[i] = 0;
      Idead[i] = 0;
      
      //step containers
      n_chrprev[i] = 0;
      double Max_ChrSize[NChr];
	
      for(int j = 0; j < NChr; j++){
	
		C[i*NChr+j] = 0;
		CMix[i*NChr+j] = 0;
		CSize[i*NChr+j] = 0;
	         
		Cprev[i*NChr+j]=0;
		CMixprev[i*NChr+j]=0;
		Max_ChrSize[j] = 0;		
	
		DivC[i*NChr+j]=0;
		DivCprev[i*NChr+j]=0;
	
      }
      
      CmaxL[i] = 0;
      CminL[i] = 0;
            
      for(int j = 0; j < NChr; j++){
      for(int k = 0; k < NChr; k++){
	
		NT[i*NChr*NChr+j*NChr+k]=0;
		NI[i*NChr*NChr+j*NChr+k]=0;
		ND[i*NChr*NChr+j*NChr+k]=0;
		NTprev[i*NChr*NChr+j*NChr+k]=0;
	    NIprev[i*NChr*NChr+j*NChr+k]=0;
	    NDprev[i*NChr*NChr+j*NChr+k]=0;
	  
	    M[i*NChr*NChr+j*NChr+k]=0;
	    MC[i*NChr*NChr+j*NChr+k]=0;
		Mprev[i*NChr*NChr+j*NChr+k]=0;
		MCprev[i*NChr*NChr+j*NChr+k]=0;
	  		
		 }
		 
	}
    
    
    
    }
}
gsl_rng_free (r);
 //stream->free_sprng();		/* free memory used to store stream state  */
  
MPI_Finalize();
  
  return(0); 
}

int read_constnts( int id, const string& filename, int &NChr, int &Npop, int &ngen, double &lowB, double &chrlength, double &minchrlen, double &Curvemax, double &minchr, double &mu_b, double &minSV, double &p_birth, double &p_death, int &Nparticles, int &MaxCellNo, int &MaxMutNo, double &Lchr_loss, double &Uchr_loss, double &Lp_tran, double &Up_tran, double &Lmu_ki, double &Umu_ki, double &Lmu_kd, double &Umu_kd,  double &Lfitness, double &Ufitness, double &Lmaxchr, double &Umaxchr, double &Lgnmdou, double &Ugnmdou ){

if (id==0){
cout << "open file" << endl;
  ifstream infile (filename.c_str());
  int counter = 0; 
  if (infile.is_open()){

    string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

		vector<string> split;
      string buf; 
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

		if(counter==0){
      
    
     NChr = atoi( split[0].c_str() ) ;
     Npop = atoi( split[1].c_str() ) ;
     ngen = atoi( split[2].c_str() ) ;
     lowB  = atof( split[3].c_str() ) ;
     chrlength = atof( split[4].c_str() ) ;
     minchrlen = atof( split[5].c_str() ) ;
     Curvemax =  atof( split[6].c_str() ) ;
     minchr = atof( split[7].c_str() ) ;
     mu_b = atof( split[8].c_str() ) ;
     minSV = atof( split[9].c_str() ) ;
     p_birth = atof( split[10].c_str() ) ;
     p_death = atof( split[11].c_str() ) ;
     Nparticles = atoi( split[12].c_str() ) ;
     MaxCellNo = atoi( split[13].c_str() ) ;
     MaxMutNo = atoi( split[14].c_str() ) ;
     
     }
     else{
      
      cout << "Line: " << split[0] << ' ' << split[1] << endl;
      Lchr_loss = atof( split[0].c_str() ) ; 
      Uchr_loss = atof( split[1].c_str() ) ; 
      Lp_tran = atof( split[2].c_str() ) ;  
      Up_tran = atof( split[3].c_str() ) ; 
      Lmu_ki = atof( split[4].c_str() ) ;  
      Umu_ki = atof( split[5].c_str() ) ;  
      Lmu_kd = atof( split[6].c_str() ) ;  
      Umu_kd = atof( split[7].c_str() ) ;  
      Lfitness = atof( split[8].c_str() ) ;  
      Ufitness = atof( split[9].c_str() ) ;  
      Lmaxchr = atof( split[10].c_str() ) ;  
      Umaxchr = atof( split[11].c_str() ) ; 
      Lgnmdou = atof( split[12].c_str() ) ;  
      Ugnmdou = atof( split[13].c_str() ) ; 
     
     }
    
      counter++;
    }
      
  }else{
    cerr << "Error: open of constant parameters file unsuccessful: " <<  filename << endl;
  }
  
  return counter;
  }
}














