

#include <iostream>
#include <fstream>
#include <cstring>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/graph_types.h>
#include <seqan/index.h>
#include "../vlmm.h"
#include "../vlmm_scoreTree.h"
#include "../vlmm_Classification.h"
#include <string>
#include <vector>


using namespace std;
using namespace seqan;

const char * EmptyString = "I'm Empty !";
const unsigned EmptyInt = UINT_MAX;

int getNfam(String<char> nameFam){
	
	ifstream family;
	family.open(toCString(nameFam),ifstream::in);
	string temp;
	int count=0;
	
	
	if(family.good()){
		while(getline(family,temp)){
			if(temp[0]=='>'){
				count++;
			}
		}
	}
	else{
		cout << "bad";
		cout << endl;
	}
	family.close();
	return(count);
	
}

StringSet<String<AminoAcid>,Owner < > >
createStringSetFromFile(String<char> Filename){
	
	StringSet<String<AminoAcid>,Owner< > > res;
	int i = 0;
	int num= getNfam(Filename);
	resize(res,num);
	fstream setFile;
	setFile.open(toCString(Filename),ios_base::in|ios_base::binary);
	String<AminoAcid> est;
	while (!setFile.eof()) {
		read(setFile,est,Fasta());
		res[i]=est;
		i++;
	}
	
	return(res);
}




//Concatene les trois chaines pour avoir le chemin complet de name:
String<char> appendPath(String<char> current,String<char> fFolder,String<char> name)
{	
	String<char> res("");
	current+=fFolder;
	current+=name;
	res+=current;
	return(res);
	
}



void gnuBuild(int time,String<char> familyName,
	      String<char> tp,String<char> cp,String<char> gp,String<char> ap,String<char> dp){
  String<char> gnuPath("/home/gl/Development/seqan-trunk/projects/library/apps/vlmm/");
  String<char> Results2("Results/Parametres/");
  String<char> tool("/");
  String<char> gnuName;
  String<char> cd("");
  cd+= gnuPath;
  cd+=Results2;
  cd+=familyName;
  if(time==1){
  gnuName="gnu.time.";
  }
  else{
  gnuName="gnu.";
  }
  gnuPath+=Results2;
  gnuPath+=familyName;
  gnuPath+=tool;
  gnuPath+=gnuName;
  gnuPath+=familyName;

  String<char> output("");
  String<char> ext(".eps");
  output+=gnuPath;
  output+=ext;
  String<char> to("'");
  
  //cout<< toCString(gnuPath)<<endl;
  
  //cout<< toCString(output)<<endl;
 
  ofstream gnuFile;
  gnuFile.open(toCString(gnuPath),ios_base::out);
  gnuFile<< "cd "<<toCString(to)<< toCString(cd)<<toCString(to)<<endl;
  gnuFile<< "set terminal postscript eps color"<<endl;
  gnuFile<< "set output " <<toCString(to) <<toCString(output)<<toCString(to) << endl;
  gnuFile<< "set key right bottom " << endl;
  gnuFile<< "plot " <<toCString(to)<< toCString(tp)<<toCString(to) << " with lines, ";
  gnuFile <<toCString(to)<< toCString(cp)<<toCString(to) << " with lines, ";
  gnuFile <<toCString(to)<< toCString(gp)<<toCString(to) << " with lines, ";
  gnuFile <<toCString(to)<< toCString(ap)<<toCString(to) << " with lines, ";
  gnuFile <<toCString(to)<< toCString(dp)<<toCString(to) << " with lines";

}

int main(int argc,char ** argv){
	//Directories:
	String<char> currentPath("/home/gl/Development/seqan-trunk/projects/library/apps/vlmm/");
	String<char> dataBasePath("Pfam/Pfam_raw_seq/");
	String<char> BaseSwiss("Pfam/swissprot_v33/");
	String<char> BaseTrain("Base/Train/");
	String<char> Results("Results/");
	String<char> Results2("Results/Parametres/");
	//Names:
	
	String<char> swissName("sprot33.dat.fasta");
	
	String<char> familyName(argv[1]);
	String<char> swissResultName("like.");
	swissResultName+=familyName;
	String<char> familyPath;
	familyPath+=Results2;
	String<char> Path("/");
	familyPath+=familyName;
	familyPath+=Path;
	String<char> train("train.");
	train+=familyName;
	
	//Results names:
	//String<char> time("time.")
	String<char> tt("time.");
	String<char> tc("time.");
	String<char> tg("time.");
	String<char> td("time.");
	String<char> ta("time.");
	String<char> thresholdN("threshold.");
	String<char> cutoffN("cutoff.");
	String<char> gammaN("gamma.");
	String<char> depthN("depth.");
	String<char> alphaN("alpha.");
	thresholdN+=familyName;
	cutoffN+=familyName;
	gammaN+=familyName;
	depthN+=familyName;
	alphaN+=familyName;
	tt+=thresholdN;
	tc+=cutoffN;
	tg+=gammaN;
	td+=depthN;
	ta+=alphaN;
	
	//Paths:
	String<char> trainPath;
	trainPath=appendPath(currentPath,BaseTrain,train);
	String<char> swissPath;
	swissPath=appendPath(currentPath,BaseSwiss,swissName);
	//String<char> swissPath;
	//SwissPath=appendPath(currentPath,BaseSwiss,swissName);
	String<char> swissResultsPath;
	swissResultsPath=appendPath(currentPath,Results,swissResultName);
	
	String<char> ttP=appendPath(currentPath,familyPath,tt);
	String<char> tcP=appendPath(currentPath,familyPath,tc);
	String<char> tgP=appendPath(currentPath,familyPath,tg);
	String<char> tdP=appendPath(currentPath,familyPath,td);
	String<char> taP=appendPath(currentPath,familyPath,ta);
	
	String<char> thresholdPath;
	thresholdPath=appendPath(currentPath,familyPath,thresholdN);
	
	String<char> cutoffPath;
	cutoffPath=appendPath(currentPath,familyPath,cutoffN);
	
	String<char> gammaPath;
	gammaPath=appendPath(currentPath,familyPath,gammaN);
	
	String<char> depthPath;
	depthPath=appendPath(currentPath,familyPath,depthN);
	
	String<char> alphaPath;
	alphaPath=appendPath(currentPath,familyPath,alphaN);
	
	String<char> resName("background.pst");
	String<char> backgroundFile;
	backgroundFile=appendPath(currentPath,Results,resName);

//PST parameters:
	float threshold,cutoff,ymin,alpha,d,max;
	max=0;
	float gthr=0;
	float gcut=0;
	float gy=0;
	float galpha=0;
	float gd=0;
	double T=10;
	double gT=0;
	//gthr=0.00001;
	//gy=0.00001;
	
	bool nl=false;
	bool ll=true;
	double sens;
	vector<double> p;
	vector<double> param;
	vector<clock_t> time;
	
	//Defining structures:
	typedef String<AminoAcid> TText;	
	typedef StringSet< TText, Owner< > > TMultiProtein;
	typedef Index<TMultiProtein, Index_Wotd< > > TWotd;
	
	String<String<AminoAcid> > familySeqs,backgroundSeqs;
	String<String<char> > familyIds,backgroundIds;
	//family is a file name of the training set.
	//lit le fichier fasta family et renvoi les sequences et les ids.
	createInputString(trainPath,familySeqs,familyIds);
	//Meme chose avec la base de données.
	//Ici background est swissProt
	createInputString(swissPath,backgroundSeqs,backgroundIds,familyIds);
	
//Read Background Model:
typedef Graph<Automaton<AminoAcid,String<AminoAcid>,
	WordGraph < VLMM < PST > > > > Vlmm;
	Vlmm backgroundModel;
	cout << "ici" <<endl;
	//Lecture du modèle de Background:
	ifstream background;
	background.open(toCString(backgroundFile),ios_base::in|ios_base::binary);
	cout << "ici2" << endl;
	String<char> entry;
	_scanNextEntry(background,entry);
	_scanNextEntry(background,entry);
	_scanNextEntry(background,entry);
	readGraph(backgroundModel,background);
	cout << "ici3" <<endl;
	TMultiProtein w;
	w=createStringSetFromFile(trainPath);
	TWotd Tesa(w);
	
	cout<< "ici4" <<endl;
	
	cout << "Test Parameters for:";
	cout << toCString(familyName)<<endl;
	cout<< "Step 1 : Threshold !!!" <<endl;
	int cpt=0;
	
	/*Prendre une sorte de moyenne pour la meilleur valeur de
	 * temps des parametres alpha depth et cutoff
	 */
	 typedef Graph<Automaton<AminoAcid,String<AminoAcid>,
		WordGraph<VLMM<ScoreTree> > > > TScoreTree;
	
		//-------------------Cutoff-------------------
	cout<< "Step 3 : Cutoff !!!" <<endl;
	for (cutoff=0; cutoff<=10; cutoff+=0.5) {
	clock_t  temps_initial=clock();	
	
		
		Vlmm pst;
		//cutoff=1.05;
		threshold=0.00000001;
		ymin=0.00000001;
		alpha=0;
		d=20;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		

	//Creation du ScoreTree:
		
	
		TScoreTree scoreTree;
		buildUnionScoreTree(scoreTree,backgroundModel,pst);
		//Equivalence number criterion:
		
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		
									
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);							
									
									
		p.push_back(cutoff);
		param.push_back(sens);
		
		if(sens>=max){
			max=sens;
			gcut=cutoff;
			
		}
		
	}
	
	
	
	//Save results (gnuplot):
	ofstream testcutoff;
	testcutoff.open(toCString(cutoffPath),ios::out);
	
	for(int i=0;i<p.size();i++){
		testcutoff << p[i];
		testcutoff << "\t";
		testcutoff << param[i];
		testcutoff << endl;
		
		
	}
	
	ofstream fc;
	fc.open(toCString(tcP),ios::out);
	for(int i=0;i<p.size();i++){
		fc << p[i];
		fc << "\t";
		fc << time[i];
		fc << endl;
		
	}
	
	
	
	//Vider les vecteurs p et param:
	p.clear();
	param.clear();
	time.clear();
	int c=0;
	max=0;
	
	
	
	//-----------------Depth:------------------
	cout<< "Step 5 : Depth !!!" <<endl;
	for (d=6; d<=20; d+=1) {
	  clock_t  temps_initial=clock();
	
		
		Vlmm pst;
		cutoff=gcut;
		threshold=0.00000001;
		ymin=0.00000001;
		alpha=0;
		//d=gd;
	cout << "1";	
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
	cout << "2";
		//Equivalence number criterion:
		//Creation du ScoreTree:

	TScoreTree scoreTree;
	cout << "3";
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
	cout << "4";
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
								familyIds,swissResultsPath,nl,ll,d,sens);
	cout << "5";
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);
		
		p.push_back(d);
		param.push_back(sens);
		//>=?
		if(sens>max){
			max=sens;
			gd=d;
			
		}
		if(gd==0){gd=10;}
		
	}
	
	//Save results (gnuplot):
	ofstream testd;
	testd.open(toCString(depthPath),ios::out);
	
	for(int i=0;i<p.size();i++){
		testd << p[i];
		testd << "\t";
		testd << param[i];
		testd << endl;
		
		
	}
	
	ofstream fd;
	fd.open(toCString(tdP),ios::out);
	for(int i=0;i<p.size();i++){
		fd << p[i];
		fd << "\t";
		fd << time[i];
		fd << endl;
		
	}
		
			//Vider les vecteurs p et param:
	p.clear();
	param.clear();
	time.clear();	
	max=0;
		

  	//-----------------Alpha:------------------
	cout<< "Step 4 : Alpha !!!" <<endl;
	for (alpha=0;alpha<=5; alpha+=0.2) {
	  clock_t  temps_initial=clock();
	
		
		Vlmm pst;
		cutoff=gcut;
		threshold=0.00000001;
		ymin=0.00000001;
		//alpha=0;
		d=gd;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:

	
	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);							
									
									
		p.push_back(alpha);
		param.push_back(sens);
		
		if(sens>=max){
			max=sens;
			galpha=alpha;
			
		}
		//c++;
		//cout<<c<<endl;
	}
	
	//Save results (gnuplot):
	ofstream testalpha;
	testalpha.open(toCString(alphaPath),ios::out);
	
	for(int i=0;i<p.size();i++){
		testalpha << p[i];
		testalpha << "\t";
		testalpha << param[i];
		testalpha << endl;
		
		
	}
	
	ofstream fa;
	fa.open(toCString(taP),ios::out);
	for(int i=0;i<p.size();i++){
		fa << p[i];
		fa << "\t";
		fa << time[i];
		fa << endl;
		
	}
	
	
	
	
	
	//Vider les vecteurs p et param:
	p.clear();
	param.clear();
	time.clear();
	max=0;
	
		
	
				//---------------ymin:--------------
	cout<< "Step 2 : GammaMin !!!" <<endl;
	for (ymin=0.00000001; ymin<0.01; ymin*=2.5) {
	clock_t  temps_initial=clock();	
	
		
		Vlmm pst;
		cutoff=gcut;
		threshold=0.00000001;
		//ymin=0.000000001;
		alpha=galpha;
		d=gd;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:

	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		
									
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);							
									
									
		p.push_back(8+log10(ymin));
		param.push_back(sens);
		
		if(sens>=max){
			max=sens;
			gy=ymin;
			
		}
		
		//cout<<c<<endl;
	}
	
	//Save results (gnuplot):
	ofstream testgamma;
	testgamma.open(toCString(gammaPath),ios::out);
	
	for(int i=0;i<p.size();i++){
		testgamma << p[i];
		testgamma << "\t";
		testgamma << param[i];
		testgamma << endl;
		
		
	}
	
	ofstream fg;
	fg.open(toCString(tgP),ios::out);
	for(int i=0;i<p.size();i++){
		fg << p[i];
		fg << "\t";
		fg << time[i];
		fg << endl;
		
	}
	
	
	
	//Vider les vecteurs p et param:
	p.clear();
	param.clear();
	time.clear();
	//c=0;
	max=0;
	

	
	cout << "Pmin" <<endl;
	//-------------------Threshold:---------------------
	for (threshold=0.00000001; threshold<=0.01; threshold*=2.5) {
		clock_t  temps_initial=clock();
	  
		
		Vlmm pst;
		cutoff=gcut;
		ymin=gy;
		alpha=galpha;
		d=gd;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:

	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);
		
		
		p.push_back(8+log10(threshold));
		param.push_back(sens);
		
		if(sens>=max){
			max=sens;
			gthr=threshold;

		}
		cpt++;
		
	}
	cout << "CPT: "<< cpt << endl;
	//Save results (gnuplot):
	ofstream testThreshold;
	testThreshold.open(toCString(thresholdPath),ios::out);
	
	for(int i=0;i<p.size();i++){
		testThreshold << p[i];
		testThreshold << "\t";
		testThreshold << param[i];
		testThreshold << endl;
		
	}
	
	ofstream ft;
	ft.open(toCString(ttP),ios::out);
	for(int i=0;i<p.size();i++){
		ft << p[i];
		ft << "\t";
		ft << time[i];
		ft << endl;
		
	}
	
	

	//Vider les vecteurs p et param:
	p.clear();
	param.clear();
	time.clear();
	max=0;
	
	
	
	
	
	

//---------------------------------------------------------------------------------------------------------------------	
//Calculate the sensibility for the first set of parameters found.	
	float goodthr=gthr;
	float goody=gy;
	float goodalpha=galpha;
	float goodcut=gcut;
	float goodd=gd;
	
		
		Vlmm pst;
	buildPST(Tesa,pst,goodcut,goodthr,goody,goodalpha,goodd);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:
	

	
	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
	
	float sensStep1=sens;
	
	
//---------------------------------------------------------------------------------------------------------------------	
	
	
		//-------------------Cutoff-------------------
	cout<< "Step 3 : Cutoff !!!" <<endl;
	for (cutoff=0; cutoff<=10; cutoff+=0.5) {
	clock_t  temps_initial=clock();	
	  
		
		Vlmm pst;
		//cutoff=1.05;
		threshold=gthr;
		ymin=gy;
		alpha=galpha;
		d=gd;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:
	
	
	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		
									
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);							
									
									
		p.push_back(cutoff);
		param.push_back(sens);
		
		if(sens>=max){
			max=sens;
			gcut=cutoff;
			
		}
		
	}
	


	
	cout<< "Step 4 : Alpha !!!" <<endl;
	for (alpha=0;alpha<=7; alpha+=1) {
	  clock_t  temps_initial=clock();
		
		Vlmm pst;
		cutoff=gcut;
		threshold=gthr;
		ymin=gy;
		//alpha=0;
		d=gd;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:

	
	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);							
									
									
		p.push_back(alpha);
		param.push_back(sens);
		
		if(sens>=max){
			max=sens;
			galpha=alpha;
			
		}
		//c++;
		//cout<<c<<endl;
	}
	
	for (threshold=0.00000001; threshold<=0.01; threshold*=3) {
		clock_t  temps_initial=clock();
	
		
		Vlmm pst;
		cutoff=gcut;
		ymin=gy;
		alpha=galpha;
		d=gd;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:
	
	
	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);
		
		
		p.push_back(8+log10(threshold));
		param.push_back(sens);
		
		if(sens>=max){
			max=sens;
			gthr=threshold;

		}
		//cpt++;
		
	}
	
	
	cout<< "Step 2 : GammaMin !!!" <<endl;
	for (ymin=0.0000001; ymin<0.01; ymin*=3) {
	clock_t  temps_initial=clock();	
	  
		
		Vlmm pst;
		cutoff=gcut;
		threshold=gthr;
		//ymin=0.000000001;
		alpha=galpha;
		d=gd;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:
	
	
	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		
									
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);							
									
									
		p.push_back(8+log10(ymin));
		param.push_back(sens);
		
		if(sens>=max){
			max=sens;
			gy=ymin;
			
		}
	}
	
		cout<< "Step 5 : Depth !!!" <<endl;
	for (d=6; d<=20; d+=2) {
	  clock_t  temps_initial=clock();
		
		
		Vlmm pst;
		cutoff=gcut;
		threshold=gthr;
		ymin=gy;
		alpha=galpha;
		//d=gd;
		
		buildPST(Tesa,pst,cutoff,threshold,ymin,alpha,d);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:

	
	TScoreTree scoreTree;
	buildUnionScoreTree(scoreTree,backgroundModel,pst);
		
		equivalenceNumberCriterion2(scoreTree,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
		clock_t temps_final=clock();
		clock_t temps_cpu=(temps_final-temps_initial)/CLOCKS_PER_SEC *1000;
		time.push_back(temps_cpu);
		
		p.push_back(d);
		param.push_back(sens);
		//>=?
		if(sens>max){
			max=sens;
			gd=d;
			
		}
		if(gd==0){gd=10;}
		
	}
	

//Verify results:
	
	
		
		Vlmm pst1;
	buildPST(Tesa,pst1,gcut,gthr,gy,galpha,gd);
		
		//Equivalence number criterion:
		//Creation du ScoreTree:

	
	TScoreTree scoreTree1;
	buildUnionScoreTree(scoreTree1,backgroundModel,pst1);
		
		equivalenceNumberCriterion2(scoreTree1,backgroundSeqs,backgroundIds,familySeqs,
									familyIds,swissResultsPath,nl,ll,d,sens);
	


	float sensStep2=sens;
	
	float bestthr,besty,bestalpha,bestcut,bestd;
	if(sensStep1>sensStep2){
	  bestthr=goodthr;
	  besty=goody;
	  bestalpha=goodalpha;
	  bestcut=goodcut;
	  bestd=goodd;
	}
	else{
	  bestthr=gthr;
	  besty=gy;
	  bestalpha=galpha;
	  bestcut=gcut;
	  bestd=gd;
	  
	}

//Write the solution:

	String<char> bestParam("best.");
	bestParam+=familyName;
	String<char> bestPath;
	bestPath=appendPath(currentPath,familyPath,bestParam);
	ofstream best;
	if(!best.good()){
	  cout<<"error"<<endl;
	}
	best.open(toCString(bestPath),ios::out);
	best<<bestthr<<endl;
	best<<besty<<endl;
	best<<bestcut<<endl;
	best<<bestalpha<<endl;
	best<<bestd<<endl;
	//best<<gT<<endl;

	cout<< "End of :";
	cout<< toCString(familyName) << endl;
	
	
	//
	//gnuBuild(0,familyName,thresholdPath,cutoffPath,gammaPath,alphaPath,depthPath);	
	gnuBuild(0,familyName,thresholdN,cutoffN,gammaN,alphaN,depthN);	
	gnuBuild(1,familyName,tt,tc,tg,ta,td);
	
	return 0;
}
						
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


