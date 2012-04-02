#ifndef SEQAN_HEADER_VLMM_CLASSIFICATION_H
#define SEQAN_HEADER_VLMM_CLASSIFICATION_H



namespace SEQAN_NAMESPACE_MAIN
{

unsigned findPosition(String<double> probVector,double randomNumber){
double sum=0;
SEQAN_ASSERT(length(probVector) >0)
for(unsigned i = 0;i<length(probVector);++i){
	sum += probVector[i];
	if(randomNumber<= sum)
		return i;
}
return (length(probVector)-1);
}

void getRandomIndices(String<unsigned> &indices,unsigned size){
MTRand_closed drand((unsigned)time(NULL));
String<double> boundaries;
String<bool> control;
resize(boundaries,size);
resize(control,size);
for(unsigned i=0;i<size;++i){
	boundaries[i] = (double)1/size;
	control[i] =false;
}

//exit(0);
unsigned index = 0;
while(index != length(indices)){
	double number = (double)drand();
	unsigned pos =  findPosition(boundaries,number);
	if(!control[pos]){
		control[pos] = true;
		++index;
	}
}
//get the indices ordered
index=0;
for(unsigned i=0;i<length(control);++i){
	if(control[i]){
		indices[index] = i+1;
		++index;
	}
	if(index == length(indices))
		break;
}
}

// use the iso point method for classification
// train only on the indices provided and do the classification on all sequences
template<typename TAlphabet, typename TCargo, typename TVLMMSpec,typename TParams>
inline void
equivalenceNumberCriterion(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > & vlmm,
		 TParams &parameters,
		 String<String<TAlphabet> > &database,
		 String<String<char> > &databaseIDs,
		 String<String<TAlphabet> > &trainingSet,
		 String<String<char> > &trainingIDs,
		 String<unsigned> &trainingIndices,
		 String<pair<double,unsigned > > &falsePos,
		 String<pair<double,unsigned > > &falseNeg)
{
// train the vlmm on the training set
typedef StringSet< String<TAlphabet>, Owner<> > TStringSet;
Index<TStringSet, Index_ESA<> > esa;
createIndexFromIndicesOfDatabase(esa,trainingSet,trainingIndices);
build(esa,vlmm,parameters);
// build two priority queues one for the likelihoods/scores one the trainingIndices one for the rest of the database
std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
estimateNormalizedLikelihoodOnSequences(vlmm,PQDatabase,database);
estimateNormalizedLikelihoodOnSequences(vlmm,PQTraining,trainingSet);
// estimate the equivalence number criterion
unsigned falsePositives=0,falseNegatives=length(trainingSet);
// output statistics about it
while(falsePositives != falseNegatives){
	if(PQDatabase.top().first <= PQTraining.top().first){
		--falseNegatives;
		PQTraining.pop();
	}
	else{
		// this is a false Positive which needs to be recognized
		++falsePositives;
		appendValue(falsePos,PQDatabase.top());
		PQDatabase.pop();
	}
}
// now we have estimated the equivalence number criterion
// all values left in the priority queue PQTraining are false negatives as they
// are not detected, so put them in falseNeg
while(!PQTraining.empty())
{
	appendValue(falseNeg,PQTraining.top());
	PQTraining.pop();
}
SEQAN_ASSERT(length(falseNeg) == length(falsePos))
}

// use the iso point method for classification
// train only on the indices provided and do the classification on all sequences
template<typename TAlphabet, typename TCargo>
inline void
equivalenceNumberCriterion(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM <ScoreTree > > > > & scoreTree,
		 String<String<TAlphabet> > &database,
		 String<String<char> > &databaseIDs,
		 String<String<TAlphabet> > &trainingSet,
		 String<String<char> > &trainingIDs,
		 String<char> &output,
		 bool ns,
		 bool ls,
		 unsigned maxDepth)
{
// build two priority queues one for the likelihoods/scores one the trainingIndices one for the rest of the database

fstream outFile;
if(ns){
	std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
	estimateNormalizedScoreOnSequences(scoreTree,PQTraining,trainingSet);
	estimateNormalizedScoreOnSequences(scoreTree,PQDatabase,database);
	openAndAppend(output,".ns.enc",outFile);
	printEquivalenceNumberCriterion(outFile,PQTraining,PQDatabase,databaseIDs,trainingIDs,length(database),length(trainingSet));
	outFile.close();
	std::cout << "wrote classification with normalized score for sequences into: "<<output<<".ns.enc"<<endl;
}
if(ls){
	std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
	//cout << "LSTraining: ";
	estimateLocalScoreOnSequences(scoreTree,PQTraining,trainingSet,maxDepth);
	//cout << "LSSwiss: ";
	estimateLocalScoreOnSequences(scoreTree,PQDatabase,database,maxDepth);
	
	openAndAppend(output,".ls.enc",outFile);
	printEquivalenceNumberCriterion(outFile,PQTraining,PQDatabase,databaseIDs,trainingIDs,length(database),length(trainingSet));
	outFile.close();
	std::cout << "wrote classification with local score for sequences into: "<<output<<".ls.enc"<<endl;
}

}


// use the iso point method for classification
// train only on the indices provided and do the classification on all sequences
template<typename TAlphabet, typename TCargo,typename TVLMMSpec>
inline void
equivalenceNumberCriterion(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM <TVLMMSpec > > > > & vlmm,
		 String<String<TAlphabet> > &database,
		 String<String<char> > &databaseIDs,
		 String<String<TAlphabet> > &trainingSet,
		 String<String<char> > &trainingIDs,
		 String<char> &output,
		 bool nl,
		 bool wl,
		 unsigned windowSize)
{
// build two priority queues one for the likelihoods/scores one the trainingIndices one for the rest of the database


if(nl){
	fstream outFile;
	std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
	estimateNormalizedLikelihoodOnSequences(vlmm,PQDatabase,database);
	estimateNormalizedLikelihoodOnSequences(vlmm,PQTraining,trainingSet);
	openAndAppend(output,".nl.enc",outFile);
	printEquivalenceNumberCriterion(outFile,PQTraining,PQDatabase,databaseIDs,trainingIDs,length(database),length(trainingSet));
	outFile.close();
	std::cout << "wrote classification with normalized likelihood for sequences into: "<<output<<".nl.enc"<<endl;
}
if(wl){
	fstream outFile;
	std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
	estimateLikelihoodWindowOnSequences(vlmm,PQTraining,trainingSet,windowSize);
	estimateLikelihoodWindowOnSequences(vlmm,PQDatabase,database,windowSize);
	openAndAppend(output,".wl.enc",outFile);
	printEquivalenceNumberCriterion(outFile,PQTraining,PQDatabase,databaseIDs,trainingIDs,length(database),length(trainingSet));
	outFile.close();
	std::cout << "wrote classification with window likelihood for sequences into: "<<output<<".wl.enc"<<endl;
}

}


template<typename TAlphabet, typename TCargo, typename TVLMMSpec,typename TParams>
inline void classify(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > & vlmm,
		 TParams &parameters,
		 String<char> &databaseFile,
		 String<char> &trainingFile,
		 String<char> &outputFile,
		 String<unsigned> &indices,
		 double trainingRatio)
{
	String<String<TAlphabet> > database,trainingSet;
	String<String<char> > databaseIDs,trainingIDs;
	createInputString(trainingFile,trainingSet,trainingIDs);
	// this function puts only sequences in the string database, which are not trainingIDs
	createInputString(databaseFile,database,databaseIDs,trainingIDs);
	String<pair<double,unsigned> > falsePos,falseNeg;
	// how many values of the trainingSet are used for the training
	/*unsigned numberIndices = (unsigned)floor(trainingRatio * length(trainingSet));
	String<unsigned> indices;
	resize(indices,numberIndices);
	getRandomIndices(indices,length(trainingSet));*/
	equivalenceNumberCriterion(vlmm,parameters,database,databaseIDs,trainingSet,trainingIDs,indices,falsePos,falseNeg);
	// print out the 
	printEquivalenceNumberCriterion( outputFile, falsePos, falseNeg,databaseIDs,trainingIDs,length(database),length(trainingSet));
}

// the indices should be relative to the database given +1
// i.e pos 1 is not 0 like in the array
template<typename TAlphabet,typename TType,typename TIndex>
void createIndexFromIndicesOfDatabase( Index<TType ,TIndex> &esa,
									   String<String<TAlphabet> >  &database,
									   String<unsigned> &indices)
{
	unsigned index = 0;
	for(unsigned i=0;i<length(database);++i){
		if(indices[index] == i+1){
			appendValue(indexText(esa),database[i]);
			++index;
			if(index == length(indices))
				break;
		}
	}
	SEQAN_ASSERT(index == length(indices))
}

// the indices should be relative to the database given +1
// i.e pos 1 is not 0 like in the array
template<typename TAlphabet,typename TType,typename TIndex>
void createIndexWithoutIndicesOfTraining( Index<TType ,TIndex> &esa,
									   String<String<TAlphabet> >  &database,
									   String<unsigned> &trainingIndices)
{
	unsigned index = 0;
	for(unsigned i=0;i<length(database);++i){
		if(trainingIndices[index] == i+1){
			appendValue(indexText(esa),database[i]);
			++index;
			if(index == length(trainingIndices))
				break;
		}
	}
	SEQAN_ASSERT(index == length(trainingIndices))
}




/*
* Print the classification statistics
*/
// calculate the equivalence number criterion and output that into a file
void printEquivalenceNumberCriterion(fstream &outputFile,
									 std::priority_queue<pair<double,unsigned> > &PQTraining,
									 std::priority_queue<pair<double,unsigned> > &PQDatabase,
									 String<String<char> > databaseIDs,
									 String<String<char> > trainingIDs,
									 unsigned /*lengthDatabase*/,
									 unsigned lengthTrainingSet)
{
	outputFile.precision(4);
	outputFile << "ID\tScore\tFP\t#FP\t#FN\n";

// estimate the equivalence number criterion
unsigned falsePositives=0,falseNegatives=lengthTrainingSet;
// output statistics about it
while(falsePositives != falseNegatives){
	if(PQDatabase.top().first <= PQTraining.top().first){
		--falseNegatives;
		outputFile <<  trainingIDs[PQTraining.top().second]<<"\t"<<PQTraining.top().first<<"\t"<<0<<"\t";
		PQTraining.pop();
	}
	else{
		// this is a false Positive which needs to be recognized
		++falsePositives;
		outputFile <<  databaseIDs[PQDatabase.top().second]<<"\t"<<PQDatabase.top().first<<"\t"<<1<<"\t";
		PQDatabase.pop();
	}
	outputFile<<falsePositives<<"\t"<<falseNegatives<<"\n";
}
double sensitivity = ((double)lengthTrainingSet-(double)falseNegatives)/(double)lengthTrainingSet*(double)100;
unsigned pos=falsePositives;

// now we have estimated the equivalence number criterion
// all values left in the priority queue PQTraining are false negatives as they
// are not detected, so put them in falseNeg
while(!PQTraining.empty())
{
		if(PQDatabase.top().first <= PQTraining.top().first){
		--falseNegatives;
		//ici pour changer la taille des ids
		outputFile <<  trainingIDs[PQTraining.top().second]<<"\t"<<PQTraining.top().first<<"\t";
		outputFile <<0<<"\t"<<falsePositives<<"\t"<<falseNegatives<<"\n";
		PQTraining.pop();
	}
	else{
		++falsePositives;
		PQDatabase.pop();
	}

}
outputFile <<"\tSensitivity\tFalsePositives\n";
outputFile << "result:\t"<<sensitivity<<"\t"<<pos;
}

void printEquivalenceNumberCriterion(String<char> &outputFile,
		String<pair<double,unsigned> > &falsePos,
		String<pair<double,unsigned> > &falseNeg,
		String<String<char> > databaseIDs,
		String<String<char> > trainingIDs,
		unsigned lengthDatabase,
		unsigned lengthTrainingSet)
{
	ofstream outFile(toCString(outputFile),ios_base::app);
	if(!outFile){
		std::cerr <<outputFile<<" could not be opened!";
		exit(-1);
	}
	outFile << "Classification Result - Equivalence Number Criterion"<<endl;
	outFile << "Sensitivy:"<<(double)(lengthTrainingSet-length(falsePos))/(double)lengthTrainingSet*(double)100<<"\t";
	outFile << "Sequences total:"<<lengthDatabase+lengthTrainingSet<<" positive examples:"<<lengthTrainingSet<<endl;
	outFile <<" False Positives("<<length(falsePos)<<"):\n";
	for(unsigned i=0;i<length(falsePos);++i)
	{
		outFile <<"Norm. Likelihood:"<< falsePos[i].first <<"\tID:"<<databaseIDs[falsePos[i].second]<<endl;
	}
	outFile <<" False Negatives("<<length(falseNeg)<<"):\n";
	for(unsigned i=0;i<length(falseNeg);++i)
	{
		outFile <<"Norm. Likelihood:"<< falseNeg[i].first <<"\tID:"<<trainingIDs[falseNeg[i].second]<<endl;
	}
}
	
	
	
	
	/*Attention ce qui suit est un ajout!!!!!!!!!!!!!*/
	
	template<typename TAlphabet, typename TCargo,typename TVLMMSpec>
	inline void
	equivalenceNumberCriterion2(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM <TVLMMSpec > > > > & vlmm,
							   String<String<TAlphabet> > &database,
							   String<String<char> > &databaseIDs,
							   String<String<TAlphabet> > &trainingSet,
							   String<String<char> > &trainingIDs,
							   String<char> &output,
							   bool nl,
							   bool wl,
							   unsigned windowSize,double &sens)
	{
		// build two priority queues one for the likelihoods/scores one the trainingIndices one for the rest of the database
		//double sens;
		
		if(nl){
			fstream outFile;
			std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
			estimateNormalizedLikelihoodOnSequences(vlmm,PQDatabase,database);
			estimateNormalizedLikelihoodOnSequences(vlmm,PQTraining,trainingSet);
			openAndAppend(output,".nl.enc",outFile);
			sens=printEquivalenceNumberCriterion2(outFile,PQTraining,PQDatabase,databaseIDs,trainingIDs,length(database),length(trainingSet));
			//outFile.close();
			std::cout << sens <<endl;
			//std::cout << "wrote classification with normalized likelihood for sequences into: "<<output<<".nl.enc"<<endl;
		}
		if(wl){
			fstream outFile;
			std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
			estimateLikelihoodWindowOnSequences(vlmm,PQTraining,trainingSet,windowSize);
			estimateLikelihoodWindowOnSequences(vlmm,PQDatabase,database,windowSize);
			openAndAppend(output,".wl.enc",outFile);
			sens=printEquivalenceNumberCriterion2(outFile,PQTraining,PQDatabase,databaseIDs,trainingIDs,length(database),length(trainingSet));
			//outFile.close();
			std::cout << sens <<endl;
			//std::cout << "wrote classification with window likelihood for sequences into: "<<output<<".wl.enc"<<endl;
		}
		
	}
	
	
	
	// use the iso point method for classification
	// train only on the indices provided and do the classification on all sequences
	template<typename TAlphabet, typename TCargo>
	inline void
	equivalenceNumberCriterion2(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM <ScoreTree > > > > & scoreTree,
							   String<String<TAlphabet> > &database,
							   String<String<char> > &databaseIDs,
							   String<String<TAlphabet> > &trainingSet,
							   String<String<char> > &trainingIDs,
							   String<char> &output,
							   bool ns,
							   bool ls,
							   unsigned maxDepth,double &sens)
	{
		// build two priority queues one for the likelihoods/scores one the trainingIndices one for the rest of the database
		
		fstream outFile;
		if(ns){
			std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
			estimateNormalizedScoreOnSequences(scoreTree,PQTraining,trainingSet);
			estimateNormalizedScoreOnSequences(scoreTree,PQDatabase,database);
			openAndAppend(output,".ns.enc",outFile);
			sens=printEquivalenceNumberCriterion2(outFile,PQTraining,PQDatabase,databaseIDs,trainingIDs,length(database),length(trainingSet));
			cout << sens<<endl;
			outFile.close();
			//std::cout << "wrote classification with normalized score for sequences into: "<<output<<".ns.enc"<<endl;
		}
		if(ls){
			cout << "local Score"<<endl;
			std::priority_queue<pair<double,unsigned> > PQDatabase,PQTraining;
			std::priority_queue<pair<double,unsigned> > PQDatabaseC,PQTrainingC;
			estimateLocalScoreOnSequences(scoreTree,PQDatabase,database,maxDepth);
			estimateLocalScoreOnSequences(scoreTree,PQTraining,trainingSet,maxDepth);
			PQDatabaseC=PQDatabase;
			PQTrainingC=PQTraining;
			fstream fi;
			fi.open("/Users/gregoirelejay/Development/seqan-trunk/projects/library/apps/vlmm/FFHistDat.txt",ios_base::out|ios_base::app);
			while (!PQDatabaseC.empty()) {
				fi << PQDatabaseC.top().first <<endl;
				PQDatabaseC.pop();
			}
			fstream fc;
			fc.open("/Users/gregoirelejay/Development/seqan-trunk/projects/library/apps/vlmm/FFHistTra.txt",ios_base::out|ios_base::app);
			while (!PQTrainingC.empty()) {
				fc << PQTrainingC.top().first <<endl;
				PQTrainingC.pop();
			}
			
			openAndAppend(output,".ls.enc",outFile);
						sens=printEquivalenceNumberCriterion2(outFile,PQTraining,PQDatabase,databaseIDs,trainingIDs,length(database),length(trainingSet));
						outFile.close();
			//std::cout << "wrote classification with local score for sequences into: "<<output<<".ls.enc"<<endl;
		}
		
	}
	//Find the iso-point:
	double printEquivalenceNumberCriterion2(fstream &outputFile,
										 std::priority_queue<pair<double,unsigned> > &PQTraining,
										 std::priority_queue<pair<double,unsigned> > &PQDatabase,
										 String<String<char> > databaseIDs,
										 String<String<char> > trainingIDs,
										 unsigned /*lengthDatabase*/,
										 unsigned lengthTrainingSet)
	{	
		outputFile.precision(4);
		outputFile << "ID\tScore\tFP\t#FP\t#FN\n";
		
		// estimate the equivalence number criterion
		unsigned falsePositives=0,falseNegatives=lengthTrainingSet;
		// output statistics about it
		while(falsePositives != falseNegatives){
			if(PQDatabase.top().first <= PQTraining.top().first){
				--falseNegatives;
				outputFile <<  trainingIDs[PQTraining.top().second]<<"\t"<<PQTraining.top().first<<"\t"<<0<<"\t";
				PQTraining.pop();
			}
			else{
				// this is a false Positive which needs to be recognized
				++falsePositives;
				outputFile <<  databaseIDs[PQDatabase.top().second]<<"\t"<<PQDatabase.top().first<<"\t"<<1<<"\t";
				PQDatabase.pop();
			}
			outputFile<<falsePositives<<"\t"<<falseNegatives<<"\n";
		}
		cout << (double)lengthTrainingSet<<endl;
		cout << (double) falseNegatives<<endl;
		double sensitivity = ((double)lengthTrainingSet-(double)falseNegatives)/(double)lengthTrainingSet*(double)100;
		unsigned pos=falsePositives;
		
		// now we have estimated the equivalence number criterion
		// all values left in the priority queue PQTraining are false negatives as they
		// are not detected, so put them in falseNeg
		while(!PQTraining.empty())
		{
			if(PQDatabase.top().first <= PQTraining.top().first){
				--falseNegatives;
				outputFile <<  trainingIDs[PQTraining.top().second]<<"\t"<<PQTraining.top().first<<"\t";
				outputFile <<0<<"\t"<<falsePositives<<"\t"<<falseNegatives<<"\n";
				PQTraining.pop();
			}
			else{
				++falsePositives;
				PQDatabase.pop();
			}
			
		}
		outputFile <<"\tSensitivity\tFalsePositives\n";
		outputFile << "result:\t"<<sensitivity<<"\t"<<pos;
		return(sensitivity);
	}
	
	
} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
